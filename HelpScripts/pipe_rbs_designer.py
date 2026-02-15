#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime script for RBS design using the Salis thermodynamic model.

Designs synthetic ribosome binding sites (RBS) via simulated annealing to match
a target translation initiation rate (TIR). Uses ViennaRNA for RNA free energy
calculations.

Reference: Salis, Mirsky & Voigt, Nature Biotechnology 2009.
           doi:10.1038/nbt.1568

The thermodynamic model (Equation 2 in the paper):

    dG_tot = dG_mRNA:rRNA + dG_start + dG_spacing - dG_standby - dG_mRNA

where:
    dG_mRNA:rRNA  : Hybridization energy between mRNA SD site and 16S rRNA
                    anti-SD (3'-AUUCCUCCACUAG, last 9 nt: ACCUCCUUA).
                    Computed using NuPACK/ViennaRNA duplex folding, scanning
                    all suboptimal configurations in the upstream window.
    dG_start      : Start codon identity: AUG = -1.194, GUG = -0.075 kcal/mol
    dG_spacing    : Penalty for non-optimal aligned spacing s between the 3'
                    end of the 16S rRNA binding site and the start codon.
                    s > 5 nt: quadratic  c1*(s-5)^2 + c2*(s-5)   (Eq. 7 in SI)
                    s < 5 nt: sigmoidal  c1/[1+exp(c2*(s-s_opt+2))]^3  (Eq. 8)
    dG_standby    : Energy to unfold the 4-nt standby site upstream of the
                    16S rRNA binding site (ref 31: de Smit & van Duin, 2003).
    dG_mRNA       : MFE of the mRNA sub-sequence S2 (n_start-35 to n_start+35).

The translation initiation rate is proportional:  r = K * exp(-beta * dG_tot)
where beta = 0.45 mol/kcal (experimentally determined) and K = 2500.

Usage: python pipe_rbs_designer.py --config <config.json>
"""

import argparse
import json
import math
import random
import sys
import time

import numpy as np
import pandas as pd

try:
    import RNA
except ImportError:
    print("ERROR: ViennaRNA Python bindings not found.")
    print("Install via conda: conda install -c bioconda viennarna")
    sys.exit(1)


# =============================================================================
# Constants from Salis et al. 2009
# =============================================================================

# Last 9 nt of E. coli 16S rRNA 3' tail (anti-Shine-Dalgarno)
# Full 3' tail is 3'-AUUCCUCCACUAG-5', the 9 nt used for hybridization:
RRNA_16S = "ACCUCCUUA"

# Model parameters (from paper main text and SI)
OPTIMAL_SPACING = 5       # optimal aligned spacing (nt) between SD and start codon
BETA = 0.45               # apparent Boltzmann factor (mol/kcal), Fig 2b
K = 2500.0                # proportionality constant for TIR, SI Fig S1

# Subsequence cutoff (nt before and after start codon), SI Section 4, Fig S3
CUTOFF = 35

# Start codon free energies (kcal/mol) — paper p.947
START_CODON_ENERGIES = {
    "AUG": -1.194,
    "GUG": -0.075,
    "UUG":  0.0,
}

# Spacing penalty coefficients — SI Equations 7 and 8
# Quadratic (s > 5): dG_spacing = c1*(s-s_opt)^2 + c2*(s-s_opt)
SPACING_C1_QUAD = 0.048     # kcal/mol/nt^2
SPACING_C2_QUAD = 0.24      # kcal/mol/nt
# Sigmoidal (s < 5): dG_spacing = c1 / [1 + exp(c2*(s - s_opt + 2))]^3
SPACING_C1_SIG = 12.2       # kcal/mol
SPACING_C2_SIG = 2.5        # nt^-1

# Standby site length (nt upstream of the 16S rRNA binding site)
STANDBY_SITE_LEN = 4        # paper p.947

# TIR presets (au, from paper Fig 2)
TIR_PRESETS = {
    "low": 100,
    "medium": 1000,
    "high": 10000,
    "maximum": 100000,
}

# Simulated annealing parameters — Online Methods
SA_MAX_ITERATIONS = 10000
SA_CONVERGENCE_THRESHOLD = 0.25   # kcal/mol, paper main text p.948
SA_RBS_MIN_LEN = 5
SA_RBS_MAX_LEN = 35
SA_MAX_UNFOLD_ENERGY = 6.0        # kcal/mol, Online Methods constraint 1
SA_MIN_BASE_PAIR_PROB = 6e-5      # Online Methods constraint 2
SA_BP_EXPONENT = -1.44            # growth model exponent, ref 34

# Adaptive temperature: paper says "T_SA is continually adjusted to maintain
# a 5-20% acceptance rate." We use initial temp with geometric cooling,
# combined with reheating when acceptance rate drops below target.
SA_INITIAL_TEMP = 10.0
SA_COOLING_RATE = 0.99
SA_MIN_TEMP = 0.001
SA_TARGET_ACCEPT_LOW = 0.05
SA_TARGET_ACCEPT_HIGH = 0.20
SA_ACCEPT_WINDOW = 50             # track acceptance over last N trials

BASES_DNA = "ATCG"


# =============================================================================
# Sequence conversion utilities
# =============================================================================

def dna_to_rna(seq):
    """Convert DNA sequence to RNA (T -> U)."""
    return seq.upper().replace("T", "U")


def rna_to_dna(seq):
    """Convert RNA sequence to DNA (U -> T)."""
    return seq.upper().replace("U", "T")


# =============================================================================
# ViennaRNA wrapper functions
# =============================================================================

def rna_mfe(sequence):
    """
    Compute minimum free energy of an RNA sequence using ViennaRNA.

    Uses the Mfold 3.0 energy parameters (Turner 2004) as bundled with ViennaRNA,
    consistent with the paper's use of NuPACK 'mfe' and Mfold 3.0.

    Args:
        sequence: RNA sequence string (AUCG).

    Returns:
        Minimum free energy in kcal/mol.
    """
    fc = RNA.fold_compound(sequence)
    _, mfe = fc.mfe()
    return mfe


def rna_duplex(seq1, seq2):
    """
    Compute hybridization (duplex) energy between two RNA sequences.

    Equivalent to the NuPACK 'subopt' calculation in the paper,
    returning the minimum free energy of intermolecular base pairing.

    Args:
        seq1: First RNA sequence.
        seq2: Second RNA sequence.

    Returns:
        Duplex free energy in kcal/mol.
    """
    result = RNA.duplexfold(seq1, seq2)
    return result.energy


# =============================================================================
# TIR / dG conversions — Equation 1 and SI Equation S7
# =============================================================================

def tir_to_dg(tir):
    """Convert translation initiation rate to total free energy.

    From r = K * exp(-beta * dG_tot), solving for dG_tot:
        dG_tot = -ln(r / K) / beta
    """
    return -math.log(tir / K) / BETA


def dg_to_tir(dg):
    """Convert total free energy to translation initiation rate.

    Equation 1: r proportional to exp(-beta * dG_tot)
    SI Equation S7: r_i = K * exp(-beta * dG_i)
    """
    return K * math.exp(-BETA * dg)


# =============================================================================
# Spacing penalty — SI Equations 7 and 8
# =============================================================================

def calc_spacing_penalty(spacing):
    """
    Calculate the spacing penalty dG_spacing.

    From SI Section 3 (Eqs. 7-8):
      s > s_opt (stretched): dG = c1*(s - s_opt)^2 + c2*(s - s_opt)
      s < s_opt (compressed): dG = c1 / [1 + exp(c2*(s - s_opt + 2))]^3

    where s_opt = 5, c1/c2 differ between the two regimes.

    Args:
        spacing: Aligned spacing s (nt) between 3' end of 16S rRNA
                 binding site and first nucleotide of start codon.

    Returns:
        Spacing penalty in kcal/mol (>= 0).
    """
    s = spacing
    if s >= OPTIMAL_SPACING:
        ds = s - OPTIMAL_SPACING
        return SPACING_C1_QUAD * ds * ds + SPACING_C2_QUAD * ds
    else:
        return SPACING_C1_SIG / (
            (1.0 + math.exp(SPACING_C2_SIG * (s - OPTIMAL_SPACING + 2))) ** 3
        )


# =============================================================================
# Thermodynamic model — Paper Equation 2, Online Methods
# =============================================================================

def calc_dg_total(mrna_rna, start_pos):
    """
    Compute the total free energy of translation initiation (dG_total).

    Implements the 5-term model from Salis et al. 2009 Equation 2:
        dG_tot = dG_mRNA:rRNA + dG_start + dG_spacing - dG_standby - dG_mRNA

    The mRNA sub-sequences (Online Methods):
        S1 = mrna[max(1, n_start - 35) : n_start]     (upstream only, for SD scan)
        S2 = mrna[max(1, n_start - 35) : n_start + 35] (for dG_mRNA folding)

    Spacing is computed as: s = n_start - n1 - n2
    where n1, n2 are the mRNA and rRNA positions of the farthest 3' base pair
    in the 16S rRNA binding site.

    dG_standby is computed by constraining the 4-nt standby site (immediately
    upstream of the 16S binding site) to be single-stranded and comparing
    the folding energy of S2 with and without this constraint.

    Args:
        mrna_rna: Full mRNA sequence in RNA (5'UTR + RBS + CDS region).
        start_pos: 0-based index of the first nucleotide of the start codon.

    Returns:
        Dictionary with all free energy terms and total dG.
    """
    # --- dG_start: start codon identity (paper p.947) ---
    start_codon = mrna_rna[start_pos:start_pos + 3]
    dg_start = START_CODON_ENERGIES.get(start_codon, 0.0)

    # --- Define sub-sequences S1 and S2 (Online Methods) ---
    s1_start = max(0, start_pos - CUTOFF)
    s2_start = max(0, start_pos - CUTOFF)
    s2_end = min(len(mrna_rna), start_pos + CUTOFF)

    # S1: upstream of start codon only (for SD scanning)
    upstream_seq = mrna_rna[s1_start:start_pos]

    # S2: the full folding window around start codon (for dG_mRNA)
    s2_seq = mrna_rna[s2_start:s2_end]

    # --- dG_mRNA:rRNA: SD / anti-SD hybridization (paper p.947) ---
    # "All possible hybridizations between the mRNA and 16S rRNA are
    # considered to find the highest affinity 16S rRNA binding site."
    # We use ViennaRNA duplexfold which considers all possible alignments.
    # The binding site that minimizes the sum dG_mRNA:rRNA + dG_spacing
    # is selected (Online Methods).
    #
    # We scan all possible binding positions in the upstream window.
    # For each position, we compute the duplex energy of the local
    # subsequence against the full anti-SD, then add the spacing penalty.

    anti_sd = RRNA_16S
    best_combined = 0.0  # dG_mRNA:rRNA + dG_spacing (want to minimize)
    best_duplex_energy = 0.0
    best_spacing = OPTIMAL_SPACING
    best_sd_end = start_pos  # position after last nt of SD site in mRNA

    if len(upstream_seq) >= 4:
        # Try all sub-windows of the upstream region against the anti-SD
        for win_len in range(4, min(len(anti_sd), len(upstream_seq)) + 1):
            for offset in range(len(upstream_seq) - win_len + 1):
                subseq = upstream_seq[offset:offset + win_len]
                energy = rna_duplex(subseq, anti_sd)

                if energy >= 0.0:
                    continue  # no favorable binding

                # The 3' end of this SD binding site in the mRNA
                sd_end_pos = s1_start + offset + win_len
                # Aligned spacing: distance from SD 3' end to start codon
                spacing = start_pos - sd_end_pos

                dg_sp = calc_spacing_penalty(spacing)
                combined = energy + dg_sp

                if combined < best_combined:
                    best_combined = combined
                    best_duplex_energy = energy
                    best_spacing = spacing
                    best_sd_end = sd_end_pos

    dg_mrna_rrna = best_duplex_energy
    spacing = best_spacing
    dg_spacing = calc_spacing_penalty(spacing)

    # --- dG_mRNA: MFE of sub-sequence S2 (Online Methods) ---
    # "the mfe configuration of sequence S2 is calculated and its free
    # energy is designated dG_mRNA"
    dg_mrna = rna_mfe(s2_seq) if len(s2_seq) >= 4 else 0.0

    # --- dG_standby: energy to unfold the standby site (paper p.947) ---
    # "dG_standby is the work required to unfold any secondary structures
    # sequestering the standby site (dG_standby < 0) after the 30S complex
    # assembly. We define the standby site as the four nucleotides upstream
    # of the 16S rRNA binding site."
    #
    # Online Methods: "The energy required to unfold the standby site is
    # determined by calculating the mfe of sequence S2 with and without
    # the standby site constrained to be single-stranded. The difference
    # between these mfes is designated dG_standby."
    #
    # To calculate the constrained mfe: the standby site (4 nt immediately
    # upstream of the 16S binding site) must be unpaired. We split S2 into
    # sub-sequences around the standby site and sum their mfes.
    dg_standby = 0.0
    standby_end = best_sd_end               # standby ends where SD starts
    standby_start = standby_end - STANDBY_SITE_LEN

    if standby_start >= s2_start and standby_end <= s2_end:
        # Positions relative to S2
        sb_rel_start = standby_start - s2_start
        sb_rel_end = standby_end - s2_start

        # Constrained folding: force standby site to be single-stranded
        # by splitting S2 into the part before standby and the part after
        # (including standby nucleotides as unpaired), then summing mfes.
        part_before = s2_seq[:sb_rel_start]
        part_after = s2_seq[sb_rel_end:]

        mfe_constrained = 0.0
        if len(part_before) >= 4:
            mfe_constrained += rna_mfe(part_before)
        if len(part_after) >= 4:
            mfe_constrained += rna_mfe(part_after)

        # dG_standby = dG_mRNA(unconstrained) - dG_mRNA(constrained)
        # This is negative when the standby site participates in structure
        dg_standby = dg_mrna - mfe_constrained

    # --- Total: Equation 2 ---
    dg_total = dg_mrna_rrna + dg_start + dg_spacing - dg_standby - dg_mrna

    return {
        "dg_total": dg_total,
        "dg_mrna_rrna": dg_mrna_rrna,
        "dg_start": dg_start,
        "dg_spacing": dg_spacing,
        "dg_mrna": dg_mrna,
        "dg_standby": dg_standby,
        "spacing": spacing,
        "start_codon": start_codon,
    }


# =============================================================================
# Sequence constraint checks — Online Methods
# =============================================================================

def has_internal_start(rbs_rna):
    """Check if an RBS RNA sequence contains an internal start codon.

    Online Methods constraint 3: "the creation of new AUG or GUG start
    codons within the RBS sequence is disallowed."
    """
    for codon in ("AUG", "GUG"):
        if codon in rbs_rna:
            return True
    return False


def check_long_range_pairs(s1_rna):
    """Check for long-range base pair interactions (Online Methods constraint 2).

    "According to a growth model for random RNA sequences, the equilibrium
    probability P of nucleotides i and j forming a base pair in solution is
    proportional to P = |i-j|^(-1.44). For each base pair in sequence S1,
    we calculate P. If the minimum P is < 6 x 10^-5, then the sequence
    is rejected."

    Args:
        s1_rna: The upstream mRNA sub-sequence S1.

    Returns:
        True if the sequence passes (no problematic long-range pairs),
        False if it should be rejected.
    """
    n = len(s1_rna)
    if n < 10:
        return True

    # Use ViennaRNA to get the MFE structure and check base pair distances
    fc = RNA.fold_compound(s1_rna)
    structure, _ = fc.mfe()

    # Parse dot-bracket to find base pairs
    stack = []
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if stack:
                j = stack.pop()
                dist = abs(i - j)
                if dist > 0:
                    prob = dist ** SA_BP_EXPONENT
                    if prob < SA_MIN_BASE_PAIR_PROB:
                        return False
    return True


# =============================================================================
# RBS sequence utilities
# =============================================================================

def random_rbs(length=15):
    """Generate a random RBS DNA sequence without internal start codons."""
    for _ in range(1000):
        seq = ''.join(random.choice(BASES_DNA) for _ in range(length))
        rna = dna_to_rna(seq)
        if not has_internal_start(rna):
            return seq
    return "A" * length


def mutate_rbs(rbs_dna):
    """
    Apply a single random mutation to an RBS DNA sequence.

    Online Methods: "the simulated annealing optimization algorithm randomly
    deletes, inserts or replaces a nucleotide in the RBS sequence."

    Returns:
        Mutated RBS DNA sequence.
    """
    seq = list(rbs_dna)
    r = random.random()

    if r < 0.70:
        # Point mutation (replace)
        pos = random.randint(0, len(seq) - 1)
        old = seq[pos]
        new_bases = [b for b in BASES_DNA if b != old]
        seq[pos] = random.choice(new_bases)
    elif r < 0.85:
        # Insertion
        if len(seq) < SA_RBS_MAX_LEN:
            pos = random.randint(0, len(seq))
            seq.insert(pos, random.choice(BASES_DNA))
    else:
        # Deletion
        if len(seq) > SA_RBS_MIN_LEN:
            pos = random.randint(0, len(seq) - 1)
            seq.pop(pos)

    return ''.join(seq)


# =============================================================================
# Simulated annealing RBS designer — Online Methods
# =============================================================================

def design_rbs(cds_dna, target_dg, pre_sequence=""):
    """
    Design an RBS sequence to achieve a target dG_total using simulated annealing.

    From Online Methods:
    "An initial RBS sequence is randomly generated and inserted in between a
    presequence and protein coding sequence to create a sequence S. The dG_tot
    of the sequence S is calculated and the objective function
    O_old = |dG_tot - dG_target| is evaluated. In an iterative procedure, the
    simulated annealing optimization algorithm randomly deletes, inserts or
    replaces a nucleotide in the RBS sequence. The dG_tot and objective function
    O_new are then recalculated. If the dG_tot calculation of S invalidates the
    sequence constraints, then the mutation is immediately rejected. Otherwise,
    the mutation is accepted with probability max(1, exp((O_old - O_new)/T_SA)),
    where T_SA is the simulated annealing temperature. The T_SA is continually
    adjusted to maintain a 5-20% acceptance rate."

    "The procedure continues until the synthetic sequence has a predicted dG_tot
    to within 0.25 kcal/mol of the target."

    Args:
        cds_dna: Coding DNA sequence (starts with ATG/GTG).
        target_dg: Target total free energy (kcal/mol).
        pre_sequence: Optional fixed 5'UTR DNA prepended before the RBS.

    Returns:
        Dictionary with designed RBS and thermodynamic results.
    """
    t_start = time.time()

    cds_rna = dna_to_rna(cds_dna)
    pre_rna = dna_to_rna(pre_sequence) if pre_sequence else ""

    # Initialize with random RBS
    best_rbs_dna = random_rbs(15)
    best_rbs_rna = dna_to_rna(best_rbs_dna)

    full_mrna = pre_rna + best_rbs_rna + cds_rna
    start_pos = len(pre_rna) + len(best_rbs_rna)

    best_result = calc_dg_total(full_mrna, start_pos)
    best_cost = abs(best_result["dg_total"] - target_dg)

    current_rbs_dna = best_rbs_dna
    current_cost = best_cost

    temperature = SA_INITIAL_TEMP

    # Track acceptance rate for adaptive temperature (paper: 5-20%)
    recent_accepted = 0
    recent_total = 0
    total_rejected_constraint = 0

    # Progress logging interval
    log_interval = max(SA_MAX_ITERATIONS // 10, 100)

    print(f"\n        Initial: dG={best_result['dg_total']:.3f} (target={target_dg:.3f}), "
          f"|error|={best_cost:.3f} kcal/mol", flush=True)

    for iteration in range(SA_MAX_ITERATIONS):
        # Check convergence
        if best_cost < SA_CONVERGENCE_THRESHOLD:
            elapsed = time.time() - t_start
            print(f"        Converged at iteration {iteration}: |error|={best_cost:.3f} kcal/mol "
                  f"({elapsed:.1f}s)", flush=True)
            break

        # Progress logging
        if iteration > 0 and iteration % log_interval == 0:
            elapsed = time.time() - t_start
            print(f"        Iter {iteration}/{SA_MAX_ITERATIONS}: "
                  f"best dG={best_result['dg_total']:.3f}, |error|={best_cost:.3f}, "
                  f"T={temperature:.4f}, rejected={total_rejected_constraint} "
                  f"({elapsed:.1f}s)", flush=True)

        # Mutate
        candidate_rbs_dna = mutate_rbs(current_rbs_dna)
        candidate_rbs_rna = dna_to_rna(candidate_rbs_dna)

        # Constraint 3: no internal start codons
        if has_internal_start(candidate_rbs_rna):
            total_rejected_constraint += 1
            continue

        # Build candidate mRNA
        candidate_mrna = pre_rna + candidate_rbs_rna + cds_rna
        candidate_start = len(pre_rna) + len(candidate_rbs_rna)

        # Constraint 1: unfolding energy of the 16S rRNA binding site
        # "The first constraint calculates the energy required to unfold
        # the 16S rRNA binding site on the mRNA sequence and rejects the
        # ones that require >6 kcal/mol to unfold."
        fold_start = max(0, candidate_start - CUTOFF)
        fold_end = min(len(candidate_mrna), candidate_start + CUTOFF)
        fold_region = candidate_mrna[fold_start:fold_end]
        if len(fold_region) >= 4:
            if abs(rna_mfe(fold_region)) > SA_MAX_UNFOLD_ENERGY:
                total_rejected_constraint += 1
                continue

        # Constraint 2: long-range base pair interactions
        upstream_s1 = candidate_mrna[max(0, candidate_start - CUTOFF):candidate_start]
        if len(upstream_s1) >= 10:
            if not check_long_range_pairs(upstream_s1):
                total_rejected_constraint += 1
                continue

        # Evaluate thermodynamic model
        candidate_result = calc_dg_total(candidate_mrna, candidate_start)
        candidate_cost = abs(candidate_result["dg_total"] - target_dg)

        # Metropolis criterion: accept with probability max(1, exp((O_old - O_new)/T_SA))
        delta = candidate_cost - current_cost
        accepted = False
        if delta <= 0:
            accepted = True
        elif temperature > 0:
            accept_prob = math.exp(-delta / temperature)
            if random.random() < accept_prob:
                accepted = True

        recent_total += 1

        if accepted:
            current_rbs_dna = candidate_rbs_dna
            current_cost = candidate_cost
            recent_accepted += 1

            if candidate_cost < best_cost:
                best_rbs_dna = candidate_rbs_dna
                best_cost = candidate_cost
                best_result = candidate_result

        # Adaptive temperature adjustment (Online Methods):
        # "T_SA is continually adjusted to maintain a 5-20% acceptance rate"
        if recent_total >= SA_ACCEPT_WINDOW:
            accept_rate = recent_accepted / recent_total
            if accept_rate < SA_TARGET_ACCEPT_LOW:
                temperature *= 1.5   # reheat to increase acceptance
            elif accept_rate > SA_TARGET_ACCEPT_HIGH:
                temperature *= SA_COOLING_RATE  # cool to decrease acceptance
            else:
                temperature *= SA_COOLING_RATE  # gentle cooling in target range
            recent_accepted = 0
            recent_total = 0

            if temperature < SA_MIN_TEMP:
                temperature = SA_MIN_TEMP
    else:
        # Loop completed without convergence
        elapsed = time.time() - t_start
        print(f"        Max iterations reached ({SA_MAX_ITERATIONS}): "
              f"|error|={best_cost:.3f} kcal/mol ({elapsed:.1f}s)", flush=True)

    # Final evaluation with best RBS
    best_rbs_rna = dna_to_rna(best_rbs_dna)
    final_mrna = pre_rna + best_rbs_rna + cds_rna
    final_start = len(pre_rna) + len(best_rbs_rna)
    final_result = calc_dg_total(final_mrna, final_start)
    elapsed = time.time() - t_start

    return {
        "rbs_dna": best_rbs_dna,
        "rbs_rna": best_rbs_rna,
        "dg_total": final_result["dg_total"],
        "tir_predicted": dg_to_tir(final_result["dg_total"]),
        "spacing": final_result["spacing"],
        "dg_mrna_rrna": final_result["dg_mrna_rrna"],
        "dg_start": final_result["dg_start"],
        "dg_spacing": final_result["dg_spacing"],
        "dg_mrna": final_result["dg_mrna"],
        "dg_standby": final_result["dg_standby"],
        "iterations": iteration + 1,
        "elapsed_s": elapsed,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="RBS design using Salis thermodynamic model")
    parser.add_argument("--config", type=str, required=True, help="Path to config JSON file")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = json.load(f)

    sequences_csv = config["sequences_csv"]
    target_tir = float(config["tir"])
    pre_sequence = config.get("pre_sequence", "")
    rbs_output = config["rbs_output"]
    info_output = config["info_output"]

    target_dg = tir_to_dg(target_tir)

    total_start = time.time()

    print("=" * 70)
    print("RBS Designer - Salis Thermodynamic Model")
    print("Reference: Salis, Mirsky & Voigt, Nat. Biotechnol. 27, 946-950 (2009)")
    print("=" * 70)
    print()
    print(f"  Input:          {sequences_csv}")
    print(f"  Target TIR:     {target_tir:.0f} au")
    print(f"  Target dG_tot:  {target_dg:.3f} kcal/mol")
    print(f"  Pre-sequence:   {pre_sequence if pre_sequence else '(none)'}")
    print(f"  Convergence:    |dG_tot - target| < {SA_CONVERGENCE_THRESHOLD} kcal/mol")
    print(f"  Max iterations: {SA_MAX_ITERATIONS} per sequence")
    print()

    sequences_df = pd.read_csv(sequences_csv)

    if "dna_sequence" in sequences_df.columns:
        seq_col = "dna_sequence"
    elif "sequence" in sequences_df.columns:
        seq_col = "sequence"
    else:
        raise ValueError("Input CSV must have 'dna_sequence' or 'sequence' column")

    if "id" not in sequences_df.columns:
        raise ValueError("Input CSV must have 'id' column")

    n_seqs = len(sequences_df)
    print(f"Designing RBS for {n_seqs} sequence{'s' if n_seqs != 1 else ''}...")
    print("-" * 70)

    results = []
    for idx, row in sequences_df.iterrows():
        seq_id = row["id"]
        cds_dna = row[seq_col].upper()

        print(f"\n  [{idx + 1}/{n_seqs}] {seq_id} ({len(cds_dna)} bp)", flush=True)

        result = design_rbs(cds_dna, target_dg, pre_sequence=pre_sequence)

        full_gene = pre_sequence + result["rbs_dna"] + cds_dna

        results.append({
            "id": seq_id,
            "dna_sequence": cds_dna,
            "rbs_sequence": result["rbs_dna"],
            "full_gene": full_gene,
            "dg_total": round(result["dg_total"], 4),
            "tir_predicted": round(result["tir_predicted"], 2),
            "target_tir": target_tir,
            "target_dg": round(target_dg, 4),
            "spacing": result["spacing"],
            "dg_mrna_rrna": round(result["dg_mrna_rrna"], 4),
            "dg_start": round(result["dg_start"], 4),
            "dg_spacing": round(result["dg_spacing"], 4),
            "dg_mrna": round(result["dg_mrna"], 4),
            "dg_standby": round(result["dg_standby"], 4),
        })

        converged = abs(result["dg_total"] - target_dg) < SA_CONVERGENCE_THRESHOLD
        status = "CONVERGED" if converged else "best found"
        print(f"        Result ({status}): RBS={result['rbs_dna']} "
              f"({len(result['rbs_dna'])} nt)", flush=True)
        print(f"          dG_tot={result['dg_total']:.3f} kcal/mol, "
              f"TIR={result['tir_predicted']:.0f} au, "
              f"spacing={result['spacing']} nt, "
              f"iters={result['iterations']}, "
              f"time={result['elapsed_s']:.1f}s", flush=True)
        print(f"          dG_mRNA:rRNA={result['dg_mrna_rrna']:.3f}, "
              f"dG_start={result['dg_start']:.3f}, "
              f"dG_spacing={result['dg_spacing']:.3f}, "
              f"dG_mRNA={result['dg_mrna']:.3f}, "
              f"dG_standby={result['dg_standby']:.3f}", flush=True)

    results_df = pd.DataFrame(results)
    results_df.to_csv(rbs_output, index=False)

    total_elapsed = time.time() - total_start

    # Summary
    print()
    print("-" * 70)
    print("Summary")
    print("-" * 70)
    n_converged = sum(1 for r in results if abs(r["dg_total"] - target_dg) < SA_CONVERGENCE_THRESHOLD)
    dg_errors = [abs(r["dg_total"] - target_dg) for r in results]
    print(f"  Sequences:  {n_seqs}")
    print(f"  Converged:  {n_converged}/{n_seqs}")
    if dg_errors:
        print(f"  |dG error|: mean={np.mean(dg_errors):.3f}, "
              f"max={np.max(dg_errors):.3f} kcal/mol")
    print(f"  Total time: {total_elapsed:.1f}s "
          f"({total_elapsed / n_seqs:.1f}s/sequence)")
    print(f"  Output:     {rbs_output}")
    print()

    info_text = f"""RBS Design Information
======================

Method: Salis Thermodynamic Model with Simulated Annealing
Reference: Salis, Mirsky & Voigt, Nat. Biotechnol. 27, 946-950 (2009)
           doi:10.1038/nbt.1568

Model (Equation 2):
  dG_tot = dG_mRNA:rRNA + dG_start + dG_spacing - dG_standby - dG_mRNA

Parameters:
  16S rRNA anti-SD  : {RRNA_16S} (last 9 nt of E. coli 16S rRNA 3' tail)
  Optimal spacing   : {OPTIMAL_SPACING} nt
  Beta              : {BETA} mol/kcal (experimentally measured, Fig. 2b)
  K                 : {K} (proportionality constant, SI Fig. S1)
  Subsequence cutoff: {CUTOFF} nt (before and after start codon, SI Section 4)

Start codon energies:
  AUG = {START_CODON_ENERGIES['AUG']} kcal/mol
  GUG = {START_CODON_ENERGIES['GUG']} kcal/mol

Spacing penalty (SI Eqs. 7-8):
  s > 5 nt: dG = {SPACING_C1_QUAD}*(s-5)^2 + {SPACING_C2_QUAD}*(s-5)
  s < 5 nt: dG = {SPACING_C1_SIG} / [1 + exp({SPACING_C2_SIG}*(s-3))]^3

Target TIR: {target_tir:.0f}
Target dG_total: {target_dg:.3f} kcal/mol
Pre-sequence: {pre_sequence if pre_sequence else '(none)'}

Simulated Annealing (Online Methods):
  Convergence: |dG_tot - dG_target| < {SA_CONVERGENCE_THRESHOLD} kcal/mol
  Max iterations: {SA_MAX_ITERATIONS}
  Adaptive temperature: 5-20% acceptance rate
  Sequence constraints:
    1. Unfolding energy of mRNA < {SA_MAX_UNFOLD_ENERGY} kcal/mol
    2. No long-range base pairs (P = |i-j|^{SA_BP_EXPONENT} < {SA_MIN_BASE_PAIR_PROB})
    3. No internal AUG/GUG start codons in RBS

Total sequences designed: {len(results_df)}
"""
    with open(info_output, "w") as f:
        f.write(info_text)
    print(f"Saved design info to: {info_output}")
    print("\nRBS design completed successfully!")


if __name__ == "__main__":
    main()
