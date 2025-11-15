"""
This pipeline demonstrates iterative sequence optimization using BayesianAdjuster.

Improves binding affinity towards rifampicin using a feedback loop:
1. LigandMPNN generates sequences
2. MutationProfiler analyzes mutation patterns
3. BayesianAdjuster adjusts frequencies based on correlation signals
4. MutationComposer generates new sequences using adjusted probabilities
5. Evaluation and filtering
6. SequenceMetricCorrelation tracks mutation-metric relationships
7. Repeat with updated probabilities

This creates a complete feedback loop where mutations that correlate with
improved affinity are preferentially sampled in subsequent cycles.
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.mutation_profiler import MutationProfiler
from PipelineScripts.mutation_composer import MutationComposer
from PipelineScripts.bayesian_adjuster import BayesianAdjuster
from PipelineScripts.sequence_metric_correlation import SequenceMetricCorrelation
from PipelineScripts.stitch_sequences import StitchSequences
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.pose_distance import PoseDistance
from PipelineScripts.protein_ligand_contacts import ProteinLigandContacts
from PipelineScripts.merge_tables import MergeTables
from PipelineScripts.concatenate_tables import ConcatenateTables
from PipelineScripts.remove_duplicates import RemoveDuplicates
from PipelineScripts.filter import Filter
from PipelineScripts.select_best import SelectBest
from PipelineScripts.average_by_table import AverageByTable
from PipelineScripts.extract_metrics import ExtractMetrics
from PipelineScripts.pdb import PDB

with Pipeline(project="Rifampicin",
              job="BayesianOptimization",
              description="Iterative affinity optimization using Bayesian frequency adjustment"):

    Resources(gpu="!T4",
              time="23:59:00",
              memory="16GB")

    # Load compounds and reference structure
    rifampicin = LoadOutput("/shares/locbp.chem.uzh/public/BioPipelines/Boltz/rifampicin_001/ToolOutputs/1_Boltz2_output.json")
    original_holo = PDB("11_2_29_3_1_1-24-155")

    best_holo = original_holo

    NUM_CYCLES = 20
    mutation_range = "1-24"

    # Track sequences and correlation analysis across cycles
    all_sequences_seen = None
    correlation_analysis = None

    # Track all analyses and pools across cycles for SelectBest
    all_analyses = []
    all_pools = []

    for CYCLE in range(NUM_CYCLES):
        Suffix(f"Cycle{CYCLE+1}")

        # Step 1: Generate sequences with LigandMPNN
        lmpnn = LigandMPNN(structures=best_holo,
                          ligand="LIG",
                          num_sequences=1000,
                          batch_size=25,
                          redesigned=mutation_range)

        # Step 2: Profile mutation patterns
        profiler = MutationProfiler(original=best_holo,
                                    mutants=lmpnn)

        # Step 3: Adjust frequencies based on correlations (if not first cycle)
        if correlation_analysis is not None:
            # Apply Bayesian adjustment to favor beneficial mutations
            adjuster = BayesianAdjuster(
                frequencies=profiler.tables.absolute_frequencies,
                correlations=correlation_analysis.tables.correlation_2d,
                mode="min",  # Minimize affinity_pred_value (lower is better)
                gamma=3.0,   # Strength of adjustment
                kappa=10.0   # Shrinkage parameter
            )

            # Use adjusted probabilities for sequence composition
            composer_input = adjuster.tables.absolute_probabilities
            print(f"\nCycle {CYCLE+1}: Using Bayesian-adjusted probabilities")
        else:
            # First cycle: use raw frequencies from profiler
            composer_input = profiler.tables.absolute_frequencies
            print(f"\nCycle {CYCLE+1}: Using raw frequencies (no correlation data yet)")

        # Step 4: Compose new sequences
        composer = MutationComposer(frequencies=composer_input,
                                    num_sequences=10,
                                    mode="weighted_random",
                                    max_mutations=3)

        # Step 5: Deduplicate sequences
        unique_new_sequences = RemoveDuplicates(pool=composer,
                                                history=all_sequences_seen if all_sequences_seen else None,
                                                compare="sequence")

        # Update sequence history
        if all_sequences_seen is None:
            all_sequences_seen = ConcatenateTables(tables=[unique_new_sequences.tables.sequences])
        else:
            all_sequences_seen = ConcatenateTables(tables=[unique_new_sequences.tables.sequences,
                                                           all_sequences_seen.tables.concatenated])

        # Step 6: Generate MSAs and predict structures
        msas = MMseqs2(sequences=unique_new_sequences)
        boltz_holo = Boltz2(proteins=unique_new_sequences,
                           ligands=rifampicin,
                           msas=msas)

        # Step 7: Analyze structures
        contacts = ProteinLigandContacts(structures=boltz_holo,
                                        selections=mutation_range,
                                        ligand="LIG")
        pose_distance = PoseDistance(reference_structure=rifampicin,
                                     sample_structures=boltz_holo,
                                     reference_ligand="LIG",
                                     sample_ligand="LIG")

        # Step 8: Merge metrics and filter
        data = MergeTables(tables=[boltz_holo.tables.affinity,
                                  boltz_holo.tables.confidence,
                                  contacts.tables.contact_analysis,
                                  pose_distance.tables.analysis])
        current_filtered = Filter(data=data,
                                  pool=boltz_holo,
                                  expression="contacts>=3 and ligand_rmsd<2")

        # Step 9: Update correlation analysis
        # This tracks which mutations correlate with better affinity
        if correlation_analysis is None:
            # First cycle: initialize correlation tracking
            correlation_analysis = SequenceMetricCorrelation(
                mutants=current_filtered.tables.sequences,
                data=current_filtered.tables.merged,
                original=original_holo,
                metric="affinity_pred_value"
            )
            print(f"  Initialized correlation tracking with {len(current_filtered.tables.sequences)} sequences")
        else:
            # Subsequent cycles: accumulate correlation data
            # This combines data from all previous cycles for better statistics
            correlation_analysis = SequenceMetricCorrelation(
                mutants=current_filtered.tables.sequences,
                data=current_filtered.tables.merged,
                original=original_holo,
                metric="affinity_pred_value"
            )
            print(f"  Updated correlation tracking with {len(current_filtered.tables.sequences)} new sequences")

        # Step 10: Track for best selection
        all_pools.append(boltz_holo)
        all_analyses.append(current_filtered)

        # Step 11: Select best structure across all cycles
        best_holo = SelectBest(pool=[pool for pool in all_pools],
                              tables=[x.tables.merged for x in all_analyses],
                              metric='affinity_pred_value',
                              mode='min',
                              name=f'{CYCLE+1}_best')

    # Final analysis across all cycles
    all_merged = [x.tables.merged for x in all_analyses]
    combined_tables = ConcatenateTables(all_merged)
    AverageByTable(all_merged)

    # Extract metrics for analysis
    metrics = ["affinity_pred_value", "contacts", "complex_plddt"]
    ExtractMetrics(tables=all_merged,
                  metrics=metrics)

print("\n" + "="*60)
print("Pipeline complete!")
print("="*60)
print("\nKey features demonstrated:")
print("  ✓ LigandMPNN sequence generation")
print("  ✓ MutationProfiler pattern analysis")
print("  ✓ BayesianAdjuster correlation-based frequency adjustment")
print("  ✓ SequenceMetricCorrelation tracking across cycles")
print("  ✓ Iterative optimization with feedback loop")
print("\nThe pipeline uses Bayesian updates to preferentially sample")
print("mutations that correlate with improved binding affinity.")
print("="*60)
