# Example Pipelines

Complete example pipelines demonstrating BioPipelines capabilities. All examples are available in the `ExamplePipelines/` folder.

---

## Inverse Folding: Ubiquitin

**File:** `ubiquitin.py`

Retrieve a structure from PDB, design new sequences with ProteinMPNN, validate with AlphaFold, generate DNA sequences, and filter results based on confidence and structural similarity.

```python
from biopipelines.pipeline import *
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.conformational_change import ConformationalChange
from biopipelines.panda import Panda
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="Ubiquitin", job="InverseFolding"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")
    ubiquitin = PDB("4LCD", chain="E")
    sequences = ProteinMPNN(structures=ubiquitin,
                            num_sequences=10,
                            soluble_model=True)
    folded = AlphaFold(proteins=sequences)
    dna = DNAEncoder(sequences=sequences, organism="EC")
    conf_change = ConformationalChange(reference_structures=ubiquitin,
                                       target_structures=folded)
    filtered_sequences = Panda(tables=[folded.tables.confidence,
                                       conf_change.tables.changes],
                               operations=[Panda.merge(),
                                           Panda.filter("RMSD < 1.5 and plddt > 80")],
                               pool=sequences)
    dna = DNAEncoder(sequences=filtered_sequences, organism="EC")
```

---

## RFdiffusion3 + ProteinMPNN + LigandMPNN + Boltz2

**File:** `rfd3_pmpnn_ligandmpnn_boltz2.py`

Enzyme redesign using RFdiffusion3 for backbone generation, split sequence design (ProteinMPNN for positions far from the ligand, LigandMPNN for positions near the ligand), stitching, and multi-condition Boltz2 folding.

```python
from biopipelines.pipeline import *
from biopipelines.rfdiffusion3 import RFdiffusion3
from biopipelines.distance_selector import DistanceSelector
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.stitch_sequences import StitchSequences
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="Examples",
              job="RFD3-ProteinMPNN-LigandMPNN-Boltz",
              description="Redesign of a portion of adenilate kinase"):

    Resources(gpu="any", time="4:00:00", memory="16GB")

    adenylate_kinase = PDB("3BE4")
    atp = Ligand("ATP")
    amp = Ligand("AMP")
    adp = Ligand("ADP")
    ap5 = Ligand("AP5")

    # Boltz2 prediction as clean starting structure
    adenylate_kinase_boltz = Boltz2(proteins=adenylate_kinase, ligands=ap5)
    adenylate_kinase_boltz_renamed = PDB(adenylate_kinase_boltz,
                                         PDB.Rename("LIG", ":L:"))

    rfd3 = RFdiffusion3(pdb=adenylate_kinase_boltz_renamed,
                        ligand_code=':L:',
                        contig='A1-121,50-70,A170-214',
                        num_designs=3)

    # Split sequence design: far from ligand vs near ligand
    distances = DistanceSelector(structures=rfd3, ligand=":L:", distance=5,
                                  restrict_to=rfd3.tables.structures.designed)
    pmpnn = ProteinMPNN(structures=rfd3, num_sequences=2,
                        redesigned=distances.tables.selections.beyond)
    lmpnn = LigandMPNN(structures=rfd3, ligand=":L:", num_sequences=2,
                       redesigned=distances.tables.selections.within)

    sequences = StitchSequences(template=rfd3,
                                substitutions={
                                    distances.tables.selections.beyond: pmpnn,
                                    distances.tables.selections.within: lmpnn
                                })

    # Multi-condition folding
    boltz_apo = Boltz2(proteins=sequences)
    boltz_atp = Boltz2(proteins=sequences, ligands=atp, msas=boltz_apo)
    boltz_amp = Boltz2(proteins=sequences, ligands=amp, msas=boltz_apo)
    boltz_adp = Boltz2(proteins=sequences, ligands=adp, msas=boltz_apo)

    # Compare binding affinities
    merged = Panda(
        tables=[boltz_atp.tables.affinity, boltz_amp.tables.affinity,
                boltz_adp.tables.affinity],
        operations=[
            Panda.merge(on="id", prefixes=["atp_", "amp_", "adp_"]),
            Panda.sort("atp_affinity_pred_value", ascending=True)
        ]
    )
```

---

## Compound Library Screening

**File:** `compound_screening.py`

Protein–DNA–ligand co-folding with combinatorial compound library screening. Uses `CompoundLibrary` with a ChemDraw R-group enumeration file, `Bundle` for the homodimer, and `Each` to screen each compound individually.

```python
from biopipelines.pipeline import *
from biopipelines.compound_library import CompoundLibrary
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="TrpRepressor", job="CompoundLibraryScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    TrpR = Sequence("MAQQSPYSAAMA...EVLLKSD", ids="TrpR")
    DNA = Sequence("TGTACTAGTTAACTAGTAC", ids="dna")
    library = CompoundLibrary("./ExamplePipelines/compound_library_TrpR.cdxml")
    cofolded = Boltz2(proteins=Bundle(TrpR, TrpR),
                      dsDNA=DNA,
                      ligands=Each(library))
    merged = Panda(tables=[library.tables.compounds, cofolded.tables.affinity],
                   operations=[Panda.merge(),
                               Panda.calculate({"aff_uM": "10**affinity_pred_value"})])
    Plot(Plot.Scatter(data=merged.tables.result,
                      x="R1", y="aff_uM",
                      title="Predicted Affinity by R1 group",
                      xlabel="R1",
                      ylabel="Predicted Affinity [uM]"))
```

---

## Boltz2 with Combinatorics

**File:** `boltz2.py`

Comprehensive examples of Boltz2 with various input types and combinatorics (Bundle/Each). Demonstrates sequences, PDB structures, ligands, MSAs, and nested combinatorics.

```python
from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.combinatorics import Bundle, Each

with Pipeline(project="Examples", job="Boltz2",
              description="Boltz2 with various inputs"):

    Resources(gpu="A100", time="24:00:00", memory="16GB")

    protein_a = Sequence("MVLSPADKT...", ids="ProteinA")
    protein_b = Sequence("MNIFEMLRI...", ids="ProteinB")

    ligand_library = Ligand(['aspirin', 'caffeine', 'ibuprofen'])

    # Default: Each protein x Each ligand = 6 predictions
    boltz_each = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=ligand_library
    )

    # Bundle ligands: Each protein with all ligands = 2 predictions
    boltz_bundle = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # Nested: Each ligand bundled with a cofactor
    # Affinity calculated for the library ligand (first in bundle)
    atp = Ligand("ATP")
    boltz_nested = Boltz2(
        proteins=protein_a,
        ligands=Bundle(Each(ligand_library), atp)
    )
```

---

## Iterative Design with LigandMPNN + MutationComposer

**File:** `iterative_binding_optimization.py`

Multi-cycle iterative optimization: generate sequences with LigandMPNN, profile mutations, compose new sequences, fold with Boltz2, and select the best across cycles. Demonstrates the Panda tool for merging, filtering, and tracking results.

```python
from biopipelines.pipeline import *
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.mutation_composer import MutationComposer
from biopipelines.boltz2 import Boltz2
from biopipelines.distance import Distance
from biopipelines.panda import Panda

with Pipeline(project="Examples",
              job="LigandMPNN-MutationComposer-Cycle",
              description="Iterative optimization of binding affinity"):

    Resources(gpu="any", time="4:00:00", memory="16GB")

    HaloTag = Sequence("MAEIGTGFPF...", ids="HT")
    benzamide = Ligand(smiles=r"O=C(C1=CC=CC=C1)NCCOCCOCCCCCCCl", ids="BENZAMIDE")
    triazole = Ligand(smiles=r"ClCCCCCCOCCOCCN1C=C(C)N=N1", ids="TRIAZOLE")

    # Initial predictions
    original_benzamide = Boltz2(proteins=HaloTag, ligands=benzamide)
    original_triazole = Boltz2(proteins=HaloTag, ligands=triazole,
                               msas=original_benzamide)

    best_benzamide = original_benzamide

    for CYCLE in range(3):
        Suffix(f"Cycle{CYCLE+1}")

        lmpnn = LigandMPNN(structures=best_benzamide, ligand="LIG",
                           num_sequences=1000, batch_size=25,
                           redesigned="141+143+145+...")

        profiler = MutationProfiler(original=best_benzamide, mutants=lmpnn)
        composer = MutationComposer(frequencies=profiler.tables.absolute_frequencies,
                                    num_sequences=3, mode="weighted_random",
                                    max_mutations=3)

        boltz_benzamide = Boltz2(proteins=composer, ligands=benzamide)
        boltz_triazole = Boltz2(proteins=composer, ligands=triazole,
                                msas=boltz_benzamide)

        # Merge and filter results, select best across all cycles
        analysis = Panda(
            tables=[boltz_benzamide.tables.affinity, boltz_triazole.tables.affinity],
            operations=[
                Panda.merge(on="id", prefixes=["benzamide_", "triazole_"]),
                Panda.calculate({"affinity_delta":
                    "benzamide_affinity_pred_value - triazole_affinity_pred_value"})
            ]
        )
        # ... select best for next cycle
```

---

## De Novo Binder Design with BoltzGen

**File:** `denovo_dopamine_boltzgen_part1.py`

De novo protein binder design against a small molecule using BoltzGen. Runs parallel batches for high-throughput design.

```python
from biopipelines.pipeline import *
from biopipelines.ligand import Ligand
from biopipelines.boltzgen import BoltzGen

for batch in range(10):
    with Pipeline(project="Examples",
                  job="Dopamine_BoltzGen_1000designs",
                  description="De novo protein binder design against dopamine"):
        Resources(gpu="80GB|96GB", time="24:00:00", memory="16GB")
        Suffix(f"batch{batch}")
        designs = BoltzGen(
            ligand=Ligand(lookup="dopamine", ids="dopamine", codes="LIG"),
            binder_spec="140-180",
            protocol="protein-small_molecule",
            num_designs=1000,
            steps=["design", "inverse_folding", "folding",
                   "design_folding", "affinity"]
        )
```

---

## Tool Installation

**File:** `install_tools.py`

Install external tools and their environments via a single pipeline.

```python
from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.boltz2 import Boltz2
from biopipelines.pymol import PyMOL

with Pipeline(project="Setup", job="InstallTools",
              description="Install external tools and environments"):

    Resources(time="8:00:00", memory="32GB")

    RFdiffusion.install()      # Creates SE3nv, clones repo, downloads weights
    ProteinMPNN.install()      # Clones repo (uses SE3nv from RFdiffusion)
    AlphaFold.install()        # Downloads LocalColabFold installer
    Boltz2.install()           # Creates Boltz2Env
    PyMOL.install()            # Creates ProteinEnv
```
