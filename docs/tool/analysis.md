# Analysis

[← Back to Tool Reference](../ToolReference.md)



### Distance

Measures distances between specific atoms and residues in structures. Useful for tracking ligand-protein interactions or structural features.

**Environment**: `biopipelines`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atom`: str (required) - Atom selection (e.g., 'LIG.Cl', 'name CA', 'A10.CA')
- `residue`: str (required) - Residue selection (e.g., 'D in IGDWG', '145', 'resn ALA')
- `method`: str = "min" - Distance calculation method (min, max, mean, closest)
- `metric_name`: str = None - Custom name for distance column in output (default: "distance")
- `unit`: str = "angstrom" - Output unit: "angstrom" (default) or "nm"

**Tables**:
- `distances`:

  | id | source_structure | {metric_name} | unit |
  |----|------------------|---------------|------|

**Example**:
```python
from biopipelines.distance import Distance

distances = Distance(
    structures=boltz,
    atom="LIG.Cl",
    residue="D in IGDWG",
    method="min",
    metric_name="chlorine_distance"
)
```


---

### Angle

Calculates angles between atoms in structures. Three modes are supported, selected by the shape of the `atoms` tuple. Useful for backbone phi/psi analysis, side chain rotamers, ligand geometry verification, and measuring relative orientations between structural elements.

**Environment**: `biopipelines`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atoms`: tuple (required) - Atom selections; form determines the angle mode (see below)
- `metric_name`: str = None - Custom name for angle column (default: `"angle"`, `"torsion"`, or `"vector_angle"` depending on mode)
- `unit`: str = "degrees" - Output unit: `"degrees"` or `"radians"`

**Selection Syntax**:
Same as Distance, with additional support for `residue.atom` format:
- `'10.CA'` - Alpha carbon of residue 10
- `'-1.C'` - Carbonyl carbon of last residue (C-terminus)
- `'LIG.C1'` - Atom C1 of ligand residue LIG
- `'D in IGDWG'` - Aspartic acid in sequence context (uses centroid if multiple atoms)

**Tables**:
- `angles`:

  | id | source_structure | {metric_name} | unit |
  |----|------------------|---------------|------|

**Angle Modes**:

| Form | Mode | Description | Range |
|------|------|-------------|-------|
| `(a, o, b)` | Bond | Angle at vertex o (a–o–b) | 0–180° |
| `(a, x1, x2, b)` | Torsional | Dihedral along axis x1–x2 | −180° to 180° |
| `((a1, a2), (b1, b2))` | Vector | Angle between vectors a1→a2 and b1→b2 | 0–180° |

Each selection may resolve to multiple atoms; centroids are used in that case.

**Example**:
```python
from biopipelines.angle import Angle

# Bond angle at CA (N-CA-C)
bond_angle = Angle(
    structures=boltz,
    atoms=('10.N', '10.CA', '10.C'),
    metric_name="nca_angle"
)

# Phi angle (C-N-CA-C)
phi = Angle(
    structures=boltz,
    atoms=('9.C', '10.N', '10.CA', '10.C'),
    metric_name="phi"
)

# Psi angle (N-CA-C-N)
psi = Angle(
    structures=boltz,
    atoms=('10.N', '10.CA', '10.C', '11.N'),
    metric_name="psi"
)

# Chi1 angle
chi1 = Angle(
    structures=boltz,
    atoms=('50.N', '50.CA', '50.CB', '50.CG'),
    metric_name="chi1"
)

# Ligand geometry
ligand_angle = Angle(
    structures=boltz,
    atoms=('LIG.C1', 'LIG.C2', 'LIG.C3')
)

# Angle between two structural vectors (e.g. helix axis proxies)
orientation = Angle(
    structures=boltz,
    atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
    metric_name="orientation"
)

# Output in radians (for use with cos/sin in Panda.calculate)
orientation_rad = Angle(
    structures=boltz,
    atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
    metric_name="orientation",
    unit="radians"
)
```

---

### DistanceSelector

Selects protein residues based on proximity to ligands or other reference points. Generates position specifications for downstream design tools.

**Installation**: It requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for distance reference
- `distance`: float = 5.0 - Distance cutoff in Angstroms
- `reference_type`: str = "ligand" - Type of reference (ligand, atoms, residues)
- `reference_selection`: str = "" - Specific PyMOL selection if not using ligand

**Tables**:
- `selections`:

  | id | pdb | within | beyond | distance_cutoff | reference_ligand |
  |----|-----|--------|--------|-----------------|------------------|

**Example**:
```python
from biopipelines.distance_selector import DistanceSelector

selector = DistanceSelector(
    structures=boltz,
    ligand="ATP",
    distance=8.0
)
```

---

### ConformationalChange

Quantifies structural changes between reference and target structures using PyMOL's alignment RMSD.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structures`: Union[DataStream, StandardizedOutput] (required) - Reference structures. Can be one or the same number as targets.
- `target_structures`: Union[DataStream, StandardizedOutput] (required) - Target structures to compare
- `selection`: Optional[str] = None - Residue range (e.g., '10-20+30-40'). None = whole structure.
- `alignment`: str = "align" - Alignment method (align, super, cealign). Rule of thumb: sequence similarity > 50% -> align; otherwise cealign.
- `atoms`: str = "all" - Which atoms to use for alignment:
  - `"all"` (default): all atoms
  - `"CA"`: alpha-carbon only
  - `"backbone"`: backbone atoms (CA+C+N+O)
  - Any `+`-separated atom names, e.g. `"CA+CB"`

**Tables**:
- `changes`:

  | id | reference_structure | target_structure | selection | num_aligned_atoms | RMSD |
  |----|---------------------|------------------|-----------|-------------------|------|

**Example**:
```python
from biopipelines.conformational_change import ConformationalChange

# All-atom RMSD on a specific region
conf_change = ConformationalChange(
    reference_structures=apo_structures,
    target_structures=holo_structures,
    selection="10-50",
    alignment="super"
)

# CA-only RMSD
conf_change = ConformationalChange(
    reference_structures=design,
    target_structures=refolded,
    atoms="CA"
)

# Backbone RMSD on redesigned region
conf_change = ConformationalChange(
    reference_structures=kinase,
    target_structures=refolded,
    selection=backbones.tables.structures.fixed,
    atoms="backbone"
)
```

---

### Contacts
**Environment**: `ProteinEnv`

Analyzes contacts between selected protein regions and ligands. For each selected residue, calculates the minimum distance to any ligand atom. Returns contact counts and distance statistics.

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `selections`: Union[str, ToolOutput] = None - Protein region selections (string format: '10-20+30-40', table reference, or None for all protein)
- `ligand`: str (required) - Ligand residue name (3-letter code, e.g., 'LIG', 'ATP', 'GDP')
- `contact_threshold`: float = 5.0 - Distance threshold for counting contacts (Angstroms)
- `contact_metric_name`: str = None - Custom name for contact count column (default: "contacts")

**Tables**:
- `contacts`:

  | id | source_structure | selections | ligand | contacts | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized |
  |----|------------------|------------|--------|----------|--------------|--------------|---------------|-------------------------------|

**Output Columns**:
- `id`: Structure identifier
- `source_structure`: Path to input structure file
- `selections`: Protein residues analyzed (e.g., '10-20+30-40' or 'all_protein')
- `ligand`: Ligand residue name
- `contacts`: Number of residues within contact_threshold distance
- `min_distance`: Minimum distance from any selected residue to ligand (Å)
- `max_distance`: Maximum distance from any selected residue to ligand (Å)
- `mean_distance`: Mean distance from selected residues to ligand (Å)
- `sum_distances_sqrt_normalized`: Sum of distances divided by √(number of residues)

**Example**:
```python
from biopipelines.contacts import Contacts

# Analyze contacts with specific protein regions
contacts = Contacts(
    structures=rfdaa,
    selections=rfdaa.tables.structures.designed,
    ligand="LIG",
    contact_threshold=5.0
)

# Use fixed selection for all structures
contacts = Contacts(
    structures=boltz,
    selections='10-20+30-40',
    ligand="ATP",
    contact_threshold=4.0
)

# Analyze all protein residues
contacts = Contacts(
    structures=boltz,
    ligand="GDP"
)
```

---

### PoseBusters

Validates computationally generated molecule poses by checking bond lengths, bond angles, internal steric clashes, volume overlap with protein, and more. Supports `dock` mode (ligand + protein) and `redock` mode (+ reference ligand for RMSD comparison).

**Environment**: `posebusters`

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Protein-ligand complexes (PDB/CIF)
- `ligand`: str (required) - 3-letter residue code for the ligand (e.g., 'LIG', 'ATP')
- `reference_ligand`: Union[DataStream, StandardizedOutput, None] = None - Reference ligand structure for redock mode
- `reference_ligand_code`: Optional[str] = None - Residue code in reference structure (defaults to `ligand`)
- `mode`: str = "dock" - Validation mode: 'dock' or 'redock' (auto-set to 'redock' if reference_ligand provided)
- `check`: Union[str, List[str]] = "pose" - Which checks contribute to `all_pass`. All check columns are still written to the CSV — excluded ones are placed to the right of `all_pass` so the audit trail is preserved.

  Presets:
  - `"pose"` (default): every check except `no_radicals`. That one check fails spuriously on predicted poses because RDKit's bond perception from coordinates flags unsaturated atoms as radicals.
  - `"all"`: every available check (PoseBusters default behaviour).
  - `"loading"`: only that the molecule loaded.
  - `"geometry"`: loading + bond lengths/angles/clashes/flatness/energy.
  - `"chemistry"`: loading + sanity + radicals + connectivity.

  List form: `check=["all_atoms_connected", "bond_lengths"]` to compute `all_pass` over a custom subset.

**Tables**:
- `posebusters` (dock mode):

  | id | source_structure | mol_pred_loaded | mol_cond_loaded | sanitization | inchi_convertible | all_atoms_connected | bond_lengths | bond_angles | internal_steric_clash | aromatic_ring_flatness | non-aromatic_ring_non-flatness | double_bond_flatness | internal_energy | protein-ligand_maximum_distance | minimum_distance_to_protein | minimum_distance_to_organic_cofactors | minimum_distance_to_inorganic_cofactors | minimum_distance_to_waters | volume_overlap_with_protein | volume_overlap_with_organic_cofactors | volume_overlap_with_inorganic_cofactors | volume_overlap_with_waters | all_pass | no_radicals |
  |----|------------------|-----------------|-----------------|--------------|-------------------|---------------------|--------------|-------------|-----------------------|------------------------|--------------------------------|----------------------|-----------------|--------------------------------|-----------------------------|---------------------------------------|-----------------------------------------|----------------------------|-----------------------------|---------------------------------------|-----------------------------------------|----------------------------|----------|-------------|

  Columns to the left of `all_pass` are the checks selected by `check=` (default preset `"pose"` shown). Columns to the right are still computed but excluded from `all_pass`. In `redock` mode, additional columns `mol_true_loaded`, `molecular_formula`, `molecular_bonds`, `double_bond_stereochemistry`, `tetrahedral_chirality`, `stereochemistry_preserved`, and `rmsd_≤_2å` are also included.

**Output Columns**:
- `id`: Structure identifier (suffixed with `_lig1`, `_lig2` etc. when multiple ligand copies exist)
- `source_structure`: Path to input structure file
- Boolean check columns: Each PoseBusters test (True = pass, False = fail)
- `all_pass`: True only if all checks selected by `check=` pass

**Example**:
```python
from biopipelines.posebusters import PoseBusters

# Default — "pose" preset, drops electronics so all_pass reflects geometry/placement only
validation = PoseBusters(
    structures=boltz_holo,
    ligand="LIG"
)

# Custom subset — all_pass over only these two columns
validation = PoseBusters(
    structures=boltz_holo,
    ligand="LIG",
    check=["all_atoms_connected", "bond_lengths"]
)

# Reproduce the original PoseBusters strict behaviour
validation = PoseBusters(
    structures=boltz_holo,
    ligand="LIG",
    check="all"
)

# Redock mode — validate against crystal reference
xrc = PDB("1ABC")
validation = PoseBusters(
    structures=boltz_holo,
    ligand="ATP",
    reference_ligand=xrc,
    reference_ligand_code="ATP",
    mode="redock"
)
```

---

### PoseChange

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: str (required) - Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
- `sample_ligand`: Optional[str] = None - Ligand residue name in sample structures (default: same as reference_ligand)
- `reference_alignment`: Optional[str] = None - PyMOL selection for reference structure alignment (default: "not resn {reference_ligand}")
- `target_alignment`: Optional[str] = None - PyMOL selection for target structure alignment (default: "not resn {sample_ligand}")

**Tables**:
- `changes`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | alignment_rmsd | num_ligand_atoms | reference_alignment | target_alignment |
  |----|------------------|---------------------|-------------|-------------------|----------------|------------------|---------------------|------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand
- `reference_alignment`: PyMOL selection used for reference structure alignment
- `target_alignment`: PyMOL selection used for target structure alignment

**Example**:
```python
from biopipelines.pose_change import PoseChange
from biopipelines.pdb import PDB
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda

# Compare designed structures to XRC reference
xrc = PDB(pdbs="4ufc", ids="reference")
designed = Boltz2(proteins=sequences, ligands="CCO")

pose_analysis = PoseChange(
    reference_structure=xrc,
    sample_structures=designed,
    reference_ligand="ATP",
    sample_ligand="LIG",
    # By default aligns on "not resn ATP" (reference) and "not resn LIG" (target)
    # Override individually if needed:
    # reference_alignment="chain A and not resn ATP",
    # target_alignment="chain A and not resn LIG"
)

# Filter structures with RMSD < 2.0 Å
good_poses = Panda(
    tables=pose_analysis.tables.changes,
    pool=designed,
    operations=[Panda.filter("ligand_rmsd < 2.0")]
)
```

---

### CABSflex

Fast coarse-grained Monte Carlo simulation of protein structure flexibility. Produces conformational ensembles and per-residue RMSF profiles. Optionally rebuilds models to all-atom representation using MODELLER.

**WARNING**: MODELLER requires a license key (free for academics). Get one at https://salilab.org/modeller/registration.html and set `KEY_MODELLER` before running.

**Environment**: `CABSflex` (Python 2.7 with `modeller` and `cabs`)

**Installation**: `conda create -n CABSflex python=2.7 && conda install -c salilab modeller -c lcbio cabs`. On Colab: `pip install modeller cabs`.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Input protein structures
- `num_models`: int = 10 - Number of cluster medoids / final models per input structure
- `mc_cycles`: int = 50 - Monte Carlo cycles between trajectory frames
- `mc_steps`: int = 50 - Monte Carlo steps
- `mc_annealing`: int = 20 - Temperature annealing cycles
- `temperature`: Optional[str] = None - Temperature range as "TINIT TFINAL" (default: "1.4 1.4")
- `flexibility`: Optional[str] = None - Residue flexibility: float (0=flexible, 1=stiff), 'bf', 'bfi', 'bfg', or filename
- `filtering_count`: int = 1000 - Number of low-energy models for clustering
- `aa_rebuild`: bool = True - Rebuild to all-atom with MODELLER
- `restraints`: str = "ss2" - Restraint mode: 'all', 'ss1', 'ss2'
- `restraints_gap`: int = 3 - Min gap along chain for restraints
- `restraints_min`: float = 3.8 - Min distance in Å for restraints
- `restraints_max`: float = 8.0 - Max distance in Å for restraints
- `weighted_fit`: str = "gauss" - Fit method: 'gauss', 'flex', 'ss', 'off', or filename

**Streams**:
- `structures`: PDB ensemble models (`num_models` per input structure)
- `images`: SVG plots (RMSF, RMSD, energy) per input structure
- `rmsf`: Per-residue RMSF CSVs (`per-residue-values-csv` format, columns: id, chain, resi, rmsf), one file per input structure

**Tables**:
- `structures`:

  | id | file | structures.id |
  |----|------|---------------|

- `rmsf_all` (merged across all input structures):

  | id | chain | resi | rmsf |
  |----|-------|------|------|

**Example**:
```python
from biopipelines.cabsflex import CABSflex

# Basic flexibility simulation
protein = PDB("1AHN")
flex = CABSflex(structures=protein, num_models=10)

# Access merged RMSF table
flex.tables.rmsf_all  # all structures combined

# Access per-structure RMSF stream
for rmsf in flex.streams.rmsf:
    print(f"{rmsf.ids[0]}: {rmsf.files[0]}")

# Access ensemble models
flex.streams.structures  # 10 PDB models per input

# Custom simulation parameters
flex = CABSflex(
    structures=protein,
    num_models=5,
    mc_cycles=100,
    mc_annealing=30,
    temperature="1.4 2.0",
    flexibility=0.5,
    aa_rebuild=True
)
```

### SASA

Solvent-accessible surface area analysis. Computes the change in SASA of a ligand when bound to vs separated from its protein partner — `delta_SASA = SASA(ligand alone) − SASA(ligand in complex)`. A larger delta indicates more ligand surface buried by the binder, typically interpreted as tighter packing.

**Environment**: `ProteinEnv` (installed alongside `PyMOL`).

**Installation**: no extra step beyond `PyMOL.install()`; SASA reuses the same env.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Protein-ligand complexes.
- `ligand`: str (required) — Ligand residue identifier in the structure (e.g., `"LIG"`, `"AMX"`, `":X:"`).
- `dot_density`: int = 4 — PyMOL `get_area` dot density (1–4, higher = more accurate, slower).

**Tables**:
- `sasa`:

  | id | structure | sasa_ligand_alone | sasa_ligand_complex | delta_sasa |
  |----|-----------|-------------------|---------------------|------------|

**Example**:
```python
from biopipelines.sasa import SASA

complex = PDB("9RTM", convert="pdb")
sasa = SASA(structures=complex, ligand="LIG")
sasa.tables.sasa  # delta-SASA per structure
```

---

### ADMETAI

ADMET endpoint predictions for compound libraries via the ADMET-AI Chemprop-RDKit model. Reads SMILES from a compounds stream and writes one row per input compound with all upstream-reported ADMET properties (~40+ regression and classification scores covering absorption, distribution, metabolism, excretion, and toxicity).

**Environment**: `admet_ai`

**Installation**: `ADMETAI.install()` creates a dedicated conda env (Python 3.10) and installs `admet-ai` via pip. Verification instantiates `ADMETModel`, which downloads the bundled weights and warms the cache.

**Parameters**:
- `compounds`: Union[DataStream, StandardizedOutput] (required) — Input compounds. SMILES are read from the `smiles` column of the compounds map_table.

**Tables**:
- `admet`:

  | id | smiles | <ADMET endpoint columns from upstream model> |
  |----|--------|----------------------------------------------|

  The endpoint column set is determined by the installed `admet-ai` version's `ADMETModel.predict()` output and is preserved verbatim in the CSV.

**Example**:
```python
from biopipelines.admet_ai import ADMETAI
from biopipelines.compound_library import CompoundLibrary

with Pipeline("Project", "ADMET", description="Library ADMET screen"):
    Resources(memory="8GB", time="1:00:00")
    library = CompoundLibrary("my_library.csv")
    admet = ADMETAI(compounds=library)

# Filter for high oral bioavailability candidates downstream
from biopipelines.panda import Panda
filtered = Panda(
    tables=admet.tables.admet,
    operations=[Panda.filter("HIA_Hou > 0.8")],
)
```

**Reference**: Swanson et al. (2024) ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. *Bioinformatics* 40, btae416. https://github.com/swansonk14/admet_ai
