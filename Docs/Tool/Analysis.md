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

**Tables**:
- `distances`:

  | id | source_structure | {metric_name} |
  |----|------------------|---------------|

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

Calculates bond angles (3 atoms) or torsional/dihedral angles (4 atoms) between specified atoms in structures. Useful for backbone phi/psi analysis, side chain rotamers, and ligand geometry verification.

**Environment**: `biopipelines`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atoms`: List[str] (required) - List of 3 or 4 atom selections
- `metric_name`: str = None - Custom name for angle column (default: "angle" or "torsion")

**Selection Syntax**:
Same as Distance, with additional support for `residue.atom` format:
- `'10.CA'` - Alpha carbon of residue 10
- `'-1.C'` - Carbonyl carbon of last residue (C-terminus)
- `'LIG.C1'` - Atom C1 of ligand residue LIG
- `'D in IGDWG'` - Aspartic acid in sequence context (uses centroid if multiple atoms)

**Tables**:
- `angles`:

  | id | source_structure | {metric_name} |
  |----|------------------|---------------|

**Angle Types**:
- **3 atoms (A-B-C)**: Bond angle at B in degrees (0-180°)
- **4 atoms (A-B-C-D)**: Torsional angle in degrees (-180° to 180°)

**Example**:
```python
from biopipelines.angle import Angle

# Bond angle at CA (N-CA-C angle)
bond_angle = Angle(
    structures=boltz,
    atoms=['10.N', '10.CA', '10.C'],
    metric_name="nca_angle"
)

# Phi angle (C-N-CA-C)
phi = Angle(
    structures=boltz,
    atoms=['9.C', '10.N', '10.CA', '10.C'],
    metric_name="phi"
)

# Psi angle (N-CA-C-N)
psi = Angle(
    structures=boltz,
    atoms=['10.N', '10.CA', '10.C', '11.N'],
    metric_name="psi"
)

# Chi1 angle for a residue
chi1 = Angle(
    structures=boltz,
    atoms=['50.N', '50.CA', '50.CB', '50.CG'],
    metric_name="chi1"
)

# Ligand geometry
ligand_angle = Angle(
    structures=boltz,
    atoms=['LIG.C1', 'LIG.C2', 'LIG.C3']
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
- `selection`: Optional[str] = None - Residue range (e.g., '10-20+30-40'). None = all CA atoms.
- `alignment`: str = "align" - Alignment method (align, super, cealign). Rule of thumb: sequence similarity > 50% -> align; otherwise cealign.

**Tables**:
- `changes`:

  | id | reference_structure | target_structure | selection | num_aligned_atoms | RMSD |
  |----|---------------------|------------------|-----------|-------------------|------|

**Example**:
```python
from biopipelines.conformational_change import ConformationalChange

conf_change = ConformationalChange(
    reference_structures=apo_structures,
    target_structures=holo_structures,
    selection="10-50",
    alignment="super"
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

### PoseChange

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: str (required) - Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
- `sample_ligand`: Optional[str] = None - Ligand residue name in sample structures (default: same as reference_ligand)
- `alignment_selection`: str = "protein" - PyMOL selection for protein alignment (e.g., "chain A", "backbone")
- `calculate_centroid`: bool = True - Calculate ligand centroid distance
- `calculate_orientation`: bool = False - Calculate orientation angle difference

**Tables**:
- `changes`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | alignment_rmsd | num_ligand_atoms | alignment_selection |
  |----|------------------|---------------------|-------------|-------------------|----------------|------------------|---------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand

**Example**:
```python
from biopipelines.pose_change import PoseChange
from biopipelines.pdb import PDB
from biopipelines.boltz2 import Boltz2
from biopipelines.filter import Filter

# Compare designed structures to XRC reference
xrc = PDB(pdbs="4ufc", ids="reference")
designed = Boltz2(proteins=sequences, ligands="CCO")

pose_analysis = PoseChange(
    reference_structure=xrc,
    sample_structures=designed,
    reference_ligand="ATP",
    sample_ligand="LIG",
    alignment_selection="chain A and backbone"
)

# Filter structures with RMSD < 2.0 Å
good_poses = Filter(
    data=pose_analysis.tables.changes,
    pool=designed,
    expression="ligand_rmsd < 2.0"
)
```
