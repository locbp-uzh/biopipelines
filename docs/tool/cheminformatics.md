# Cheminformatics

[← Back to Tool Reference](../tool_reference.md)

These tools work on the *chemistry* of small molecules: converting molecule files between formats, building 3-D coordinates, and computing descriptors. They read a **compounds** stream (the SMILES/codes a `Ligand` or `CompoundLibrary` carries) and may also read or write a **structures** stream (the actual coordinate files). See [The Ligand Contract](../developer_manual.md#the-ligand-contract-compounds--chemistry-structures--coordinates) for why chemistry and coordinates live in two separate streams.

---

### OpenBabel

Converts a molecule stream between chemical file formats, optionally adding hydrogens (pH-aware) and generating 3-D coordinates. This is the standard way to turn a `Ligand` (which carries only SMILES) into a docking-ready 3-D file: `OpenBabel(compounds=lig, convert_3d="sdf")`.

**References**: https://github.com/openbabel/openbabel

**Environment**: `biopipelines` (on Colab, install the `openbabel` pip extra).

**Parameters**:
- `compounds`: DataStream | StandardizedOutput = None — Source compounds (SMILES). Mutually exclusive with `structures`.
- `structures`: DataStream | StandardizedOutput = None — Source coordinate files. Mutually exclusive with `compounds`.
- `convert_3d`: str = None — Target 3-D coordinate format (`"sdf"`, `"mol2"`, `"mol"`, `"pdb"`, `"pdbqt"`, `"xyz"`). Produces a `structures` stream.
- `convert_1d`: str | List[str] = None — Line/notation format(s) (`"smi"`, `"inchi"`, `"cml"`). Added as columns on the `compounds` stream.
- `add_hydrogens`: bool = False — Add explicit hydrogens.
- `pH`: float = None — Add hydrogens at the given pH (implies `add_hydrogens=True`).
- `gen3d`: bool = False — Generate 3-D coordinates from a SMILES (compounds input only).
- `gen3d_quality`: str = "medium" — Embedding effort: `"fastest"`, `"fast"`, `"medium"`, `"better"`, `"best"`.
- `minimize`: bool = False — Force-field geometry minimization after embedding.
- `ff`: str = "MMFF94" — Force field for minimization (`"MMFF94"`, `"MMFF94s"`, `"UFF"`, `"GAFF"`, `"Ghemical"`).
- `minimize_steps`: int = 500 — Maximum minimization iterations.
- `use_structure_template`: bool = True — When converting from a `structures` input, use the input coordinates/connectivity as a bond-order template rather than re-perceiving from scratch.

**Streams**:
- `structures` — coordinate files (when `convert_3d` or `structures` input is used).
- `compounds` — chemistry passthrough (with extra columns when `convert_1d` is used).

**Example**:
```python
from biopipelines.openbabel import OpenBabel
from biopipelines.entities import Ligand

# Turn a SMILES-only Ligand into a 3-D SDF for docking-adjacent tools
aspirin = Ligand("aspirin")
sdf = OpenBabel(compounds=aspirin, convert_3d="sdf")
# sdf.streams.structures -> the SDF; sdf.streams.compounds -> chemistry passthrough

# Protonate a ligand at physiological pH
protonated = OpenBabel(compounds=aspirin, convert_3d="sdf", pH=7.4)
```

---

### RDKit

Computes per-compound cheminformatics descriptors from SMILES — molecular weight, logP, TPSA, hydrogen-bond donors/acceptors, rotatable bonds, QED, fraction sp³, and more. Useful for filtering or annotating a compound library before screening.

**References**: https://github.com/rdkit/rdkit

**Environment**: `biopipelines` (RDKit is pinned there — no extra installation).

**Parameters**:
- `compounds`: DataStream | StandardizedOutput (required) — Input compounds; SMILES are read from the compounds map_table.
- `descriptors`: List[str] = None — Subset of descriptor names to compute. `None` computes the default wide set.
- `morgan_fp`: bool = False — Also emit a Morgan fingerprint per compound.

**Tables**:
- `descriptors`:

  | id | smiles | MW | logP | TPSA | HBA | HBD | rotatable_bonds | QED | ... |
  |----|--------|----|------|------|-----|-----|-----------------|-----|-----|

**Example**:
```python
from biopipelines.rdkit_descriptors import RDKit
from biopipelines.compound_library import CompoundLibrary
from biopipelines.panda import Panda

library = CompoundLibrary("my_library.csv")
desc = RDKit(compounds=library)

# Lipinski-style filter downstream
druglike = Panda(
    tables=desc.tables.descriptors,
    operations=[Panda.filter("MW < 500 and logP < 5 and HBD <= 5 and HBA <= 10")],
)
```
