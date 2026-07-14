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

Two independent jobs, selected by which inputs are given. With `compounds=` it computes per-compound cheminformatics descriptors from SMILES — molecular weight, logP, TPSA, hydrogen-bond donors/acceptors, rotatable bonds, QED, fraction sp³, and more — for filtering or annotating a compound library before screening. With `structures=` it computes the **conformer strain** of each posed ligand: how much internal (torsional) energy the bound conformation carries relative to a relaxed one. Both may be given together.

Strain answers a question PoseBusters cannot: a pose can be clash-free and geometrically valid yet still be twisted into an energetically implausible conformation. It is *internal* strain, not interaction energy — for the latter see XTB.

**References**: https://github.com/rdkit/rdkit

**Environment**: `biopipelines` (RDKit is pinned there — no extra installation).

**Parameters**:
- `compounds`: DataStream | StandardizedOutput = None — Input compounds; SMILES are read from the compounds map_table. Produces the `descriptors` table.
- `structures`: DataStream | StandardizedOutput = None — Posed ligand coordinate files, one per id (e.g. the `structures` stream of `Ligand(code=..., structures=poses)`, which extracts a bound ligand *with* its coordinates). Produces the `strain` table.
- `descriptors`: List[str] = None — Subset of descriptor names to compute. `None` computes the default wide set.
- `morgan_fp`: bool = False — Also emit a Morgan fingerprint per compound.
- `smiles`: str | (TableInfo, column) = None — Bond-order template for the posed coordinates: a literal SMILES, or a table-column reference for a per-id template. Defaults to the `smiles` column of the `compounds` stream when one is supplied. **Mandatory for the strain path** — coordinate-only perception mis-assigns conjugated and charged systems, so there is no silent fallback; supplying neither raises.
- `restrain_bonds`: List[Tuple[str, str]] = None — Explicit torsions to restrain, as `(atom_name, atom_name)` pairs naming the two central atoms of each bond, matching the PDB atom names (e.g. `[("C64", "C72"), ("N67", "C42")]`). Default `None` restrains every rotatable bond.
- `ff`: str = "auto" — Force field: `"auto"` (MMFF94, falling back to UFF where MMFF has no parameters), `"mmff"`, or `"uff"`.

**Tables**:
- `descriptors` (when `compounds=` is given):

  | id | smiles | MW | logP | TPSA | HBA | HBD | rotatable_bonds | QED | ... |
  |----|--------|----|------|------|-----|-----|-----------------|-----|-----|

- `strain` (when `structures=` is given):

  | id | smiles | e_pose | e_relaxed | strain | ff_engine |
  |----|--------|--------|-----------|--------|-----------|

  `e_pose` — force-field energy at the pose, after a torsion-restrained minimisation (kcal/mol). `e_relaxed` — energy after a free local minimisation from the same starting coordinates. `strain` — `e_pose - e_relaxed`. `ff_engine` — `MMFF` or `UFF`, whichever actually typed the molecule.

#### Why the reference state is a torsion-restrained minimum

**Do not "simplify" this to `E(pose) - E(free minimum)`.** That naive form is wrong and produces meaningless numbers. Ligand coordinates from a neural predictor (Boltz) or a docking program carry bond lengths and angles that differ slightly from force-field ideals. Minimisation releases that mismatch as a large constant "strain" that has nothing to do with conformation — measured over 257 real poses it puts a **57–95 kcal/mol floor** (median 69) under every pose, burying the real signal.

The procedure that works:

1. Build the molecule from the pose coordinates and assign bond orders from a SMILES template (`AssignBondOrdersFromTemplate`); add hydrogens with `addCoords=True`.
2. `e_pose` — minimise a copy with the rotatable torsions **restrained at their current values**, then score the result on a force field rebuilt *without* the restraint term. Bond lengths and angles relax off the geometry mismatch while the conformation is held.
3. `e_relaxed` — free local minimisation from the same starting coordinates. Deliberately the **nearest** local minimum, not a global conformer search: a global reference would add a constant bulk offset to every pose.
4. `strain = e_pose - e_relaxed`.

Both minimisations must converge. One that hits its iteration cap has not reached a minimum, so the difference of the two energies is not a strain — that pose raises and is reported as a failure rather than emitting a plausible-looking number.

With this reference the same 257 poses give a sane distribution — 45% at ≈0–1 kcal/mol, bulk under 5, a thin tail to ~13, a clean gap, then 13 clearly-junk poses at 30–89. If an implementation reproduces a ~60 kcal/mol floor, the reference state is wrong.

**Thresholds.** Real bound ligands typically sit at 0–3 kcal/mol. The literature red-flag for an implausible pose is >10 kcal/mol, but that is a permissive outlier detector; **~5 kcal/mol is a sensible acceptance gate** (10 kcal/mol ≈ 17 RT — roughly a whole binding-energy budget spent on internal strain). Small negative values are expected and normal.

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

Strain-gating a pose set — extract the ligand *with* its coordinates, score it, then keep only poses that are both geometrically valid and unstrained:

```python
from biopipelines.entities import Ligand

posed = Ligand(code="LIG", structures=valid_poses)     # extract ligand WITH its coordinates
strain = RDKit(structures=posed, smiles=DYE_SMILES)

kept = Panda(
    tables=[pb.tables.posebusters, strain.tables.strain],
    operations=[Panda.merge(on="id"),
                Panda.filter("all_pass == True and strain < 5")],
    pool=valid_poses,
)
```
