# Cheminformatics

[← Back to Tool Reference](../tool_reference.md)

These tools work on the *chemistry* of small molecules: converting molecule files between formats, building 3-D coordinates, computing descriptors, planning syntheses, and generating new molecules. They read a **compounds** stream (the SMILES/codes a `Ligand` or `CompoundLibrary` carries) and may also read or write a **structures** stream (the actual coordinate files). See [The Ligand Contract](../developer_manual.md#the-ligand-contract-compounds--chemistry-structures--coordinates) for why chemistry and coordinates live in two separate streams.

---

### OpenBabel

Converts a molecule stream between chemical file formats, optionally adding hydrogens (pH-aware) and generating 3-D coordinates. This is the standard way to turn a `Ligand` (which carries only SMILES) into a docking-ready 3-D file: `OpenBabel(compounds=lig, convert_3d="sdf")`.

**Tags**: data, refine-structure, small-molecule

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

**Tags**: measure, small-molecule, energy

**References**: https://github.com/rdkit/rdkit

**Environment**: `biopipelines` (RDKit is pinned there — no extra installation).

**Parameters**:
- `compounds`: DataStream | StandardizedOutput = None — Input compounds; SMILES are read from the compounds map_table. Produces the `descriptors` table.
- `structures`: DataStream | StandardizedOutput = None — Posed **ligand-only** coordinate files, one per id (e.g. the `structures` stream of `Ligand(structures=poses, codes="UNL")`, which carves the bound ligand out of every input structure *with* its coordinates). Produces the `strain` table. A whole complex is refused rather than minimised: bond-order templating matches the ligand substructure inside a complex without erroring, so every protein atom would be silently dragged into the minimisation.
- `descriptors`: List[str] = None — Subset of descriptor names to compute. `None` computes the default wide set.
- `morgan_fp`: bool = False — Also emit a Morgan fingerprint per compound.
- `smiles`: str | (TableInfo, column) = None — Bond-order template for the posed coordinates: a literal SMILES, or a table-column reference for a per-id template. Defaults to the `smiles` column of the `compounds` stream when one is supplied. **Mandatory for the strain path** — coordinate-only perception mis-assigns conjugated and charged systems, so there is no silent fallback; supplying neither raises.
- `restrain_bonds`: List[Tuple[str, str]] = None — Explicit torsions to restrain, as `(atom_name, atom_name)` pairs naming the two central atoms of each bond, matching the PDB atom names (e.g. `[("C64", "C72"), ("N67", "C42")]`). Default `None` restrains every rotatable bond.
- `ff`: str = "auto" — Force field: `"auto"` (MMFF94, falling back to UFF where MMFF has no parameters), `"mmff"`, or `"uff"`.
- `max_iters`: int = 2000 — Minimizer iteration cap per pose, for both the restrained and free minimizations. A pose that hits the cap has not reached a minimum, so its strain would be unreliable — it raises rather than returning a misleading number (see "Both minimisations must converge" above). Larger or more strained ligands may need more than the default; raise this rather than trusting a capped result.

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

posed = Ligand(codes="LIG", structures=valid_poses)     # extract ligand WITH its coordinates
strain = RDKit(structures=posed, smiles=DYE_SMILES)

kept = Panda(
    tables=[pb.tables.posebusters, strain.tables.strain],
    operations=[Panda.merge(on="id"),
                Panda.filter("all_pass == True and strain < 5")],
    pool=valid_poses,
)
```

---

### AiZynthFinder

Retrosynthetic route planning. Runs a Monte-Carlo tree search that recursively applies reaction templates to break each target molecule apart, stopping when every branch ends in a purchasable building block. Use it to ask whether a designed compound is actually *makeable*, and out of what — a synthetic-accessibility filter with the routes and starting materials attached, rather than a single opaque score.

The output is nested, and the tool exposes all three levels: a **target** is searched, producing several ranked **routes**, each of which ends in several **precursors**. Routes are ranked best-first by the state score (a combination of stock availability and route length) unless `route_scorer` picks another; precursors carry only an `in_stock` boolean, so they are molecules rather than scored objects.

The `precursors` stream is a **compounds** stream with the canonical `id | format | smiles | ccd` schema, so the building blocks feed directly into any compounds consumer (`ADMETAI`, `RDKit`, `OpenBabel`, a docking run). Its map_table *is* the `precursors` table, so filtering it with `Panda` and passing `pool=` carries the stream along.

**Tags**: predict-property, small-molecule

**References**: https://github.com/MolecularAI/aizynthfinder — Genheden et al. (2020) *J. Cheminform.* 12, 70; Saigiridharan et al. (2024) *J. Cheminform.* 16, 57 (AiZynthFinder 4.0).

**Environment**: `aizynthfinder` (dedicated, Python 3.11; upstream 4.x supports 3.10 to 3.12). `AiZynthFinder.install()` creates the env from the unpinned `aizynthfinder` package and runs `download_public_data`, which fetches the expansion policies (`uspto`, `ringbreaker`), reaction templates and quick-filter model from Zenodo and the ZINC stock from figshare, then writes the `config.yml` the CLI reads. An existing env is kept, so a download that failed halfway is retried by calling `install()` again; `force_reinstall=True` rebuilds both. The install is verified by starting `aizynthcli`, which imports the whole package, and by the presence of `config.yml`.

**Parameters**:
- `compounds`: DataStream | StandardizedOutput — Target molecules. SMILES are read from the `smiles` column of the compounds map_table.
- `config`: str = None — Path to an aizynthfinder `config.yml`, overriding the downloaded one. This is how your own expansion model, filter model or stock enters: define it in the file's `expansion:`, `filter:` or `stock:` section (see the upstream configuration guide), then select it by name with the parameters below.
- `stocks`: List[str] = None — Stock names to search against (default: all in the config).
- `expansion_policy`: List[str] = None — Expansion policy names (default: the first in the config, `uspto` in the public data; `ringbreaker` is the other one it defines).
- `filter_policy`: List[str] = None — Filter policy names (default: all in the config).
- `max_routes`: int = 25 — Most routes kept per target.
- `min_routes`: int = 5 — Fewest routes kept per target.
- `time_limit`: int = None — Per-target search cutoff in seconds.
- `iteration_limit`: int = None — Per-target cap on tree-search iterations.
- `nproc`: int = None — Targets searched in parallel. Upstream then writes one `aizynthcli<N>.log` per process into the working directory, which is where a skipped target's reason is printed.
- `search_algorithm`: str = None — `"mcts"` (Monte Carlo tree search, the upstream default), `"retrostar"`, `"breadth_first"` or `"dfpn"`. Left unset, the config decides.
- `max_transforms`: int = None — Deepest route searched, in reaction steps (upstream default 6).
- `route_scorer`: str = None — Scorer that ranks the routes and defines `targets.top_score` (default `"state score"`). One of the names under *Scorers* below, or one your own `config` defines.
- `scorers`: List[str] = None — Further scorers computed for every returned route, each adding a column to `routes`. Any of `"fraction in stock"`, `"average template occurrence"`, `"sum of prices"`, `"route cost"`.
- `filter_cutoff`: float = None — Feasibility threshold of the quick-filter, in (0, 1]: a reaction the filter scores below it is rejected (upstream default 0.05). Applied to every quick-filter in the config.
- `frozen_substructures`: List[str] = None — SMARTS patterns that must survive every reaction; a reaction that splits one of them is rejected. A pattern that does not occur in the molecule being expanded is ignored for that molecule. Not supported by `"dfpn"`.
- `stock_stop_criteria`: dict = None — Extra conditions a stock molecule must meet to count as purchasable: `{"weight": <max exact mass, g/mol>, "price": <max price>, "counts": {"C": <max atoms of that element>}}`.

**Search algorithms**. `mcts` is the original search. `retrostar` is a best-first AND/OR search guided by the policy priors; with the default zero molecule cost it does not weigh the starting materials, and a learned Retro* cost model needs a config file. `breadth_first` expands every open molecule one level at a time and so is exhaustive but grows quickly with depth. `dfpn` is a depth-first proof-number search and does not apply a filter policy. Each reports the same `targets` columns, but `number_of_routes` counts what that algorithm found, which is not comparable across algorithms.

**Scorers**. Every route is always scored by `state score`, `number of reactions`, `number of pre-cursors` and `number of pre-cursors in stock`; these fill `score`, `number_of_steps`, `number_of_precursors` and `number_of_precursors_in_stock`. `route_scorer` may be any of them or of the four extra scorers, and the direction is upstream's: lowest first for `number of reactions`, `number of pre-cursors`, `sum of prices` and `route cost`, highest first for the rest. The extra scorers are `fraction in stock` (share of the precursors that are purchasable), `average template occurrence` (mean frequency of the templates used, higher means better-precedented chemistry), `sum of prices` and `route cost` (Badowski et al. 2019: reaction cost and yield on top of the precursor prices). Prices come from the stock, and the public ZINC stock carries none, so the two price scorers fall back to a default cost per precursor; give a `config` whose stock entry sets `price_col` for meaningful values. A scorer that does not look at whether the route is solved (`number of reactions`, `number of pre-cursors`) can rank an unsolved route first, so filter `routes` on `is_solved` when ranking by one. Each extra scorer adds a `routes` column named after it with spaces as underscores (`route cost` becomes `route_cost`); `score` is always the state score, so ranking by another scorer reorders the `route` ranks but leaves `score` as is. The scorers that need extra model or data files (route similarity, expert-augmented, reaction class rank) and the search rewards of `mcts` (single, weighted or multi-objective) are only reachable through your own `config`.

**Filters and stock**. The quick-filter is the one in the public data (`uspto`); `filter_policy` selects among the filters the config defines, and `frozen_substructures` adds one more, named `frozen_substructure`, that is always active even when `filter_policy` lists others. `stock_stop_criteria` is upstream's `stop_criteria` stock setting: a molecule found in a stock only counts as purchasable if it also passes every criterion (`weight` at most, `counts` at most per element, `price` at most; `price` needs a stock with prices). It changes which molecules end a branch, so it changes `in_stock` in `precursors` too.

**Expansion policies and your own model**. `expansion_policy` selects which policies drive the search, and several may be listed. A user-supplied template-based model is a `config` entry with its model and template files; any other upstream expansion strategy, including combining policies with weights, is the same route through the `expansion:` section. The Chemformer and ModelZoo strategies ship in the upstream repository's `plugins/` folder rather than in the pip package and need their own service or packages, so they are not set up by `install()`. Not exposed at all: atom-index based `break_bonds` / `freeze_bonds` (the indices are per molecule, which a batched run cannot give), route clustering (it needs upstream's `all` extra, which pulls in torch) and aizynthcli's checkpoint restart.

**Streams**:
- `routes` — one JSON per route holding its reaction tree. ids are `<compound>_<route>`.
- `precursors` — compounds stream of route leaf molecules. ids are `<compound>_<route>_<precursor>`.

**Tables**:
- `targets` — `id | smiles | is_solved | top_score | number_of_routes | number_of_returned_routes | number_of_solved_routes | number_of_nodes | search_time | first_solution_time | first_solution_iteration`. One row per input compound; a target with no route is kept with `is_solved=False`. `number_of_routes` is how many routes the search found, `number_of_returned_routes` how many `max_routes` let through — the latter is how many route ids exist.
  `top_score` is the `route_scorer` score of the best route. `first_solution_time` and `first_solution_iteration` are 0 when the search never flagged a solution, as upstream reports them; `breadth_first` flags one only when every branch it expanded ends in stock, so it can report 0 for a target it does solve. A SMILES RDKit cannot sanitize is never searched but keeps its row, with everything but `id` and `smiles` empty.
- `routes` — `id | compounds.id | route | score | is_solved | number_of_steps | number_of_precursors | number_of_precursors_in_stock`, then one column per extra scorer (`scorers`, plus `route_scorer` when it is not one of the four always-computed ones). `score` is the state score.
- `precursors` — `id | format | smiles | ccd | compounds.id | route | precursor | in_stock`.
- `missing` — `id | removed_by | kind | cause`. Targets with no route (and compounds with no SMILES) land here as `kind="filter"`; a SMILES that cannot be sanitized, or a target `aizynthcli` returned no row for, as `kind="failure"`.

**Example**:
```python
from biopipelines import AiZynthFinder, Panda, CompoundLibrary

library = CompoundLibrary("designs.csv")
azf = AiZynthFinder(compounds=library, max_routes=10, time_limit=120, nproc=8)

# Keep only designs that are synthesizable at all
makeable = Panda(
    tables=azf.tables.targets,
    operations=[Panda.filter("is_solved == True"),
                Panda.sort("top_score", ascending=False)],
    pool=library,
)

# The buyable building blocks, as a compounds stream ready for the next tool
buyable = Panda(
    tables=azf.tables.precursors,
    operations=[Panda.filter("in_stock == True"),
                Panda.drop_duplicates("smiles")],
    pool=azf,
)
```

Shortest-route-first selection uses the `routes` table, whose `route` column is the rank (`route == 1` is the best-scoring route for each target):

```python
best_routes = Panda(
    tables=azf.tables.routes,
    operations=[Panda.filter("route == 1 and number_of_steps <= 3")],
    pool=azf,
)
```

Choosing the search, the ranking and the constraints (all optional; every one left out keeps the config's value):

```python
azf = AiZynthFinder(
    compounds=library,
    search_algorithm="retrostar",
    max_transforms=5,
    route_scorer="route cost",                  # rank by the cheapest route; adds a route_cost column
    scorers=["fraction in stock"],              # adds a fraction_in_stock column
    filter_cutoff=0.02,                         # accept reactions the quick-filter scores low
    frozen_substructures=["O=c1ccc2ccccc2o1"],  # no route may break the coumarin
    stock_stop_criteria={"weight": 300},        # purchasable only up to 300 g/mol
)

cheapest = Panda(
    tables=azf.tables.routes,
    operations=[Panda.filter("route == 1"), Panda.sort("route_cost")],
    pool=azf,
)
```

---

### REINVENT4

Generative design of small molecules with the pretrained priors of REINVENT 4. Four generators cover the usual design questions: invent molecules from scratch (`reinvent`), decorate a scaffold with R-groups (`libinvent`), find the linker between two fragments (`linkinvent`), or propose analogues of a molecule within a chosen kind of similarity (`mol2mol`). Only REINVENT's `sampling` run mode is wrapped, so the molecules come from the prior as it is and are not optimized against any score. Rank or filter them downstream, with `samples.nll`, `RDKit`, `ADMETAI` or `AiZynthFinder`.

The output is a **compounds** stream with the canonical `id | format | code | smiles | ccd` schema, so it feeds any compounds consumer, and a `samples` table with REINVENT's negative log-likelihood for every molecule. The generators that start from input compounds run REINVENT once per input, so each molecule keeps the id of the compound it came from (`<input id>_<n>`, the input's id repeated in `compounds.id`) and an input that REINVENT rejects fails alone instead of aborting the others. The prior is loaded again for each input.

**Tags**: generate-molecule, small-molecule

**References**: https://github.com/MolecularAI/REINVENT4 (Apache-2.0) — Loeffler et al. (2024) *J. Cheminform.* 16, 20.

**Environment**: `REINVENT4` (dedicated, Python 3.11). `REINVENT4.install()` clones MolecularAI/REINVENT4 (shallow, default branch, unpinned) into `folders.repositories.REINVENT4`, runs upstream's `install.py <processor> -d none` in the env (core dependencies only: no OpenEye, Chemprop, iSIM or PUMAS), and downloads ten priors into `<repo>/priors/`, checking the md5 of each. `install.py` does not check pip's exit status, so the install is verified by `reinvent --version`, which imports PyTorch. An existing env is kept, so a download that failed halfway is retried by calling `install()` again; `force_reinstall=True` rebuilds the env and reruns `install.py` on the existing clone (delete the repository folder to refresh the clone).
- `processor`: str = `"cu126"` — install keyword: the PyTorch build tag upstream's `install.py` takes (`cu126`, `cu128`, `cpu`, `rocm6.4`, `xpu`, `mac`). `REINVENT4.install(processor="cpu")` for a machine without a GPU.

The priors come from Zenodo record 15641297 (version 1 of the public priors, the one upstream's documentation links), pinned by record number: the concept DOI 10.5281/zenodo.15641296 follows the newest version, and the newest one (June 2026) holds different files (`reinvent_pubchem.prior` instead of `reinvent.prior`, and none of the `mol2mol_*` priors). A GPU is optional: REINVENT runs on `cuda:0` when the interpreter in the env reports CUDA available and on the CPU otherwise, so `Resources(gpu="none")` gives a CPU run. Upstream states that about 8 GiB of main and GPU memory is enough for many design tasks.

**Parameters**:
- `generator`: str = `"reinvent"` — `"reinvent"` (de novo, no input), `"libinvent"` (scaffold decoration), `"linkinvent"` (linker design) or `"mol2mol"` (molecule optimization).
- `compounds`: DataStream | StandardizedOutput = None — The input compounds, read from the `smiles` column of the compounds map_table. Required for `libinvent`, `linkinvent` and `mol2mol`; refused for `reinvent`. What a SMILES has to look like depends on the generator:
  - `libinvent`: a scaffold with its attachment points marked `[*:0]`, `[*:1]`, ..., e.g. `[*:0]Cc2ccc1cncc(C[*:1])c1c2`. The R-groups at the points are generated.
  - `linkinvent`: a pair of warheads, each with one attachment point `*`, joined by `|`, e.g. `*c1ccc(F)cc1|*c1ccncc1` (REINVENT also reads `.` as the separator). The linker between them is generated.
  - `mol2mol`: any molecule. REINVENT standardizes it (largest fragment, uncharged) before sampling.

  These strings are not molecules RDKit can embed, so give them to the tool as written with `CompoundLibrary({"w1": "*c1ccc(F)cc1|*c1ccncc1"})` rather than building 3-D coordinates for them with `Ligand`. REINVENT validates the input against the vocabulary of the prior, and one it cannot read fails that input alone.
- `num_compounds`: int = 100 — Molecules sampled per input compound, or in total for `reinvent` (upstream `num_smiles`). REINVENT removes invalid SMILES and duplicates, so fewer can come back. With `sample_strategy="beamsearch"` it is the beam size.
- `similarity`: str = None — `mol2mol` only: which prior, that is which kind of similarity the analogues have to their input. Left unset it is `"medium_similarity"`, the one upstream's example configuration uses. Definitions are those of Table 4 of the REINVENT 4 paper:
  - `"similarity"`: pairs with Tanimoto >= 0.5
  - `"medium_similarity"`: 0.5 <= Tanimoto < 0.7
  - `"high_similarity"`: Tanimoto >= 0.7
  - `"scaffold"`: pairs sharing the same Murcko scaffold
  - `"scaffold_generic"`: pairs sharing the same unlabelled Murcko scaffold
  - `"mmp"`: matched molecular pairs, extracted with mmpdb
  - `"pubchem"`: trained on PubChem pairs with Tanimoto >= 0.5 on ECFP4 count fingerprints, with a ranking loss that ties the NLL to the similarity

  The first six are trained on ChEMBL.
- `sample_strategy`: str = None — `mol2mol` only: `"multinomial"` (stochastic; upstream's default) or `"beamsearch"` (deterministic, always unique).
- `temperature`: float = None — `mol2mol` only, multinomial sampling only: the factor the logits are divided by (upstream default 1.0). Below 1 is more conservative, above 1 more random. Refused together with `"beamsearch"`, which does not use it.
- `seed`: int = None — Random seed passed to REINVENT (`--seed`) so a run can be repeated; an integer >= 1, because upstream treats 0 as no seed. Unset leaves it random.

**Streams**:
- `compounds` — the generated molecules, value-based. ids are `<input id>_<n>` for `libinvent`, `linkinvent` and `mol2mol`, and `<pipeline name>_<n>` for `reinvent`, with `n` counting from 1 in the order REINVENT wrote the molecules.

**Tables**:
- `compounds` — `id | format | code | smiles | ccd | compounds.id`. The stream's map_table: `format` is `smiles`, `code` is `LIG` (as `Ligand` does for a SMILES), `ccd` is empty, `compounds.id` is the input compound the molecule came from and is absent for `reinvent`.
- `samples` — One row per generated molecule, with the columns REINVENT reports for the generator. `nll` is the negative log-likelihood of the molecule under the prior, rounded to two decimals by REINVENT: lower means the prior finds it more probable. The columns by generator:
  - `reinvent`: `id | smiles | nll`
  - `libinvent`: `id | compounds.id | smiles | nll | scaffold | r_groups`
  - `linkinvent`: `id | compounds.id | smiles | nll | warheads | linker`
  - `mol2mol`: `id | compounds.id | smiles | nll | input_smiles | tanimoto`

  `scaffold` is the scaffold as given to the model and `r_groups` the R-groups generated for it; `warheads` and `linker` likewise for `linkinvent`; `input_smiles` is the input as given to the model and `tanimoto` REINVENT's Tanimoto similarity of the analogue to it. `compounds.id` is the input compound the molecule came from. REINVENT writes isomeric SMILES for `mol2mol` only: for the other generators `isomeric_smiles` is at its upstream default (false), so the molecules carry no stereochemistry.
- `inputs` — `id | smiles | num_compounds`. `libinvent`, `linkinvent` and `mol2mol` only; one row per input compound, with the number of molecules it yielded (empty when its run failed or it had no SMILES).
- `missing` — `id | removed_by | kind | cause`. An input with no SMILES, or that yielded no valid molecule, is `kind="filter"`; an input REINVENT could not run (its message is the `cause`) is `kind="failure"`.

Every REINVENT run leaves its configuration, input file, raw CSV and log in the step's `execution/sampling/` folder, named after the input id (`denovo` for `reinvent`).

**Not wrapped**: transfer learning, staged learning (reinforcement learning with scoring components), the `scoring` run mode, the PepInvent generator, the transformer variants of LibInvent and LinkInvent, and the optional sampling filters and checks (`filter.smarts`, `target_smiles_path`, `isomeric_smiles`, `randomize_smiles`, `tb_logdir`).

**Example**:
```python
from biopipelines import REINVENT4, CompoundLibrary, Panda

# De novo: 500 molecules from the prior
novel = REINVENT4(generator="reinvent", num_compounds=500, seed=7)

# Decorate a scaffold with R-groups at both attachment points
scaffolds = CompoundLibrary({"quinoline": "[*:0]Cc2ccc1cncc(C[*:1])c1c2"})
decorated = REINVENT4(generator="libinvent", compounds=scaffolds, num_compounds=200)

# Link two warheads
warheads = CompoundLibrary({"w1": "*c1ccc(F)cc1|*c1ccncc1"})
linked = REINVENT4(generator="linkinvent", compounds=warheads, num_compounds=100)

# Analogues of a molecule that stay above Tanimoto 0.7
reference = CompoundLibrary({"celecoxib": "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"})
analogues = REINVENT4(generator="mol2mol", compounds=reference, similarity="high_similarity",
                      num_compounds=100)

# The 50 decorated molecules the prior finds most probable, with the stream carried along
likely = Panda(tables=decorated.tables.samples,
               operations=[Panda.sort("nll"), Panda.head(50)],
               pool=decorated)
```
