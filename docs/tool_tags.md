# Tool tags

[← Back to Tool Reference](tool_reference.md)

---

Every tool section in `docs/tool/*.md` carries a `**Tags**:` line. Tags exist for one job: **finding a tool when you do not know its name.** A search for "covalent" over the tool index returns nothing without them, because a capability like covalent linkage lives in a parameter, never in the one-line summary.

Tags are a *closed* vocabulary of 38 terms in four facets. Closed on purpose: an open list drifts into `dock` / `docking` / `ligand-docking` within a month, and then neither a person nor an agent can guess which spelling to search. A tag earns its place only if it changes the answer to at least two distinct questions.

## How to search

```
bp_tools(tags=["dock", "covalent"])                      # tools that do both
bp_tools(tags=["protein", "binding"], exclude=["small-molecule"])   # → Prodigy
```

Within a facet the terms are **OR**-ed; across facets they are **AND**-ed. So `tags=["dock", "predict-structure", "covalent"]` means *(docking or structure prediction) and covalent*, which is almost always what you meant. `exclude=` removes any tool carrying any excluded tag, and it is what makes the near-universal tags useful: `protein` alone matches three quarters of the catalog, but `protein` minus `small-molecule` is a sharp question.

**Tags are the intent layer, not the metric layer.** There is deliberately no `rmsd`, `plddt` or `confidence` tag — those are answered exactly, from the tools' own source, by

```
bp_tools(outputs="rmsd")     # any stream, table or column name, case-insensitive substring
bp_tools(inputs="structures")
```

Use `outputs=` when you know what the number is called and `tags=` when you only know what you want to achieve.

## ACTION — what the tool does

Every tool carries one or two.

| Tag | Applies when | Not when |
|---|---|---|
| `generate-backbone` | Invents coordinates that did not exist | Predicting a structure for a known sequence |
| `inverse-folding` | Structure in, sequence out | Constructing a sequence without a backbone |
| `design-sequence` | Builds or enumerates sequences directly | Reading a sequence off a structure |
| `predict-structure` | Sequence in, coordinates out — including co-folding a complex | Placing an existing ligand into an existing pocket |
| `dock` | Places a ligand into a **receptor structure you supply** | The complex is generated from sequence; that is co-folding, and it gets `predict-structure` + `complex` |
| `sample-ensemble` | Produces many conformers of one entity | Producing many *designs* |
| `refine-structure` | Minimises, protonates or repacks without changing identity | Changing the sequence |
| `predict-property` | A **learned or empirical model** returns a number | The value is a deterministic function of the coordinates |
| `measure` | A **deterministic function** of coordinates or chemistry returns a number | A model was trained to produce it |
| `validate` | Answers pass/fail on physical plausibility | Returning a score to rank by |
| `detect` | Finds a site that was not given | Measuring proximity to a site you already know |
| `build-msa` | Builds or serves a multiple sequence alignment | — |
| `fetch` | Reaches the **network at runtime** to bring data in | Weights are downloaded once at install time — that is nearly universal and says nothing |
| `data` | Reshapes ids, tables, selections or files between stages | — |
| `visualize` | Renders an image or a plot | — |

`measure` versus `predict-property` is deterministic-versus-learned, not regression-versus-classification. You ask for a property; how the tool arrives at it is a separate question, and method-shaped tags (`classify`, `gnn`) have never once been the way anyone searched.

`fetch` doubles as the offline filter. `exclude=["fetch"]` answers "what still runs on an isolated compute node". AlphaFold and Boltz2 are deliberately **not** `fetch`: querying an MSA server is incidental to what they do, and both accept a supplied MSA instead.

## SUBJECT — what it operates on

| Tag | Applies when |
|---|---|
| `protein` | Takes or emits a protein |
| `small-molecule` | Takes or emits a ligand or compound |
| `nucleic-acid` | Handles DNA or RNA as a modelled entity |
| `complex` | Two or more entities, or the interface *is* the object |
| `pocket` | The pocket or binding site is the object — detected, selected, or designed |
| `residues` | Operates on, selects or reports per-residue positions |
| `msa` | Consumes or produces alignments |
| `ensemble` | Many conformers of one entity |

`residues` is one tag on purpose. Someone who wants to work at residue level should not have to search `select-residues`, then `per-residue`, then `residue-scores`.

There is no `library` tag: every tool here works on a stream of many ids, which is the framework itself. A tag true of everything by construction cannot discriminate.

## READOUT — what number comes back

`binding` · `stability` · `solubility` · `fitness` · `sasa` · `interactions` · `flexibility` · `energy`

A READOUT tag is carried only where it unifies **differently named** columns. `binding` covers `affinity_pred_value`, `cnn_affinity`, `pkd_pred` and a ΔG in kcal/mol, which no literal search can reach at once.

> **`binding` does not mean the numbers are comparable.** The tools it groups answer the same *question* on incompatible *scales*, and two of them run in opposite directions: Boltz2's `affinity_pred_value` is log10(IC50), where **lower is stronger**, while GEMS's `pkd_pred` is a pKd, where **higher is stronger**. Others are uncalibrated rankings — RTMScore is a pose-scoring function, not an affinity regressor. Read the tool's `**Units.**` note before comparing or ranking across tools; `bp_tools(outputs="affinity")` prints it alongside each match.

## CAPABILITY — what it can do that the headline never says

`covalent` · `symmetry` · `motif-scaffolding` · `binder-design` · `sequence-only` · `all-atom` · `sidechain`

This is the facet the summaries lose, and the reason tags exist at all. `sidechain` means *changes or optimises sidechains without moving the backbone*.

**A capability tag describes the wrapper, not the upstream project.** If upstream can do something our wrapper does not expose, the tag is off: promising a capability you cannot reach from BioPipelines turns an honest miss into a confident wrong answer. Gnina upstream has `--covalent_rec_atom` and ProteinMPNN has `--symmetry_residues`; neither is reachable here, so neither carries the tag. Those gaps are tracked as feature requests, not papered over with a tag.

## Editing the tags

The `**Tags**:` line in the tool's own section is the single source of truth — there is no separate table to keep in sync. `biopipelines/tool_tags.py` holds the vocabulary and the lint; `tests/test_tool_tags.py` fails if a section declares a tag outside the vocabulary, or if a section's body argues for a capability tag it does not declare. Adding a *new* tag means editing `VOCABULARY` and saying, in the commit, which two questions it newly answers.

`Scripting` carries no tags. It is the escape hatch — it can do anything, so no tag helps anyone find it.
