# Multi-chain designs in a `sequences` stream

Decision record for the change that split a multi-chain inverse-folding record into one `sequences` row per chain, and added the `Grouped` axis that puts them back together. Written 2026-09-21.

## The failure

Reported from a binder-design campaign: four 100-design runs produced high-confidence complexes of the wrong molecule, and nothing warned.

On a two-chain backbone ProteinMPNN emits both chains in a single `sequence` value joined by `/` — its own native convention. The `sequences` stream passed that string through unchanged, and `Boltz2` wrote it verbatim into one YAML chain entry. Boltz mapped the unrecognized `/` to a `UNK` residue and folded on, producing a chain containing the designed sequence, a `UNK`, and a second fused copy of the fixed chain. Reported `iptm` was 0.94–0.97 for that malformed complex. The only tool that objected was `Prodigy`, which validates residue names and rejected all 393 designs with `Unsupported non-standard amino acid found: UNK` — and the step still reported COMPLETED, which is the second defect below (§4).

## Why this is a contract violation, not a missing feature

The framework already expresses a complex as **several rows combined by an axis**: `Bundle` for rows that always travel together. AlphaFold colon-joins a bundle's chains for ColabFold; Boltz2 allocates one chain id per bundled row through a shared counter. Both read the *axis*, never the cell.

A separator inside one cell bypasses that entirely, and nothing in `biopipelines/` or `pipe_scripts/` has ever split on one — a grep for a chain-separator split returns only path handling. So a fused value reaches all 33 sequence-consuming pipe scripts as if it were a single chain. Boltz2 is where it was caught, not where it stops: an MSA built from a `/`-bearing string is garbage, and `get_msa_file`'s by-sequence lookup could never match a per-chain MSA.

The rule is therefore: **one `sequences` row is one polymer chain.** It is recorded in `developer_manual.md` as the Chain Contract.

## Options considered

| Option | Why not |
|---|---|
| Declare `/` a first-class separator and teach every consumer to split | Creates a second multi-chain vocabulary beside `Bundle`; requires touching ~33 consumers; silently disables Boltz MSA recycling, whose sequence key would be the joined string |
| One id, several map_table rows, disambiguated by a `chain` column | `select_ids`, `write_filtered_map_table`, `get_value` and Panda's slicers all assume id→row is 1:1. The framework's existing "one id, many sub-entities" precedent (`resi-csv`) puts the sub-entities in a per-id *file*, not in the map_table |
| Refuse only — error on a separator, no correct pipeline supported | Removes the silence but supports no multi-chain work |
| **Split at the producer, regroup with an axis** | **Chosen.** Makes all 33 single-sequence consumers correct with no change to any of them, keeps MSA recycling working, and makes "fold only the designed chain" expressible |

The chosen option's one real cost: with one row per chain, a `Panda` filter on a per-design metric can keep one chain and drop another, splitting a complex. Group-aware filtering is not implemented; see *Not done*.

## What was built

### 1. The producer splits — `pipe_fa_to_csv_fasta.py`

Shared by **both** MPNN wrappers, so one fix covers ProteinMPNN (`/`) and LigandMPNN (`:`).

A new `chains=` parameter says which chains become rows. It follows `PDB.split_chains`'s existing convention (`pdb.py:988-997`) rather than inventing one: an explicit chain list gives deterministic ids, `"all"` gives the lazy bracket.

| `chains=` | `sequences` ids | Lazy? |
|---|---|---|
| `None` (default) | `<structure>_<n>` — today's ids | no |
| `"A"` / `["A"]` | `<structure>_<n>` | no |
| `["A","B"]` | `<structure>_<n>_<chain>` | no |
| `"all"` | `<structure>_<n>[_<?>]` | yes |

Only a setting that can emit more than one chain adds a suffix, so **every existing single-chain pipeline is byte-identical**. Laziness is opt-in and confined to the one spelling whose chain letters genuinely cannot be known at configuration time.

Under the default, a multi-chain backbone is recorded in `missing.csv` with `kind=failure` naming the chains found. A local failure row is not excused by the completion check, so the step is reported FAILED — loud and per-design, rather than aborting a whole batch on one bad structure.

**Chain letters** are resolved by matching each segment's length against the backbone's per-chain residue counts, trying alphabetical order (what `protein_mpnn_run.py` concatenates in) and then the PDB's own record order. Counts rather than sequences, because the model's three-to-one letter mapping and `pdb_parser`'s disagree on non-standard residues while the counts do not. A design that cannot be reconciled is dropped with a cause naming both, never mislabelled.

Deduplication moved from the chain to the **design**: deduplicating per chain would delete the second chain of a homodimer.

### 2. `designs` — the grouping key

Both wrappers now emit a `designs` stream carrying one row per design, with **exactly the ids `sequences` used to have**. It is deterministic whatever `chains=` does, which is what keeps a grouped consumer's output ids deterministic even when the member stream is lazy. `designs` was already in `KNOWN_STREAM_NAMES`, so no registry change was needed.

### 3. `Grouped` — the axis

`Each` makes one output per row; `Bundle` makes one output from every row. Neither can say "these rows belong in one prediction". `Grouped` iterates the **groups** of a stream, taking its partition from a stream of parents exactly as `Consensus`'s `groups=` does.

```python
folded = Boltz2(proteins=Bundle(Grouped(pmpnn), tag))
```

`groups=` takes any stream whose ids are the keys, at any ancestor level — the chain rows' design, or the backbone they all came from. Omitting it defaults to the source's `designs`, reached either from a tool output's own streams or, for a bare stream, through its `_producer` back-reference. Only a stream built by hand, with no producer, raises.

Two properties, both verified:

- **Cardinality is a configuration-time fact, membership a runtime one.** Output ids are the group ids, which come from the group stream.
- **Membership is the framework's id matching**, `get_mapped_ids(group_ids, member_ids, unique=False)`, using the child-match tier. This is why ids as unrelated-looking as `chihuahua+labradudor_1_A` still find their design: `+` is an ordinary character inside a segment, and only `_` separates parent from child. Do **not** pass `closest_siblings_only=True` — it restricts to the sibling tier and returns nothing here.

Runtime support is in `pipe_boltz_config_unified.py` (a group becomes consecutive chain entries off the shared chain-id counter) and `pipe_alphafold_queries.py` (a group's members are colon-joined, and a static partner then appends one more chain).

### 4. The completion check

One defect, and it is what made the corruption above easy to miss.

- **A step that exited non-zero could be marked COMPLETED.** The check judged declared paths by existence alone and wrote its marker *before* the footer tested `BP_MAIN_RC`; `pipeline.py`'s guard then found a COMPLETED marker and declined to write FAILED. The exit status is now passed in with `--main-rc` and decides first.

## Verified

Producer, both consumers and the completion check are covered by `tests/test_multichain_sequences.py` and `tests/test_completion_exit_status.py`. The Boltz2 test asserts the exact shape the report showed as malformed: three separate chain entries, no separator in any of them, and the partner present once rather than twice.

Not verified: nothing has run on a GPU or a cluster. This is a local, configuration-time and pipe-script-level change set.

## Not done

- **Group-aware filtering in Panda.** A filter on a per-design metric can keep one chain of a complex and drop another. This is the chosen option's one structural cost and it is open.
- **Chain provenance back to a bundled origin.** For a structure folded from `Bundle(chihuahua, labradudor)` the mapping chain A ↔ chihuahua exists upstream (`proteins.1` / `proteins.2` on the folded structure's map_table, chain letters allocated in item order). ProteinMPNN deliberately records only the PDB chain letter and leaves that join to the user: the correspondence is order-dependent, `static_first` is computed at runtime, and the developer manual already warns against decomposing a `+` id by position.
- **A monomer's `chain` column is empty** under the default, because the fast path does not parse the PDB. Filling it would mean parsing every backbone for a label nobody asked for.

## RFdiffusion3

RFdiffusion3 splits its `sequences` rows per chain and carries the `chain` column, so it satisfies the row contract. It does **not** emit a `designs` stream, so the rows are not groupable by the default axis: `Grouped(rfd3)` has nothing to partition by, and refolding a multi-chain design as one complex needs the partition named explicitly — `Grouped(rfd3.streams.sequences, groups=rfd3.streams.structures)`.

Its declared ids follow the contig rather than a `chains` parameter, since the chain count is a property of the contig: a `/0` chain break declares the lazy `<structure>[_<?>]`, anything else declares one id per structure. A contig arriving per structure at runtime cannot be read at configuration time and is treated as multi-chain, because declaring too few ids fails silently and declaring a lazy one costs only precision.

Giving RFdiffusion3 a `designs` stream would complete the contract and is not done.
