# Per-ID ancestry — which design came from which structure

`bp_lineage` answers *how many* IDs each step produced and dropped. It does not answer *which design came from which input*, and that is the question a referee, a supervisor and a pharma partner actually ask. This document settles how the answer is assembled.

## What already exists on disk

Every step writes `<stream>_map.csv` with an `id` column and, where the framework knows it, one or more `<axis>.id` provenance columns. A real eight-step campaign (`BinderDesign/tagged_binder_001`) shows all three shapes the reader has to cope with:

| Step | `id` | provenance columns | shape |
|---|---|---|---|
| `003_RFdiffusion` | `3kzy_A_1` | `structures.id = 3kzy_A` | concrete parent |
| `004_ProteinMPNN` | `3kzy_A_1_1` | `structures.id = 3kzy_A_1` | concrete parent |
| `005_Boltz2` | `3kzy_A_100_1+tag` | `proteins.id = 3kzy_A_<1..100>_<1..2>` | **a pattern, not a parent** |
| `009_Panda` | `9_Panda_1` | `structures.id = 3kzy_A_50_2+tag`, `pool.id` | concrete parent |

A fourth shape exists and is invisible here: `prune_redundant_provenance_columns` deliberately **drops** a plain `<axis>.id` column when every cell is recoverable from the id alone. So a missing provenance column is not missing lineage — it is lineage the id itself carries.

This is why ancestry cannot be a CSV join. Two of the four shapes have no concrete parent to join on.

## The rule

For each produced id, parents are resolved by descending tiers, and **the tier is recorded on the edge**. A reader who cannot see how an edge was established cannot audit it.

1. **`recorded`** — a `<axis>.id` cell holds a concrete id that some earlier step produced. Taken as-is; nothing is inferred.
2. **`matched`** — no concrete cell resolved, so the child id is matched against the ids earlier steps produced, using `id_map_utils.get_mapped_ids_with_tiers`. The edge carries the matcher's own tier (`parent`, `sibling`, …). This is the same ladder the framework uses at runtime to wire one step's output into the next step's input, so an edge found this way is the wiring that actually happened, not a guess about it.
3. **unresolved** — neither tier answered. Recorded as a root with no parent, and counted, because a campaign where many ids have no traceable parent is a finding, not an empty result.

Provenance cells containing `<` are patterns and never satisfy tier 1. They are combinatorial declarations (`3kzy_A_<1..100>_<1..2>`), not identities.

## Ancestry is a DAG, not a tree

A Boltz2 complex has two parents — the designed chain and the tag — and the `+` convention is what lets the matcher find both:

```
get_mapped_ids_with_tiers(['3kzy_A_100_1+tag'], ['3kzy_A_100_1'])  ->  ('3kzy_A_100_1', 'parent')
get_mapped_ids_with_tiers(['3kzy_A_100_1+tag'], ['tag'])           ->  ('tag',          'parent')
```

Every id therefore has zero or more parents, and the walk-back is breadth-first over a set, not a chain. A renderer that assumes one parent silently drops the tag axis, which is exactly the half of a binder campaign a reviewer wants to check.

## Cost

Ancestry needs the `id` and `<axis>.id` columns of every map table — `bp_lineage` needs only line counts. Reading whole map tables over ssh means megabytes of file-path columns for a large campaign.

So the projection happens **on the far side**: one batched call returns only the id-ish columns of every map table under the run. Two constraints on that snippet, both learned the hard way:

- **Stdlib only, no `biopipelines` import.** A helper that imports the framework makes the reader demand that the cluster run the laptop's version; that is the `ImportError` against 1.4.0 that the MCP server already hit once.
- **`csv`, not `awk` or `cut`.** Column order varies per table and quoted cells can contain commas. Field-position arithmetic is wrong on the first Panda table that carries a `value` column.

## What it does not do

It does not reconstruct ancestry for runs written before the map tables carried provenance columns. Where the column is absent *and* the id does not carry its parent, tier 3 applies and the id is reported as a root. Saying so is the point; inventing an edge would make the artifact worthless for the one use it has.
