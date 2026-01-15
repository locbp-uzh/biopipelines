# Ligands

This folder is for storing ligand PDB files used by pipelines.

- **Not committed**: Your ligand files here are gitignored (except this README)
- **Auto-cached**: Downloaded ligands are automatically cached here for reuse

## File Organization

```
Ligands/
├── README.md         # This file (committed)
├── ATP.pdb           # Cached ligands (not committed)
├── ATP.csv           # Metadata (not committed)
├── aspirin.pdb       # Cached ligands (not committed)
├── aspirin.csv       # Metadata (not committed)
└── ...
```

## CSV Metadata

Each cached ligand has a companion `.csv` file with metadata columns:
- `id`, `code`, `lookup`, `source`, `ccd`, `cid`, `cas`, `smiles`, `name`, `formula`
