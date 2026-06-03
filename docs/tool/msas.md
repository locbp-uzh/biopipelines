# MSAs

[← Back to Tool Reference](../tool_reference.md)

Multiple-sequence-alignment (MSA) tools generate and convert the alignments that structure-prediction models (AlphaFold, Boltz2) consume. Most users never touch these directly — AlphaFold and Boltz2 fetch their own MSAs from a public server by default. Reach for these tools when you want to **generate an MSA once and reuse it** across predictions, or **recycle** an MSA from one tool into another.

---

### MMseqs2

Generates multiple sequence alignments for structure prediction by querying a local MMseqs2 server. The client auto-starts the server (mode from `tool_overrides.mmseqs2server.mode`, default CPU) if one is not already running, so you don't need to launch `MMseqs2Server` yourself. Feed its `msas` output to `Boltz2(msas=...)`/`AlphaFold(msas=...)` to avoid the public MSA server's rate limits when folding many sequences.

**References**: https://github.com/soedinglab/MMseqs2

**Environment**: `biopipelines`

> **Platform note**: Only the partial uniref30 database is available on HPC (migrating to `colabfold_search`); the full sequence databases are too large for Colab. For most workflows, letting AlphaFold/Boltz2 use their built-in MSA server is simpler.

**Parameters**:
- `sequences`: str | List[str] | DataStream | StandardizedOutput (required) — Input sequences.
- `output_format`: str = "csv" — Output format (`"csv"`, `"a3m"`).
- `timeout`: int = 3600 — Server timeout in seconds.
- `mask`: str | tuple = "" — Optional region of each sequence to mask out of the MSA query (PyMOL-style selection string or `(TableInfo, column)`).

**Streams**: `msas`

**Tables**:
- `msas`: | id | sequences.id | sequence | msa_file |

**Example**:
```python
from biopipelines.mmseqs2 import MMseqs2

msas = MMseqs2(sequences=lmpnn, timeout=7200)
```

---

### MMseqs2Server

Starts and manages a local MMseqs2 server process (CPU or GPU mode) so repeated MSA queries hit a warm local server instead of the public endpoint. It only manages server infrastructure — it does not process sequences itself.

**Environment**: `biopipelines`

**Parameters**:
- `mode`: str = "cpu" — Server mode (`"cpu"` or `"gpu"`).
- `database`: str = "uniref30_2302_db" — Database to use.
- `max_seqs`: int = 10000 — Maximum sequences returned per query.
- `threads`: int = None — Number of threads (auto-detect if None).
- `poll_interval`: int = 10 — Job polling interval in seconds.
- `gpus`: int = 1 — Number of GPUs for the GPU server (1 or 2). With 2, the UniRef30 and environmental gpuservers are pinned to separate GPUs so their prefilters don't share one device. Only meaningful for `mode="gpu"`.
- `idle_timeout`: int = 1800 — Seconds of inactivity before the server shuts itself down.

**Note**: This is an infrastructure helper for advanced HPC setups. Typical pipelines do not need it.

---

### MSA

Converts MSA files between CSV (Boltz2 / public-server format) and A3M (AlphaFold/ColabFold format). Enables MSA recycling between prediction tools.

**Environment**: `biopipelines`

**Parameters**:
- `msas`: StandardizedOutput (required) — Tool output with an `msas` stream (e.g. from Boltz2, AlphaFold, MMseqs2). Pass the whole tool output, not `.msas`.
- `convert`: str (required) — Target format: `"a3m"` or `"csv"`.

**Streams**: `msas`

**Tables**:
- `msas`: | id | sequences.id | sequence | msa_file |

**MSA Recycling Compatibility**:

| Direction | Works? | Notes |
|-----------|--------|-------|
| AlphaFold → A3M → CSV → Boltz2 | Yes | A3M-to-CSV conversion works; Boltz2 accepts the converted CSV for recycling. |
| Boltz2 → CSV → A3M → AlphaFold | No | ColabFold ignores the converted A3M (it lacks the original headers ColabFold expects) and re-queries the MMseqs2 server. |

**Example**:
```python
from biopipelines.msa import MSA

# Convert AlphaFold A3M MSAs to CSV so Boltz2 can recycle them
af2_result = AlphaFold(proteins=seq)
csv_msas = MSA(af2_result, convert="csv")
boltz = Boltz2(proteins=seq, ligands=lig, msas=csv_msas)
```
