# MSAs

[← Back to Tool Reference](../tool_reference.md)

Multiple-sequence-alignment (MSA) tools generate and convert the alignments that structure-prediction models (AlphaFold, Boltz2) consume. Most users never touch these directly — AlphaFold and Boltz2 fetch their own MSAs from a public server by default. Reach for these tools when you want to **generate an MSA once and reuse it** across predictions, or **recycle** an MSA from one tool into another.

---

### MMseqs2

Generates multiple sequence alignments for structure prediction by querying a local MMseqs2 server. The client auto-starts the server (mode from `tool_overrides.mmseqs2server.mode`, default CPU) if one is not already running, so you don't need to launch `MMseqs2Server` yourself. Feed its `msas` output to `Boltz2(msas=...)`/`AlphaFold(msas=...)` to avoid the public MSA server's rate limits when folding many sequences.

**Tags**: build-msa, fetch, msa, protein

**References**: https://github.com/soedinglab/MMseqs2

**Environment**: `biopipelines`

> **Platform note**: Only the partial uniref30 database is available on HPC (migrating to `colabfold_search`); the full sequence databases are too large for Colab. For most workflows, letting AlphaFold/Boltz2 use their built-in MSA server is simpler.

**Parameters**:
- `sequences`: str | List[str] | DataStream | StandardizedOutput (required) — Input sequences.
- `output_format`: str = "csv" — Format requested *from the search* (`"csv"`, `"a3m"`). The step always emits `<id>.csv` either way; an a3m result is converted on arrival, so this selects the wire format, not the output format. A `server_url` server returns a3m regardless.
- `timeout`: int = 3600 — Server timeout in seconds.
- `mask`: str | tuple = "" — Optional region of each sequence to mask out of the MSA query (PyMOL-style selection string or `(TableInfo, column)`).
- `server_url`: str = "" — Query a remote ColabFold-protocol MSA server over HTTP instead of starting a local one. Any `output_format` is accepted (the step converts the a3m the protocol returns). Set `MMSEQS2_SERVER_USER` / `MMSEQS2_SERVER_PASSWORD` for basic auth; credentials are refused over plain HTTP. This queries a **shared public service** when pointed at `api.colabfold.com` — see the note below before running a large batch through it.

**Using the public ColabFold server.** `server_url="https://api.colabfold.com"` queries a **free community service** run by the ColabFold authors, not lab infrastructure. The client identifies itself, checks the submit response for `RATELIMIT` / `MAINTENANCE` instead of polling a ticket that was never queued, and backs off with jitter as it waits. None of that makes a large batch appropriate: many pipeline steps running concurrently each poll independently, and there is no cross-step throttle. For anything beyond a handful of sequences, run a local `MMseqs2Server` — that is what it exists for. If a batch fails entirely the step now exits non-zero rather than writing an empty MSA table, which previously let every downstream fold run single-sequence with nothing saying so.

**Streams**: `msas`

**Tables**:
- `msas`: | id | sequences.id | sequence | file |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.mmseqs2 import MMseqs2

msas = MMseqs2(sequences=lmpnn, timeout=7200)
```

---

### MMseqs2Server

Starts and manages a local MMseqs2 server process (CPU or GPU mode) so repeated MSA queries hit a warm local server instead of the public endpoint. It only manages server infrastructure — it does not process sequences itself.

**Tags**: data, build-msa, msa, protein

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

**Tags**: data, msa

**Environment**: `biopipelines`

**Parameters**:
- `msas`: StandardizedOutput (required) — Tool output with an `msas` stream (e.g. from Boltz2, AlphaFold, MMseqs2). Pass the whole tool output, not `.msas`.
- `convert`: str (required) — Target format: `"a3m"` or `"csv"`.

**Streams**: `msas`

**Tables**:
- `msas`: | id | sequences.id | sequence | file |

**MSA Recycling Compatibility**:

| Direction | Works? | Notes |
|-----------|--------|-------|
| AlphaFold → A3M → CSV → Boltz2 | Yes | A3M-to-CSV conversion works; Boltz2 accepts the converted CSV for recycling. |
| AlphaFold → A3M → AlphaFold | Yes | No conversion involved — pass the producing tool straight to `msas=`. |
| MMseqs2 → A3M → AlphaFold | Yes | Use `MMseqs2(..., output_format="a3m")`, which emits the search result verbatim. No `MSA()` step, and the UniRef headers survive. **Prefer this** over any route through CSV. |
| Boltz2 → CSV → A3M → AlphaFold | Monomers only | The CSV form is `key, sequence` with no headers, so `MSA(convert="a3m")` has to invent `>100, >101, …`. ColabFold's `parse_fasta` ignores the `#` line and reads only the `>` text, so a monomer folds fine. A **complex** does not: chain pairing reads species off UniRef headers, and those are gone. |

Headers are the whole difference between these rows. Once an MSA has been through the CSV form its identifiers are unrecoverable, so take A3M out of the producer when the consumer is AlphaFold.

**Example**:
```python
from biopipelines.msa import MSA

# Convert AlphaFold A3M MSAs to CSV so Boltz2 can recycle them
af2_result = AlphaFold(proteins=seq)
csv_msas = MSA(af2_result, convert="csv")
boltz = Boltz2(proteins=seq, ligands=lig, msas=csv_msas)
```
