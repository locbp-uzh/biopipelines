# Cluster usage

**Check your tool list before reading any further.** If it holds `bp_submit`, `bp_status`, `bp_logs`, `bp_runs`, `bp_table`, `bp_lineage` and `bp_fetch`, the `bp-mcp` server is registered and *those* are how you reach the cluster. The `ssh` and `scp` idioms in this file are the fallback for a session without them — running one anyway means losing the operations record the tool writes, guessing which scheduler file exists, and re-deriving status the tool already reports.

| What you want | With `bp-mcp` | Without it |
| --- | --- | --- |
| First connection | `bp_setup(host="s3it", repo="~/biopipelines")` | the ssh alias below, by hand |
| Submit a pipeline | `bp_submit(script="my_pipelines/foo.py", upload="my_pipelines/foo.py")` | `ssh s3it "./submit ..."` |
| Run status | `bp_status(job="<Job>_NNN")` | `ls` markers + `squeue` |
| One step's log | `bp_logs(job=..., step="007_Boltz2")` | `ssh s3it "tail Logs/..."` |
| A result table | `bp_table(job=..., step=..., table=...)` | `scp` the CSV |
| Ids produced and dropped | `bp_lineage(job=...)` | read `*_map.csv` / `missing.csv` |
| List projects and runs | `bp_runs()` | `ls` the output root |
| Pull a file back | `bp_fetch(remote_path=...)` | `scp` |
| Resume after a failure | `bp_resubmit(job=...)` | `ssh s3it "./resubmit ..."` |
| Stop a run | `bp_cancel(job=...)`, then `confirm=True` | `ssh s3it "scancel <ids>"` |
| Show a step's results | `bp_visualize(job=..., step=...)` | `bp-visualize` then `scp` the page |

The MCP tools already know where the output root is: they ask the cluster, because the config resolves `<username>` on the machine that owns it. You do not need to look it up or pass it.

`bp_resubmit`, `bp_cancel` and `bp_visualize` cover resuming, stopping and showing a run, so the only things left for a shell are `sinfo` and `squeue`. Say which tool you checked first.

## One-time setup

1. **Add an ssh alias** for your cluster in `~/.ssh/config`. **Name it after the site, not after what kind of machine it is.** `cluster` reads fine until you have a second one — and Daint is also a cluster, so the name stops telling you which machine you are on. Site names stay true:

   ```
   Host s3it
       HostName cluster.s3it.uzh.ch
       User <your-username>
       IdentityFile ~/.ssh/id_ed25519

   Host daint
       HostName daint.alps.cscs.ch
       User <your-username>
   ```

   Verify with `ssh s3it echo ok`.

   The alias is a **place**; the config variant is a **shape**. They are different namespaces and should not be expected to match: S3IT and any other SLURM/conda site both use `variant: cluster`, while Daint has its own variant because aarch64, the Container Engine and per-project billing genuinely differ — not because it is a particular machine. Renaming an alias later costs one line here plus one `bp_setup`; `Host s3it cluster` keeps the old name working while you switch.

   > **Safety note.** The snippet above is safe to share. `HostName`, `User`,
   > and the *path* to your key file are not secrets — `~/.ssh/config` is
   > meant to be readable by you and is fine to keep in a dotfiles repo.
   > The secret is the **private key file itself** (`~/.ssh/id_ed25519`):
   > keep it at permissions `600` (`chmod 600 ~/.ssh/id_ed25519`), protect
   > it with a passphrase, never commit or paste its contents, and only
   > share the matching `.pub` file. Do **not** silence host-key checks
   > (`StrictHostKeyChecking no`) or enable `ForwardAgent yes` to hosts you
   > don't fully trust — both weaken the protections this setup relies on.

2. **Clone the repo on the cluster** at a stable path, e.g.
   `/home/<user>/biopipelines`.

## Using the cluster without the MCP server

Everything below applies only when the `bp_*` tools are absent from your tool list. Run any command (ssh, scp, etc.) directly:

```bash
ssh s3it 'squeue -u $(whoami)'
```

Common idioms (substitute your actual remote repo path):

```bash
# Sync the remote checkout to a branch (destructive — confirm first).
ssh s3it "cd ~/biopipelines && git fetch && git reset --hard origin/<branch>"

# Copy a personal pipeline or input file to the cluster.
scp my_pipelines/foo.py s3it:~/biopipelines/my_pipelines/

# Submit a pipeline.
ssh s3it "cd ~/biopipelines && ./submit my_pipelines/foo.py"

# See which job scripts and scheduler outputs the run produced.
ssh s3it "ls <RunTime>"

# Resume after a cancel or failure (single-batch: slurm.sh; multi-batch: slurm_batch<N>.sh).
# Dependency directives are stripped by default -- the script on disk names the
# original run's job ids. Add --keep-dependencies only if the parent is still queued.
ssh s3it "cd ~/biopipelines && ./resubmit <RunTime>/slurm_batch1.sh"

# Inspect a scheduler log (single-batch: slurm.out; multi-batch: job_batch<N>.out).
ssh s3it "tail -n 100 <RunTime>/job_batch1.out"

# Render one step's outputs as a self-contained page, then pull it and open it.
# Runs on the login node; needs no scheduler, works while the job is still going.
ssh s3it "cd ~/biopipelines && bp-visualize <Job>/<NNN>_<Tool> --descending <table>.<column> --max-items 5"
scp s3it:<Job>/<NNN>_<Tool>/_extras/<NNN>_<Tool>_view.html ./

# Inspect a single tool's log.
ssh s3it "tail -n 100 <Job>/Logs/<NNN>_<ToolName>.log"

# Cancel a job.
ssh s3it "scancel <jobid>"
```

**The file names depend on how many batches the pipeline has.** Every `Resources()` call opens a new batch and each batch is submitted as its own job, so `submit` writes `<RunTime>/slurm_batch<N>.sh` → `<RunTime>/job_batch<N>.out` per batch. A pipeline that ended up with exactly one batch gets `<RunTime>/slurm.sh` → `<RunTime>/slurm.out` instead. `ls <RunTime>` first; do not assume `slurm.out` exists.

`<RunTime>` is `<Job>_NNN/RunTime`, alongside `<Job>_NNN/Logs` (per-tool `<NNN>_<ToolName>.log` files) and the `<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING` completion markers, which sit directly in `<Job>_NNN/`.

## Logging

A run's operational log lives **with the run**, not in a dated session file: `<Job>_NNN/_operations.jsonl`, alongside `Logs/` and the completion markers. It records what was done to that run — submitted, resubmitted, cancelled, fetched — so the record travels with the results when they are copied off the cluster, and exists identically on Colab and container backends.

The older `llm/log.sh` wrapper, which tee'd every `ssh`/`scp` to `llm/logs/YYYY-MM-DD.log`, has been retired: it logged the transport rather than the work, was anchored to the date rather than to the run, covered only the cluster backend, and had no mechanism behind it. Existing files under `llm/logs/` are historical and untracked; nothing writes there any more.
