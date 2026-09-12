# Cluster usage

The LLM workflow uses plain `ssh` and `scp` to talk to the cluster. The
optional `llm/log.sh` wrapper records every invocation to a dated log.

## One-time setup

1. **Add an ssh alias** for your cluster in `~/.ssh/config`. Example for UZH S3IT:

   ```
   Host cluster
       HostName cluster.s3it.uzh.ch
       User <your-username>
       IdentityFile ~/.ssh/id_ed25519
   ```

   Verify with `ssh cluster echo ok`.

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

## Using the cluster

Run any command (ssh, scp, etc.) directly. For the LLM-driven workflow,
prefix calls with `llm/log.sh` so the command and its output are appended to
`llm/logs/YYYY-MM-DD.log`:

```bash
# Equivalent — second form is logged.
ssh cluster 'squeue -u $(whoami)'
llm/log.sh ssh cluster 'squeue -u $(whoami)'
```

Common idioms (substitute your actual remote repo path):

```bash
# Sync the remote checkout to a branch (destructive — confirm first).
llm/log.sh ssh cluster "cd ~/biopipelines && git fetch && git reset --hard origin/<branch>"

# Copy a personal pipeline or input file to the cluster.
llm/log.sh scp my_pipelines/foo.py cluster:~/biopipelines/my_pipelines/

# Submit a pipeline.
llm/log.sh ssh cluster "cd ~/biopipelines && ./submit my_pipelines/foo.py"

# See which job scripts and scheduler outputs the run produced.
llm/log.sh ssh cluster "ls <RunTime>"

# Resume after a cancel or failure (single-batch: slurm.sh; multi-batch: slurm_batch<N>.sh).
# Dependency directives are stripped by default -- the script on disk names the
# original run's job ids. Add --keep-dependencies only if the parent is still queued.
llm/log.sh ssh cluster "cd ~/biopipelines && ./resubmit <RunTime>/slurm_batch1.sh"

# Inspect a scheduler log (single-batch: slurm.out; multi-batch: job_batch<N>.out).
llm/log.sh ssh cluster "tail -n 100 <RunTime>/job_batch1.out"

# Render one step's outputs as a self-contained page, then pull it and open it.
# Runs on the login node; needs no scheduler, works while the job is still going.
llm/log.sh ssh cluster "cd ~/biopipelines && bp-visualize <Job>/<NNN>_<Tool> --descending <table>.<column> --max-items 5"
llm/log.sh scp cluster:<Job>/<NNN>_<Tool>/_extras/<NNN>_<Tool>_view.html ./

# Inspect a single tool's log.
llm/log.sh ssh cluster "tail -n 100 <Job>/Logs/<NNN>_<ToolName>.log"

# Cancel a job.
llm/log.sh ssh cluster "scancel <jobid>"
```

**The file names depend on how many batches the pipeline has.** Every `Resources()` call opens a new batch and each batch is submitted as its own job, so `submit` writes `<RunTime>/slurm_batch<N>.sh` → `<RunTime>/job_batch<N>.out` per batch. A pipeline that ended up with exactly one batch gets `<RunTime>/slurm.sh` → `<RunTime>/slurm.out` instead. `ls <RunTime>` first; do not assume `slurm.out` exists.

`<RunTime>` is `<Job>_NNN/RunTime`, alongside `<Job>_NNN/Logs` (per-tool `<NNN>_<ToolName>.log` files) and the `<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING` completion markers, which sit directly in `<Job>_NNN/`.

## Logging

Each `log.sh` invocation appends:

- a timestamped header with the (shell-quoted) argv,
- the combined stdout/stderr of the command,
- an `=== exit=N ===` footer,

to `llm/logs/YYYY-MM-DD.log`. The `logs/` folder is gitignored.
