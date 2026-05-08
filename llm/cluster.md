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

# Resume after a cancel or failure.
llm/log.sh ssh cluster "cd ~/biopipelines && ./resubmit <RunTime>/slurm.sh"

# Inspect the slurm log.
llm/log.sh ssh cluster "tail -n 100 <RunTime>/slurm.out"

# Pull a result artifact back locally.
llm/log.sh scp cluster:<RunTime>/slurm.out ./

# Cancel a job.
llm/log.sh ssh cluster "scancel <jobid>"
```

## Logging

Each `log.sh` invocation appends:

- a timestamped header with the (shell-quoted) argv,
- the combined stdout/stderr of the command,
- an `=== exit=N ===` footer,

to `llm/logs/YYYY-MM-DD.log`. The `logs/` folder is gitignored.
