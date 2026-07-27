# The CSCS Alps / Daint backend (`daint` variant)

Daint is a SLURM cluster like the `cluster` variant, but three of its properties break the assumptions that variant makes:

- **aarch64 (GH200).** x86-64 conda builds and `.sif` images do not apply.
- **No lmod modules, no apptainer.** `module avail` is empty; the container runtime is the CSCS Container Engine, not apptainer.
- **Per-project billing.** A job without an account is rejected.

`config.daint.yaml` answers each. Select it with `Pipeline(config="daint")` or `BIOPIPELINES_CONFIG_VARIANT=daint`.

## venv instead of conda

`machine.env_manager.name: "venv"`. Each tool gets a venv under `machine.env_manager.venv_root`, activated by path (`source <root>/<env>/bin/activate`) rather than by name — venv has no registry, so the root is what turns an `environments:` entry like `Boltz2: "Boltz2Env"` into a path.

This is the CSCS-recommended way to install Python on Alps, and it sidesteps the aarch64 conda problem for any tool whose dependencies are pip-installable.

An env yaml pinning anything beyond `python`/`pip` **raises** at install: a venv cannot install conda packages, and silently dropping them would produce an env that looks built and fails at runtime. Tools needing real conda binaries (`mkdssp`, PyMOL, mmseqs2) need those supplied by their EDF image instead.

## GPU tools: EDF + the Container Engine

The `edf:` map points a tool at an Environment Definition File; the tool's whole script then runs via `srun --environment=<edf>`:

```yaml
edf:
  Boltz2: "<group>/edf/ngc-pytorch.toml"
```

The EDF itself is TOML, and the image can come straight from NVIDIA NGC — which publishes GH200-native PyTorch, so torch and CUDA arrive prebuilt rather than being ported:

```toml
image = "nvcr.io#nvidia/pytorch:25.06-py3"
mounts = ["/capstor", "/iopsstor"]
workdir = "/capstor/scratch/cscs/<user>"
[env]
CUDA_CACHE_DISABLE = "1"
```

A tool with no `edf:` entry runs directly on the batch node — CPU-only tools should be left out so they never pay the image pull.

Two things follow from the Container Engine defining a *job step* rather than wrapping a command:

- The wrap goes around the **whole tool script**, not the inner python call. Wrapping only the command would leave the venv activation outside the container it applies to. This is why `edf:` is separate from `containers:`, which is apptainer's command-prefix mechanism and stays empty on Daint.
- A tool's venv binds to its image's python, so it **must be built inside the same container**. Install steps are wrapped through their parent tool automatically.

Reference EDFs on shared storage, never a bare name: `--environment=ngc-pytorch` resolves against each user's `~/.edf`, so a bare name only works for whoever wrote it.

## Partitions

`Resources(partition=...)` selects one. `normal` is the production default and has a deep queue; `debug` (24 nodes, max 2 per job, 30 minutes) starts immediately, so bring-up tests belong there. `xfer` is for internal store↔scratch transfers — 2 shared nodes, no GPUs, 250 GB RAM. `low` lets some projects run past their quarterly allocation, but must be enabled for the project.

Nodes are not shared in `normal`, `debug` or `low`, so a job of any size is billed for whole nodes; `xfer` nodes are shared.

Priorities are dynamic and job shape does not change them, so expect to queue in `normal`. A usage-based project shows `RawShares=1` in `sshare`, which looks like an exhausted allocation but is not — there is simply no quarterly grant to divide. Check the real budget at https://account.cscs.ch, not on the cluster.

## Billing account

`machine.billing_account: "uzh67"` → `#SBATCH --account=uzh67`. This is the *project* charged for compute, not the user (`machine.username`). It is emitted in the batch header only — `srun` job steps inherit the allocation's account. Empty on every other variant, so their scripts are unchanged.

## Storage

The project's `$STORE` quota is 1.0 TB but only **150,000 inodes**, and the inode limit is what binds. Python environments are inode-heavy and nearly weightless: 13 tools' venvs came to 101,223 files for 15 GB, and adding PyMOL failed with `Disk quota exceeded` while the filesystem was 98.6% empty by size.

| Config | Path | Why |
|---|---|---|
| `venv_root` | `<scratch>` | 1M inodes; `$STORE` cannot hold Python envs |
| `edf:`, `cache`, `python_headers` | `<group>` = `/capstor/store/cscs/uzh/uzh67` | few files, backed up daily, setgid `uzh67` |
| `biopipelines_output` | `<scratch>` | fast, large, and run outputs are transient |

The cost of envs on scratch is the 30-day access-time purge — a tool unused for a month needs reinstalling. Installs are idempotent, so this is cheaper than running out of inodes mid-install.

`$STORE` is setgid, so EDFs and data created there are group-readable: one project member installs, the rest run.
