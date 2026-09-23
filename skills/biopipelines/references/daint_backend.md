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

---

## Working on Daint day to day

Daint is a SLURM cluster, so `cluster_backend.md` applies for the ssh/scp workflow. This file covers what differs: aarch64, whole-node billing, and a container runtime that is neither apptainer nor plain podman.

Select the variant with `BIOPIPELINES_CONFIG_VARIANT=daint` or `Pipeline(config="daint")`. `config.daint.yaml` carries the paths and the environment mapping. Its `environments:` block is a **routing table, not a verification record**: each of its ~51 keys says which venv (or, when empty, which EDF image's interpreter) a tool activates here. A tool absent from it falls back to the `biopipelines` venv, which is correct for any BP-native tool but unverified for one needing its own env. Presence in the block does *not* mean the tool has been run — the verified list is the much shorter one under *What runs here* below (17 tools at last probe), and *Blocked, with the reason* lists what cannot work here at all.

## One-time setup

1. **Sign a certificate.** CSCS issues short-lived certificates, not permanent keys. They expire after **1 day**, so this is a daily step — when anything fails with `Permission denied (publickey)`, re-sign before debugging further.

   ```bash
   cscs-key sign                                            # opens a browser for OIDC + MFA
   ssh-keygen -L -f ~/.ssh/cscs-key-cert.pub | grep Valid   # check the window
   ```

   `cscs-key` from https://github.com/eth-cscs/cscs-key/releases; keypair at `~/.ssh/cscs-key`.

2. **Add the ssh aliases** (per https://docs.cscs.ch/access/ssh/). Daint is reached through the `ela` jump host:

   ```
   Host ela
       HostName ela.cscs.ch
       User <your-cscs-username>
       IdentityFile ~/.ssh/cscs-key
       IdentitiesOnly yes

   Host daint
       HostName daint.alps.cscs.ch
       User <your-cscs-username>
       ProxyJump ela
       IdentityFile ~/.ssh/cscs-key
       IdentitiesOnly yes
   ```

   Always use the alias. `ssh user@daint.alps.cscs.ch` skips the config block, offers no key, and fails with `Permission denied (publickey)`.

3. **Set your username in the overlay**, not the committed config — `config.cluster.yaml` and `config.daint.yaml` would otherwise both claim the same Unix name and auto-detection refuses to guess:

   ```bash
   bp-config set machine.username <you> --variant daint
   bp-config set machine.billing_account <project> --variant daint
   ```

   The account is mandatory: CSCS rejects jobs without `--account`.

4. **Clone and install.** Use a venv on `$SCRATCH`, not conda and not `$STORE` — see Storage below.

   ```bash
   /usr/bin/python3.11 -m venv $SCRATCH/venvs/biopipelines   # no `python` on login nodes
   source $SCRATCH/venvs/biopipelines/bin/activate
   pip install -r environments/biopipelines.pip.daint.txt    # aarch64 wheels, not the conda yaml
   pip install -e .
   ```

5. **Install tools** with `example_pipelines/install_daint.py`, or `bp-warm <Tool>` individually.

## Partitions

Per https://docs.cscs.ch/clusters/daint/:

| Partition | Nodes | Max nodes/job | Walltime | Use for |
|---|---|---|---|---|
| `debug` | 24 | 2 | 30 min | every bring-up test — it never queued in any probe |
| `normal` | unlimited | — | 24 h | production; deep queue (774→893 pending within one evening) |
| `low` | unlimited | — | 24 h | needs QoS `slowdown`/`stop`, which a plain association does not grant |
| `xfer` | 2 | 1 | 24 h | store↔scratch transfers; shared nodes, no GPUs, 250 GB |

Nodes are **not** shared in `normal`/`debug`/`low`; `xfer` nodes are. `normal` is the only general compute partition and every node in it has 4 GH200s, so a pure-CPU job still occupies and is billed for four idle GPUs.

Expect to queue. Priority is dominated by Fairshare, and a usage-based project shows `RawShares=1` in `sshare` — that looks like an exhausted allocation and is not. Check the real budget at https://account.cscs.ch, and note `sreport` defaults to **CPU-hours**: pass `--tres=node` for the billed unit, or a 288-core node overstates usage ~100×.

## Storage

| Path | Quota | Retention | Use for |
|---|---|---|---|
| `$HOME` | 50 GB, 500k inodes | backed up | repo, scripts, configs |
| `$SCRATCH` | 150 TB, 1M inodes | **purged at 30 days** | venvs, databases, job output |
| `$STORE` | 1 TB, **150k inodes** | backed up daily, setgid | EDFs, images, weights, shared data |

**The `$STORE` inode limit is what binds, not the terabyte.** Python environments are inode-heavy and nearly weightless: seven venvs cost 101,223 files for 15 GB, and adding PyMOL failed with `Disk quota exceeded` while the filesystem was 98.6% empty by size. So `venv_root` points at scratch, and the cost is the purge — a tool unused for a month needs reinstalling, which is cheaper than running out of inodes mid-install.

Scratch purges on **access** time, so a regularly queried file survives — but only the files actually opened. Anything untouched ages out independently, which can leave a half-purged set that looks present.

## Nodes are billed whole — pack them

`OverSubscribe=EXCLUSIVE` everywhere: a job asking for one GPU is allocated and billed for the entire node (`AllocTRES = billing=288,cpu=288,gres/gpu=4`). Every node is 4× GH200, 288 CPUs, 870 GB. There are **no A100s**, so `Resources(gpu="A100")` has no meaning here.

Splitting a sweep into N single-GPU jobs therefore burns N whole nodes. Fill one instead with `Parallel(pack=N)` + `Run()`, which emits concurrent `srun` job steps:

```python
with Parallel(pack=4):
    Resources(gpu="gh", time="06:00:00")      # the allocation, node-shaped
    for chunk in designs.chunks(4):
        with Run(cpus=72, gpus=1, memory="50GB"):   # one task inside it
            Boltz2(proteins=chunk, ...)
```

`Boltz2` has no `structures=` parameter — its typed inputs are `proteins`, `ssDNA`/`dsDNA`, `ssRNA`/`dsRNA`, `ligands` and `msas`, and an unknown keyword is swallowed by `**kwargs` rather than rejected. See `example_pipelines/daint_packed_design.py` for the full working shape.

Seven things about packing are not guessable, all verified on hardware:

- **The configuration batch needs its own `Resources`.** A pipeline whose first tools are cheap (`PDB`, `Ligand`, `Load`) still submits them as a separate SLURM job, and with no `Resources` in scope before them it takes the framework default of `time="24:00:00"`, `memory="15GB"`. On an exclusive site that queues a whole node against a full day for a few seconds of work. Set `Resources(time="00:30:00", memory="20GB")` at the top of the pipeline, before the `Parallel` block.
- **`--cpus-per-task` is mandatory.** A step without it gets ONE core out of 288 — correct results, crippled throughput, nothing in the log says so. `Run(cpus=)` sets it; the default divides `machine.cores_per_node` by `pack`.
- **CPU-only steps must explicitly decline GPUs, and every step needs a memory slice.** On Daint, `--gpus-per-task=0` inherits all four allocation GPUs; packed CPU steps use `--gres=none`. A step without explicit memory similarly inherits the allocation's full RAM and serializes otherwise independent siblings. The default memory share is proportional to CPUs, while `Run(memory=)` supports asymmetric layouts such as a 300 GB MMseqs2 server beside four 50 GB GPU workers. **Concurrent slices must sum inside the allocation** — 780 + 4x200 against an 849 GB ask stalls with `step creation still disabled, retrying (Requested nodes are busy)` rather than failing.
- **`--cpus-per-task` bounds the step, not the tools inside it.** `boltz predict` sizes its input-processing pool from the visible machine, `min(n_inputs, node_cores)`, ignoring `--cpus-per-task`: 158 inputs logs "158 threads" and fits in 16 GB, while 858 inputs logs "288 threads" and the OOM killer takes it at the same 16 GB — `Killed`, then a `BrokenPipeError` cascade from the collapsing pool. So the memory a packed Boltz2 task needs scales with its **input count**, not with a fixed figure. Budget ~0.5 GB per concurrent thread and keep the sum over concurrent steps under the node's 870 GB: five steps at 200 GB plus a 780 GB server asks for 1580 GB on an 849 GB allocation, and SLURM silently stalls with `step creation still disabled, retrying (Requested nodes are busy)` rather than failing. Chunk the inputs instead when the sum does not fit. Any tool reading `os.cpu_count()` rather than `SLURM_CPUS_PER_TASK` has the same shape.
- **Nested steps fail.** Wrapping a `Run()` in an outer step with per-tool `--environment=` steps inside gives `Unable to satisfy cpu bind request` (rc=192) for every task whose mask is offset into the node — 3 of 4 died, the survivor holding core 0, so a single-task probe reports success. Each tool is its own sibling step instead, which also preserves per-tool container selection.
- **`pack=N` is capacity, not placement** when tasks exceed one node: 6 tasks at pack=4 over 2 nodes landed 3/3, not 4/2. Within one node's worth, `--ntasks-per-node` and `--distribution=block` do hold them together.
- **GPUs are requested once per node, not per task** on an exclusive site. `--gpus-per-task` multiplies by task count, so a CPU-only task alongside four GPU ones asks for 5 GPUs and SLURM rejects the allocation outright.

`DataStream.chunks(count)` / `chunks(size=)` splits a stream for fan-out. Tasks cannot consume each other's output — each is a separate subshell — so gather results in a later pipeline with `Load` + `Pool`.

## Environments: venv, conda shim, or image

There are three routes, and the choice is not stylistic.

**A plain venv** covers anything pip-installable. `pip install torch` resolves a CUDA-enabled aarch64 wheel (`torch 2.13.0+cu130`, `cuda True` on a GH200), so most GPU tools need no container at all. `machine.env_manager.name: "venv"` activates by path — venv has no registry, so `venv_root` is what turns an `environments:` entry into one.

**A conda env wrapped in a venv shim** covers packages with no aarch64 wheel — PyMOL, `mkdssp`, `openbabel` all exist on conda-forge, and micromamba has an aarch64 build. Set `machine.env_manager.conda_binary` and an env yaml naming conda packages is built with conda, then wrapped in a venv created by *its* python: that inherits the packages via `--system-site-packages` and supplies the `bin/activate` the framework sources. Binaries are not inherited and need symlinking (`dssp` does this for `mkdssp`).

**A prebuilt `.sqsh` image** covers what neither can — BioEmu is the example. Five images cost 5 inodes and 28 GB against ~101k inodes for equivalent venvs, they survive the purge, and they are genuinely shareable where a venv on personal scratch is not. Recipes live in `environments/scripts/`; the images themselves belong on `$STORE`, never in git.

Two traps common to all three:

- **The login nodes have no `python3.11-dev`**, so anything without an aarch64 wheel fails to compile with `fatal: Python.h: No such file`. Copy the headers out of a `python:3.11` image once and point `folders.infrastructure.python_headers` at them; `_env_install_block` exports `CPATH` for every venv install.
- **The host python's stdlib is incomplete.** SUSE splits `dbm`, `curses` and `ctypes`-adjacent modules into separate distro packages, and what is left behind is an `_import_failed` stub that raises `ImportError` telling you to run `zypper` — which no user here can. `import dbm` fails on `/usr/bin/python3.11`, and every `--system-site-packages` venv inherits that. This is not a wheel problem and pip cannot fix it: check the *stdlib*, not just PyPI, when a pure-python dependency fails to import. `dbm` in particular is recoverable, since `dbm/__init__.py` guards each backend separately and `dbm.dumb` is pure python — ADMETAI vendors those two files when the probe fails.
- **venvs are not relocatable.** A venv bakes its absolute path into `bin/activate`, `pyvenv.cfg`, and every console-script shebang. Moving one leaves `<newpath>/bin/python` working while `source <newpath>/bin/activate` silently puts the *old* path on `PATH`. Rebuild rather than move — it is also far faster than a cross-filesystem `mv` of 100k small files.
- **Build a venv where its interpreter will run.** One built inside a container bakes the container's python into `bin/python3` and `pyvenv.cfg`, and on the host that symlink dangles — `PATH` skips the broken entry and `python` falls through to whatever env was already active, so an installed package surfaces as `ModuleNotFoundError` instead of anything naming the environment. The tell is the step log's own header, where `Environment:` (`VIRTUAL_ENV`, intent) and `Python:` (`which python`, fact) name different venvs; check that before the traceback. A venv whose python is absent from the host is correct only for a tool with an `edf:` entry, which runs inside the image — so audit `bin/python` against the `edf:` block, not on its own.

## GPU tools: the Container Engine

Not podman. It is a Slurm plugin (pyxis + enroot) invoked as `srun --environment=<file.toml>`, reading TOML Environment Definition Files. The `edf:` block maps a tool to one:

```yaml
edf:
  Boltz2: "<group>/edf/ngc-pytorch.toml"
```

```toml
image = "nvcr.io#nvidia/pytorch:25.06-py3"   # registry ref, or a local .sqsh path
mounts = ["/capstor", "/iopsstor"]
workdir = "/capstor/scratch/cscs/<user>"
[env]
PATH = "/opt/<env>/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
NVIDIA_VISIBLE_DEVICES = "all"
NVIDIA_DRIVER_CAPABILITIES = "compute,utility"
```

NGC publishes aarch64/GH200 PyTorch, so torch+CUDA arrive prebuilt rather than being ported.

- The wrap goes around the **whole tool script**, not the inner python call: the Container Engine defines a job *step*, so an inner prefix would leave the env activation outside the container. This is why `edf:` is separate from `containers:` (apptainer's command prefix), which stays empty here.
- A venv layered on an image **must be built inside that image** — its `bin/python` symlinks into the image and will not run on the host. An `environments:` entry left empty means the image supplies the interpreter and nothing is activated.
- **Reference EDFs by absolute path.** A bare name resolves against each user's private `~/.edf`, so it only works for whoever wrote it.
- `PATH` and the `NVIDIA_*` pair are load-bearing for a self-built image: the Engine *replaces* the environment rather than extending it, and enroot injects the host driver only when those are set — without them there is no `libcuda.so` and a GPU job silently runs on the CPU. NGC images set them internally, which is why `ngc-pytorch.toml` gets away without.
- `workdir` **overrides the cwd of whatever called `srun`**, so a `cd` before the call does nothing.

## MMseqs2 databases

Built and working on scratch. The numbers, measured rather than estimated:

| | |
|---|---|
| Download (`step="databases"`, `xfer`) | 234 GB, **2h33m** |
| Index build (`step="build"`, whole node) | **6h19m** |
| Footprint once built | 1.9 TB, of which 719 GB is the three `.idx` files (uniref30 220 GB, envdb 498 GB, pdb100 1.4 GB) |

`$STORE`'s 1 TB cannot hold them, so they stay on scratch and a purge costs ~9 h to rebuild. The index files are `stripe_count: 1`; `lfs setstripe -c 4` measured 737 MB/s against 491, but the rewrite costs 719 GB and 16 stripes measured *worse* (200 MB/s), so this was left alone.

**Size the server for the database it actually maps, not for the 719 GB total.** `MMseqs2Server(database=)` defaults to `uniref30_2302_db` — a 220 GB index. Measured on a live packed node: server plus four GPU tasks used **351 GB of 856**, with 271 GB of that the mmap'd index. `Run(memory="300GB")` is enough for the default server; reserving 780 GB starves the sibling steps for nothing.

**envdb does not run here.** Its 498 GB index cannot be served on this node, and several strategies have been tried without success — treat it as a closed path rather than a sizing question. UniRef30 is what the local server maps.

Server startup is not instant, which is why the server belongs in a packed block rather than a job of its own — see `example_pipelines/daint_packed_design.py`, where it warms while RFdiffusion3 and ProteinMPNN work. Clients block until it advertises readiness, so no coordination is needed.

`MMseqs2(server_url=...)` is the alternative: any ColabFold-protocol endpoint over HTTP, no local databases and no startup cost. Credentials come from `MMSEQS2_SERVER_USER`/`MMSEQS2_SERVER_PASSWORD` and are refused over plain http.

## What runs here

These 17 are verified by running, not just importing (the `environments:` block routes many more tools than this — see the note at the top): PDB, Sequence, Ligand, RDKit, Distance, Panda, SASA, PyMOL (SASA 6536.9 Å² on lysozyme, 1LYZ→1AKI at 0.532 Å), DSSP (`mkdssp 4.6.1`), Boltz2, BoltzGen, RFdiffusion3, ESMFold2, ProteinMPNN, LigandMPNN, PoseBusters, Vina.

Docking is Vina, not GNINA. `Vina` (a `Gnina` subclass with `backend="vina"`) runs on CPU from a conda-forge env; verified on 9RTM + tetramethylrhodamine, crystal-ligand autobox best −10.4 kcal/mol with the pose centroid 0.6 Å from the crystal ligand. Two bring-up bugs, both outside the tool: a venv shim inherits conda *packages* but not *executables*, so `vina`/`obabel` needed symlinking (`_link_conda_binaries`); and a 5-character CCD code (`A1EI4`) overflows the 3-column residue field, shifting every later column so `parse_pdb_file` silently dropped every ligand atom — that fix is in `pdb_parser.py` and applies repo-wide, not just here.

Blocked, with the reason:

- **GNINA** — every release asset is a prebuilt x86-64 binary (currently `gnina.cuda12.8.static`); there is no aarch64 build, and a source port would mean rebuilding its bundled Caffe fork and CUDA kernels for sm_90. `Vina` covers docking instead, without the CNN rescoring.
- **NeuralPLexer / DiffDock / DynamicBind** — all pin PyG wheels from `data.pyg.org`, which serves `linux_x86_64` only; NeuralPLexer additionally pins `gcc_linux-64`/`gxx_linux-64` and an x86-64 pytorch3d wheel.
- **RFdiffusion** — the RosettaCommons image is x86-64 only (`no image found in image index for architecture "arm64"`), and the pip route dies on DGL, which ships an aarch64 wheel that loads a prebuilt C++ library named for an exact torch version (`libgraphbolt_pytorch_2.12.1.so`) with no matching build. RFdiffusion3 covers the same ground.
- **BioEmu** — imports and sees the GPU, but at sample time `get_colabfold_embeds` builds its own ColabFold venv pinning `colabfold==1.5.4` and x86-64 `nvidia-*-cu12` wheels, and the resolve dies on `tensorflow-cpu`. The lab's `md-bioemu` image runs on a GH200, so `edf:` points there instead.
- **PLM_Sol** — pins `python=3.8` with hard-pinned deps that cannot coexist with a 3.11+ env.

An import check is not a run. BioEmu passed `import bioemu.sample` and still failed on a real sampling job — twice, on different causes.
