# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Install the ColabFold MSA databases for the MMseqs2 server / colabfold_search.

MMseqs2Server.install() takes a `step` argument selecting the stage:

- step="databases": download stage only — fetches the ColabFold tarballs and
  rsyncs the mmCIF snapshot into the ColabFoldDatabases folder root. CPU-only, no
  GPU. Mode-independent: shared by both the cpu and gpu builds.
- step="build", mode="gpu": build the indexed, GPU-padded databases (uniref30 +
  colabfold_envdb) from the downloaded tarballs into ColabFoldDatabases/gpu/.
  Needs a GPU node, so it runs in its own GPU batch. It also installs a
  GPU-capable MMseqs2 (release >=16) into the MMseqs2 folder if absent.
- step="build", mode="cpu": build the non-padded databases into
  ColabFoldDatabases/cpu/. Runs on a CPU node (no GPU); the MMseqs2 server's CPU
  mode serves from this folder, keeping the indices resident in RAM.

The downloads stay at the ColabFoldDatabases root and are symlinked into each
mode subfolder, so the large tarballs aren't duplicated. All stages are
idempotent (resumed via marker files), so re-running skips work already done.
The stages are split into separate Resources() batches because only the GPU
build needs a GPU.
"""

from biopipelines.pipeline import *
from biopipelines import MMseqs2Server

with Pipeline(project="Setup", job="InstallMMseqs2Databases",
              description="Download + build the ColabFold MSA databases"):

    # Stage 1 — download the databases (CPU only). I/O-bound (aria2c parallelizes
    # over connections, not cores), so a couple of CPUs is plenty; the long wall
    # time is for the slow mmCIF rsync, not compute.
    Resources(time="48:00:00", memory="16GB", cpus=2)
    MMseqs2Server.install(step="databases")

    # Stage 2 — build the CPU (non-padded) indexes into cpu/ (waits for stage 1).
    # No GPU needed.
    Resources(time="24:00:00", memory="64GB", cpus=16)
    MMseqs2Server.install(step="build", mode="cpu")

    # Stage 3 — build the GPU-padded indexes into gpu/ (waits for stage 2). Any
    # GPU works; the build needs a GPU only to produce the GPU-padded format, not
    # for heavy compute, so gpu="any" queues fastest.
    Resources(gpu="any", time="24:00:00", memory="64GB", cpus=16)
    MMseqs2Server.install(step="build", mode="gpu")
