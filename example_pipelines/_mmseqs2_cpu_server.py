# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
This pipeline should not be modified as it is used by the cpu server
"""

from biopipelines.pipeline import *
from biopipelines import MMseqs2Server

with Pipeline(project="MMseqs2Server",
              job="CPU",
              description="Runs MMseqs2 server on CPU for 1 against many sequence alignment"):

    # CPU mode runs the same colabfold_search pipeline as GPU mode but with no
    # GPU: the prefilter/search runs on CPU and the whole thing is bound by RAM,
    # not cores. Crucially the CPU indices are larger than the GPU ones — the GPU
    # build drops the k-mer prefilter table (--index-subset 2) because the GPU
    # does the prefilter, but CPU mode keeps it, so the .idx total is ~746GB
    # (uniref30 ~220GB + envdb ~526GB) vs the GPU's ~410GB. To match the public
    # ColabFold server's speed the whole index must stay resident in the page
    # cache (warmed once at startup via cat .idx; the rest is mmapped via
    # --db-load-mode 2). The cluster has many ~1.5TB nodes; request 800GB to hold
    # the ~746GB index plus headroom (900GB+64cpu queued ~11h on the busy big-RAM
    # partition, so keep the ask modest). The search is CPU-bound across
    # align/expandaln/result2msa, so take a moderate core count.
    Resources(gpu="none",
              time="24:00:00",
              memory="800GB",
              cpus=32)

    # idle_timeout: auto-shut-down after 30 min with an empty queue so the
    # large-RAM node is released and SLURM priority isn't penalised.
    MMseqs2Server("cpu", idle_timeout=1800)
