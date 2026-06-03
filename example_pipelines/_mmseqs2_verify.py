# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""End-to-end verification of the MMseqs2 client + colabfold_search server.

Two MMseqs2 client jobs run in parallel, each on a couple of query sequences.
The clients are CPU jobs that submit FASTAs to the MMseqs2 server queue and wait
for the resulting MSAs; the server itself (GPU) is auto-submitted on demand.

Running two clients in parallel also checks the single-server lock: only ONE
MMseqs2Server job should be submitted, not one per client.
"""

from biopipelines.pipeline import *
from biopipelines import MMseqs2

with Pipeline(project="Examples",
              job="MMseqs2Verify",
              description="Verify MMseqs2 client + colabfold_search server end to end"):

    with Parallel():
        # Client A — CPU, waits on the server.
        Resources(gpu=None, time="6:00:00", memory="16GB")
        MMseqs2(sequences=Sequence(
            ["GNSKKHNLILIGAPGSGKGTQCEFIKKEYGLAHLSTGDMLREAIKNGTKIGLEAKSIIESGNFVGDEIVLGLVKEK",
             "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAK"],
            ids=["qA1", "qA2"]))

        # Client B — separate CPU job, submitted at the same time as A.
        Resources(gpu=None, time="6:00:00", memory="16GB")
        MMseqs2(sequences=Sequence(
            ["FDLGVCVNGFVLDGFPRTIPQAEGLAKILSEIGDSLTSVIYFEIDDSEIIERISGRCTHPASGRIYHVKYNPPKQ",
             "GSHMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHS"],
            ids=["qB1", "qB2"]))
