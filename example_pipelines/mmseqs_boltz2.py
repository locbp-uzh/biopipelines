# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 


"""
This pipeline shows how to run fold the sequences with Boltz2 using an MSA computed with MMseqs2.

The MMseqs2 server requires a lot of resources (800 GB ram) and takes 20-30 min to warm up depending on the filesystem. In BP, it is implemented in a client-server way, in which the client submits a job (potentially with many queries), the server picks it up, computes the MSA, and sends it to the cluster. There are two way to run it:
1. Implicit server: Only run the client, which will automatically submit a server if one is not found (see _mmseqs2_cpu_server.py). Easier, but we have to make sure the client is given enough time so that we know for sure the server will finish by the end of it.
2. Explicit server: Run the server as a service, within the pipeline. Certain, but different pipelines have to run the server warmup independently (moving 800 GB to RAM)
"""

from biopipelines.pipeline import *
from biopipelines import MMseqs2, MMseqs2Server, Boltz2

SEQ=["GNSKKHNLILIGAPGSGKGTQCEFIKKEYGLAHLSTGDMLREAIKNGTKIGLEAKSIIESGNFVGDEIVLGLVKEKFDLGVCVNGFVLDGFPRTIPQAEGLAKILSEIGDSLTSVIYFEIDDSEIIERISGRCTHPASGRIYHVKYNPPKQPGIDDVTGEPLVWRDDDNAEAVKVRLDVFHKQTAPLVKFYEDLGILKRVNAKLPPKEVTEQIKKIL"]

# Implicit server
with Pipeline(project="Examples",
              job="MMseqs-Boltz"):
    Resources(time="24:00:00", # hopefully enough time so that the server can run 
              memory="8GB")
    sequences = Sequence(SEQ)
    msas = MMseqs2(sequences=sequences) # this will submit a server
    # Now we need a GPU. This will be a separate job, dependent on the msa job
    Resources(gpu="any", 
            time="4:00:00",
            memory="16GB")
    boltz = Boltz2(proteins=sequences,
                        msas=msas) 
    
# Explicit server
with Pipeline(project="Examples",
              job="MMseqs-Boltz"):
    with Service(): # what follows a service (Sequence, MSA) runs as soon as the service (MMseqs2Server) is running, not after its completion
        Resources(gpu="none",
                  time="24:00:00",
                  memory="800GB",
                  cpus=32)
        MMseqs2Server("cpu", idle_timeout=1800)
    Resources(time="1:00:00", # this will likely finish in < 3 min.
              memory="8GB") # very low resources so we know for sure it will run before the idle timeout
    sequences = Sequence(SEQ)
    msas = MMseqs2(sequences=sequences) # this will query the existing server
    # Now we need a GPU. This will be a separate job, dependent on the msa job
    Resources(gpu="any", 
            time="4:00:00",
            memory="16GB")
    boltz = Boltz2(proteins=sequences,
                        msas=msas) 