# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu tool: emulate a protein's equilibrium structural ensemble from
sequence with a generative deep-learning model.

For each input sequence, BioEmu samples `num_samples` statistically
independent conformers approximating the equilibrium distribution. The
wrapper splits the sampled trajectory into per-conformer PDBs (a
`structures` stream, ids <seq>_<1..N>) and also keeps the compact
trajectory (`trajectories` stream: one samples.xtc + topology.pdb per
sequence). Side-chain reconstruction is optional.

Reference: Lewis et al. (2024) BioEmu. https://github.com/microsoft/bioemu
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
    from ._weights_cache import ensure_weights_block, link_weights_block
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern
    from _weights_cache import ensure_weights_block, link_weights_block


# BioEmu's inlined ColabFold downloads the AF2 params (~3.5 GB) to a hardcoded
# ~/.cache/colabfold/params on first sample() and skips the download only if a
# params/download_finished.txt marker exists there — there is no API/env-var to
# redirect it. So at install we populate the SHARED cache (download-if-absent,
# same params set every AF2 tool uses) and point ~/.cache/colabfold/params at
# it, touching the marker so BioEmu reuses the weights instead of re-fetching.
AF2_PARAMS_URL = "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"


class BioEmu(BaseConfig):
    """
    BioEmu: generative emulation of protein equilibrium ensembles from sequence.

    Inputs:
        sequences: amino-acid sequences (sequences stream or DataStream).
                   Each id's sequence is read from the stream's `sequence`
                   column. A precomputed MSA can be supplied per sequence by
                   pointing the sequence value at an .a3m file path.
        num_samples: conformers to sample per sequence (default 10).
        batch_size: sampling batch size on the GPU (default 10).
        reconstruct_sidechains: if True, run bioemu.sidechain_relax to produce
                   full-atom models (needs the bioemu[md] extra, installed by
                   default via BioEmu.install(); slower).
                   Default False (backbone-frame output).

    Install:
        BioEmu.install()             # installs bioemu[md] (sidechain relax) by default
        BioEmu.install(md=False)     # skip the heavier MD/sidechain extra
        filter_samples: drop unphysical conformers (long bonds, clashes);
                   default True.
        msa_host_url: override the ColabFold MMseqs2 MSA server (default None).

    Outputs:
        Streams:
            structures:   per-conformer PDB, ids <seq>_<1..num_samples>.
            trajectories: one <seq>.xtc trajectory + <seq>_topology.pdb per sequence.
        Tables:
            summary:  id | sequence.id | n_samples | n_residues
            missing:  id | cause
    """

    TOOL_NAME = "BioEmu"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, md=True, **kwargs):
        """Install the bioemu env. ``md=True`` (default) also installs the
        ``bioemu[md]`` extra needed for sidechain reconstruction
        (reconstruct_sidechains=True); pass ``BioEmu.install(md=False)`` to skip
        the heavier MD/sidechain dependencies."""
        biopipelines = folders.get("biopipelines", "")
        bioemu_check = f'{env_manager} run -n bioemu python -c "import bioemu.sample" >/dev/null 2>&1'

        # Reuse the shared AF2-params cache instead of letting BioEmu download
        # its own ~3.5 GB copy. Populate the shared cache if absent (so this
        # works even when BioEmu is the first AF2 tool installed), then link
        # BioEmu's hardcoded ~/.cache/colabfold/params at it + drop the marker.
        af2_params_root = folders.get("AlphaFoldParams", "")
        params_dir = f"{af2_params_root}/params"
        cf_cache = "$HOME/.cache/colabfold"
        params_marker = f"{cf_cache}/params/download_finished.txt"
        params_block = ensure_weights_block(
            dest_dir=params_dir,
            sentinel=f"{params_dir}/params_model_3.npz",
            download_cmds=f'curl -fsSL "{AF2_PARAMS_URL}" | tar x -C "{params_dir}"',
            label="AlphaFold2 params",
        )
        link_block = link_weights_block(
            target_dir=params_dir,
            link_path=f"{cf_cache}/params",
            marker=params_marker,
            label="ColabFold params",
        )
        # Two independent halves — env and params cache — so each is repaired
        # only when its own artefact is missing. On Colab the params live on a
        # persistent Drive while the env (and the ephemeral $HOME symlink that
        # points at them) are recreated every session; the old env-only skip
        # exited before re-linking, leaving BioEmu unable to find its weights.
        # The marker (download_finished.txt) sits on the ephemeral $HOME side,
        # so it proves the *link* is in place, not just that the params exist.
        params_check = f'[ -e "{params_marker}" ]'
        skip = "" if force_reinstall else f"""# Check if already installed (env + params cache both wired)
if {bioemu_check} && {params_check}; then
    echo "BioEmu already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("bioemu", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("bioemu", env_manager, biopipelines)
        # The bioemu[md] extra (sidechain reconstruction) is installed on top of
        # the base env when md=True so it can be toggled off at install time.
        md_block = (
            f'echo "Installing bioemu[md] extra (sidechain reconstruction)"\n'
            f'{env_manager} run -n bioemu pip install "bioemu[md]"\n'
            if md else ""
        )
        # Create the env only if absent (force_reinstall removed it above). The
        # params wiring then runs unconditionally below so a persisted env with
        # a lost $HOME symlink still gets re-linked.
        env_install_block = f"""if ! {bioemu_check}; then
{env_block}

{md_block}else
    echo "bioemu environment already present, skipping creation."
fi"""
        return f"""echo "=== Installing BioEmu ==="
{skip}{remove_block}
{env_install_block}

# Wire BioEmu's ColabFold cache to the shared AF2-params cache. Runs every
# install (not gated on the env) so an ephemeral-$HOME symlink is rebuilt even
# when the env itself was untouched.
echo "Wiring BioEmu's ColabFold cache to the shared AF2-params cache"
{params_block}
{link_block}

# Verify installation
if {bioemu_check} && {params_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== BioEmu installation complete ==="
else
    echo "ERROR: BioEmu verification failed (cannot import bioemu.sample or params cache not wired)"
    exit 1
fi
"""

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    trajectories_map = Path(lambda self: self.stream_map_path("trajectories"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_bioemu.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 num_samples: int = 10,
                 batch_size: int = 10,
                 reconstruct_sidechains: bool = False,
                 filter_samples: bool = True,
                 msa_host_url: str = "",
                 **kwargs):
        self.sequences = sequences
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        self.num_samples = int(num_samples)
        self.batch_size = int(batch_size)
        self.reconstruct_sidechains = bool(reconstruct_sidechains)
        self.filter_samples = bool(filter_samples)
        self.msa_host_url = msa_host_url
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if self.num_samples < 1:
            raise ValueError("num_samples must be >= 1")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"SEQUENCES: {len(self.sequences_stream)}")
        lines.append(f"NUM SAMPLES PER SEQUENCE: {self.num_samples}")
        lines.append(f"RECONSTRUCT SIDECHAINS: {self.reconstruct_sidechains}")
        lines.append(f"FILTER SAMPLES: {self.filter_samples}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.sequences_stream.save_json(self.sequences_json)
        sidechain_arg = " --reconstruct-sidechains" if self.reconstruct_sidechains else ""
        nofilter_arg = "" if self.filter_samples else " --no-filter-samples"
        msa_arg = f' --msa-host-url "{self.msa_host_url}"' if self.msa_host_url else ""
        upstream_missing_path = self._get_upstream_missing_table_path(self.sequences)
        upstream_missing_flag = f' --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# BioEmu equilibrium-ensemble sampling script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Sampling ensembles for {len(self.sequences_stream)} sequence(s), {self.num_samples} each"
python "{self.helper_py}" \\
    --sequences-json "{self.sequences_json}" \\
    --num-samples {self.num_samples} \\
    --batch-size {self.batch_size} \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map "{self.structures_map}" \\
    --trajectories-dir "{self.stream_folder('trajectories')}" \\
    --trajectories-map "{self.trajectories_map}" \\
    --scratch-dir "{self.extras_path()}" \\
    --summary-csv "{self.summary_csv}" \\
    --missing-csv "{self.missing_csv}"{sidechain_arg}{nofilter_arg}{msa_arg}{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        suffix_pattern = f"<1..{self.num_samples}>"
        structure_ids = generate_multiplied_ids_pattern(
            self.sequences_stream.ids, suffix_pattern,
            input_stream_name="sequences",
        )
        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )
        trajectories = DataStream(
            name="trajectories",
            ids=self.sequences_stream.ids,
            files=[self.stream_path("trajectories", "<id>.xtc")],
            map_table=self.trajectories_map,
            format="xtc",
        )
        tables = {
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "sequence.id", "n_samples", "n_residues"],
                description="BioEmu per-sequence ensemble summary",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "structures": structures,
            "trajectories": trajectories,
            "tables": tables,
            "output_folder": self.output_folder,
        }
