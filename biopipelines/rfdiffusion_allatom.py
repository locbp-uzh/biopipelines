# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion-AllAtom configuration for ligand-aware protein design.

Handles RFdiffusion-AllAtom workflows with ligand contexts, PPI design,
and all-atom generation capabilities.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
    from .combinatorics import generate_multiplied_ids_pattern
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
    from .guiding_potentials import RFdiffusionAllAtomGuidingPotential, render_guiding_potentials
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern
    from input_standardization import resolve_basic_input
    from biopipelines_io import Resolve, TableReference
    from ligand import Ligand
    from guiding_potentials import RFdiffusionAllAtomGuidingPotential, render_guiding_potentials


def _normalize_selection_arg(name, value):
    """Normalize a per-PDB selection argument to a (kind, token) pair.

    ``contigs`` / ``inpaint`` / ``inpaint_str`` each accept either a plain
    string (broadcast to every input PDB) or a per-PDB column reference (a
    ``TableReference`` from ``tool.tables.X.col``, or a ``(TableInfo, "col")``
    tuple). A reference is resolved per-id at runtime; a literal is emitted
    directly. Returns ``("literal", str)`` or
    ``("table", "TABLE_REFERENCE:path:col")``.
    """
    if isinstance(value, TableReference):
        return ("table", str(value))
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return ("table", str(TableReference(value[0].info.path, value[1])))
    if isinstance(value, str):
        return ("literal", value)
    raise ValueError(
        f"{name} must be a string or a (table, column) reference, got {type(value)}"
    )


def _selection_cli_arg(kind_token) -> str:
    """CLI arg for the contigs resolver: the token, or '-' when empty."""
    kind, token = kind_token
    if kind == "literal" and not token:
        return "-"
    return token


class RFdiffusionAllAtom(BaseConfig):
    """
    RFdiffusion-AllAtom variant for ligand-aware protein design.

    Extends base RFdiffusion with support for ligand contexts and
    all-atom generation capabilities including PPI design.
    """

    TOOL_NAME = "RFdiffusionAllAtom"
    TOOL_VERSION = "1.1"

    # Typed builder for the guiding potential, e.g.
    #   RFdiffusionAllAtom.GuidingPotential.ligand_ncontacts(weight=3, r_0=8, d_0=4)
    # Restricted to ligand_ncontacts (the only type AllAtom implements).
    GuidingPotential = RFdiffusionAllAtomGuidingPotential

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("RFdiffusionAllAtom", "")
        parent_dir = os.path.dirname(repo_dir)
        biopipelines = folders.get("biopipelines", "")
        # Extra deps that layer onto SE3nv, kept declarative (not hardcoded here).
        pip_reqs = f"{biopipelines}/environments/rfdiffusion_allatom.pip.txt"
        # Importable-rfdiffusion probe: the env merely existing is not enough —
        # the AllAtom extras layer on top of the rfdiffusion package, which
        # RFdiffusion.install() puts into SE3nv. Test the real precondition.
        se3nv_ready = f'{env_manager} run -n SE3nv python -c "import rfdiffusion" >/dev/null 2>&1'
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && [ -f "{repo_dir}/RFDiffusionAA_paper_weights.pt" ] && {se3nv_ready}; then
    echo "RFdiffusion-AllAtom already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing RFdiffusion-AllAtom ==="
{skip}# Repo + weights are independent of SE3nv — fetch them first, always, so they
# are in place even if the env still needs RFdiffusion.install().
mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/baker-laboratory/rf_diffusion_all_atom.git
fi
cd {repo_dir}
wget -nc http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFDiffusionAA_paper_weights.pt
git submodule init
git submodule update

# The AllAtom extras layer onto SE3nv's rfdiffusion install. Require it to be
# importable; if not, give a clear error (weights are already downloaded above,
# so finishing is just `RFdiffusion.install()` then re-running this).
if ! {se3nv_ready}; then
    echo "ERROR: SE3nv does not have an importable 'rfdiffusion' package."
    echo "       Run RFdiffusion.install() first, then re-run RFdiffusionAllAtom.install()."
    echo "       (RFdiffusion-AllAtom repo + weights have been downloaded already.)"
    exit 1
fi

# Layer the AllAtom-specific deps into the shared SE3nv env (declarative).
if [ -f "{pip_reqs}" ]; then
    {env_manager} run -n SE3nv pip install -r "{pip_reqs}"
else
    echo "WARNING: {pip_reqs} not found; skipping AllAtom extra deps."
fi

# Verify installation
if [ -f "{repo_dir}/RFDiffusionAA_paper_weights.pt" ] && {se3nv_ready}; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion-AllAtom installation complete ==="
    echo "Container mode: configure containers.RFdiffusionAllAtom in config.yaml"
else
    echo "ERROR: RFdiffusion-AllAtom verification failed (weights missing or SE3nv not ready)"
    exit 1
fi
"""

    # Lazy path descriptors
    #   main_table — standalone TableInfo CSV (tables/structures.csv).
    #   pdb_ds_json — config-time input DataStream serialization.
    main_table = Path(lambda self: self.table_path("structures"))
    inference_py_file = Path(lambda self: "run_inference.py")
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_table.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))
    substrate_json = Path(lambda self: self.configuration_path("input_substrate.json"))
    # Per-PDB contig resolution (see RFdiffusion for the contract).
    contigs_json = Path(lambda self: self.configuration_path("contig_options.json"))
    contig_args_json = Path(lambda self: self.configuration_path("contig_args.json"))
    contigs_py = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_contigs.py"))
    contigs_reader_py = Path(lambda self: self.pipe_script_path("resolve_rfdiffusion_contigs.py"))

    def __init__(self,
                 ligand: Union[str, DataStream, StandardizedOutput],
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 contigs: Union[str, TableReference, "tuple"] = "",
                 inpaint: Union[str, TableReference, "tuple"] = "",
                 num_designs: int = 1,
                 active_site: bool = False,
                 steps: int = 50,
                 partial_steps: int = 0,
                 reproducible: bool = False,
                 design_startnum: int = 1,
                 ppi_design: bool = False,
                 ppi_hotspot_residues: List[str] = None,
                 ppi_binder_length: int = None,
                 autogenerate_contigs: bool = False,
                 model_only_neighbors: bool = False,
                 num_recycles: int = 1,
                 scaffold_guided: bool = False,
                 align_motif: bool = True,
                 deterministic: bool = False,
                 inpaint_str: str = None,
                 inpaint_seq: str = None,
                 inpaint_length: int = None,
                 guiding_potentials: Union[str, "GuidingPotential", List, None] = None,
                 guide_scale: float = None,
                 guide_decay: str = None,
                 substrate: Union[str, DataStream, StandardizedOutput, None] = None,
                 **kwargs):
        """
        Initialize RFdiffusion-AllAtom configuration.

        Args:
            ligand: Ligand as a compounds stream (Ligand(code="ZIT") or any
                    compounds-producing tool). The 3-letter residue code is
                    read from the stream's `code` column at runtime and passed
                    to RFdiffusion-AllAtom's inference.ligand= hydra flag.
            pdb: Input PDB structure as DataStream or StandardizedOutput
            contigs: Contig specification (e.g., "A1-100,10-20")
            inpaint: Inpainting specification (same format as contigs)
            num_designs: Number of designs to generate
            active_site: Use active site model for small motifs
            steps: Inference diffusion steps (default 50, the Baker-lab default —
                trained at 200 but ~50 gives equivalent quality at ~4x speed; ~20
                is also reported as in-silico-equivalent)
            partial_steps: Partial diffusion steps
            reproducible: Use deterministic sampling
            design_startnum: Starting number for design numbering (default: 1)
            ppi_design: Enable protein-protein interaction design
            ppi_hotspot_residues: List of hotspot residues for PPI (e.g., ["A116","A150"])
            ppi_binder_length: Length of PPI binder
            autogenerate_contigs: Auto-infer fixed contig segments from input PDB
            model_only_neighbors: Only remodel residues neighboring the scaffold
            num_recycles: Number of diffusion recycles (iterative refinement)
            scaffold_guided: Enforce strict adherence to input scaffold geometry
            align_motif: Pre-align any functional motif before diffusion
            deterministic: Use fixed RNG seeds for reproducible outputs
            inpaint_str: Secondary-structure pattern for inpainting (e.g., "HHHEE")
            inpaint_seq: Sequence pattern for inpainting (e.g., "ACDEFG")
            inpaint_length: Target length for each inpainted region
            guiding_potentials: Guiding potential to pull the diffusing backbone
                onto the ligand so the new region wraps it instead of drifting away
                (``potentials.guiding_potentials``). Use the typed builder
                ``RFdiffusionAllAtom.GuidingPotential.ligand_ncontacts(...)``; a raw
                hydra string is also accepted::

                    guiding_potentials=RFdiffusionAllAtom.GuidingPotential.ligand_ncontacts(
                        weight=3, r_0=8, d_0=4)

                Tune ``weight`` 1-10 (higher = pull the design harder onto the ligand);
                r_0 ~8 Å (contact switching distance), d_0 ~4 Å (always-in-contact).
            guide_scale: Global multiplier on the guiding potential
                (``potentials.guide_scale``).
            guide_decay: How potential influence decays over the trajectory
                (``potentials.guide_decay``): "constant", "linear", "quadratic", or
                "cubic". "quadratic" is a good default.
            substrate: 3-letter residue code for the substrate_contacts potential
                (``potentials.substrate``). A bare string is the code; a Ligand /
                compounds stream resolves to its code at runtime. Defaults to the
                same code as ``ligand`` when a substrate_contacts potential is used,
                so you usually don't set it. Only emitted when guiding_potentials
                contains a substrate_contacts entry.
            **kwargs: Additional parameters

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | pdb | fixed | designed | source_fixed | plddt_mean | status
        """
        # Resolve optional pdb input — store stream for runtime resolution. All
        # input PDBs are iterated at execution time (see generate_script); never
        # index ids[0] at config time (breaks under lazy IDs).
        self.pdb_stream: Optional[DataStream] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_stream = pdb.streams.structures
            elif isinstance(pdb, DataStream):
                self.pdb_stream = pdb
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Core parameters
        self.contigs = contigs
        self.inpaint = inpaint
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.design_startnum = design_startnum

        # AllAtom-specific parameters — ligand is a compounds stream; the
        # residue code is resolved from its `code` column at runtime. A bare
        # string is shorthand for an internal Ligand(code=...).
        self.ligand_stream: DataStream = resolve_basic_input(
            ligand, Ligand, "compounds", "code", allow_none=False)
        self.ppi_design = ppi_design
        self.ppi_hotspot_residues = ppi_hotspot_residues or []
        self.ppi_binder_length = ppi_binder_length
        self.autogenerate_contigs = autogenerate_contigs
        self.model_only_neighbors = model_only_neighbors
        self.num_recycles = num_recycles
        self.scaffold_guided = scaffold_guided
        self.align_motif = align_motif
        self.deterministic = deterministic
        self.inpaint_str = inpaint_str
        self.inpaint_seq = inpaint_seq
        self.inpaint_length = inpaint_length
        # Accept a GuidingPotential builder (or list), or a raw hydra string for
        # back-compat; render to the single bracketed string the shell emits.
        self.guiding_potentials = render_guiding_potentials(guiding_potentials)
        self.guide_scale = guide_scale
        self.guide_decay = guide_decay
        # substrate: the 3-letter code for substrate_contacts. A bare string is the
        # code itself; otherwise it defaults at runtime to the same code as `ligand`
        # (the usual case — pull the design around the bound ligand). A Ligand stream
        # is accepted too and its code is resolved at runtime.
        self.substrate = substrate
        self.substrate_stream: Optional[DataStream] = None
        if substrate is not None and not isinstance(substrate, str):
            self.substrate_stream = resolve_basic_input(
                substrate, Ligand, "compounds", "code")

        # Normalize the per-PDB selection args (literal broadcast or table-column
        # reference). inpaint_str defaults to None here -> treat as unset.
        self._contigs_arg = _normalize_selection_arg("contigs", contigs)
        self._inpaint_arg = _normalize_selection_arg("inpaint", inpaint)
        self._inpaint_str_arg = _normalize_selection_arg(
            "inpaint_str", inpaint_str if inpaint_str is not None else "")

        if isinstance(contigs, str) and '/' in contigs:
            print("Warning: Character '/' found in contigs. RFdiffusionAllAtom uses ','.")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion-AllAtom parameters."""
        contigs_kind, contigs_token = self._contigs_arg
        if contigs_kind == "literal" and not contigs_token:
            raise ValueError("contigs parameter is required for RFdiffusion-AllAtom")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.steps <= 0:
            raise ValueError("steps must be positive")

        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")

        if self.ppi_design and not self.ppi_hotspot_residues:
            raise ValueError("PPI design requires hotspot residues")

        if self.ppi_design and self.ppi_binder_length is None:
            raise ValueError("PPI design requires binder length")

        if self.num_recycles < 1:
            raise ValueError("num_recycles must be at least 1")

        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand (a compounds stream, e.g. Ligand(code=...)) is required and must not be empty")

        # A column reference is keyed by input-PDB id, so it needs PDBs.
        # Only literal selections are interpolated raw into bash; table
        # references reach the shell through a resolved subshell.
        for name, (kind, token) in (("contigs", self._contigs_arg),
                                     ("inpaint", self._inpaint_arg),
                                     ("inpaint_str", self._inpaint_str_arg)):
            if kind == "table" and self.pdb_stream is None:
                raise ValueError(
                    f"{name} given as a table column reference requires a pdb input"
                )
            if kind == "literal":
                _validate_freeform_string(name, token)
        _validate_freeform_string("inpaint_seq", self.inpaint_seq)
        # guiding_potentials is a structured hydra value (e.g.
        # '["type:substrate_contacts,weight:3,...","type:monomer_ROG,..."]') that
        # legitimately needs commas, brackets, colons, and double quotes — so it
        # can't go through the generic freeform check. It is emitted single-quoted
        # into the shell, so we only forbid characters that would break single
        # quoting or trigger expansion despite it: a literal single quote, $, backtick,
        # and backslash.
        if self.guiding_potentials is not None:
            if not isinstance(self.guiding_potentials, str):
                raise ValueError("guiding_potentials must be a string")
            bad = set("'`$\\") & set(self.guiding_potentials)
            if bad:
                raise ValueError(
                    f"guiding_potentials contains {sorted(bad)}, which would break the "
                    "single-quoted shell argument. Use double quotes inside the list and "
                    "avoid ' ` $ \\."
                )
        _validate_freeform_string("guide_decay", self.guide_decay)
        if self.substrate is not None and isinstance(self.substrate, str):
            _validate_freeform_string("substrate", self.substrate)
        for i, res in enumerate(self.ppi_hotspot_residues):
            _validate_freeform_string(f"ppi_hotspot_residues[{i}]", res)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion-AllAtom configuration display lines."""
        config_lines = super().get_config_display()

        mode = "container" if self.uses_container() else "environment (SE3nv)"
        config_lines.extend([
            f"MODE: {mode}",
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])

        if self.pdb_stream:
            config_lines.append(f"PDB: {', '.join(self.pdb_stream.ids)}")
        if self._inpaint_arg[1]:
            config_lines.append(f"INPAINT: {self.inpaint}")
        if self.partial_steps > 0:
            config_lines.append(f"PARTIAL STEPS: {self.partial_steps}")
        if self.reproducible:
            config_lines.append(f"REPRODUCIBLE: {self.reproducible}")

        if self.ligand_stream:
            config_lines.append("LIGAND: (code resolved from compounds stream at runtime)")
        if self.ppi_design:
            config_lines.append(f"PPI DESIGN: {self.ppi_design}")
            config_lines.append(f"PPI HOTSPOTS: {','.join(self.ppi_hotspot_residues)}")
            config_lines.append(f"PPI BINDER LENGTH: {self.ppi_binder_length}")
        if self.autogenerate_contigs:
            config_lines.append(f"AUTOGENERATE CONTIGS: {self.autogenerate_contigs}")
        if self.num_recycles > 1:
            config_lines.append(f"NUM RECYCLES: {self.num_recycles}")
        if self.scaffold_guided:
            config_lines.append(f"SCAFFOLD GUIDED: {self.scaffold_guided}")
        if not self.align_motif:
            config_lines.append(f"ALIGN MOTIF: {self.align_motif}")
        if self.deterministic:
            config_lines.append(f"DETERMINISTIC: {self.deterministic}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RFdiffusion-AllAtom execution script."""
        # Serialize input DataStream to JSON for runtime file resolution.
        # configuration/ is auto-created by the pipeline.
        if self.pdb_stream:
            self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion-AllAtom execution script\n"
        script_content += self.generate_completion_check_header()
        # e3nn 0.3.3 uses torch.load() without weights_only=False,
        # which fails on PyTorch 2.6+ where the default flipped to True
        script_content += "export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1\n" 
        # activate_environment() activates the tool's configured env on the
        # host (SE3nv by default, or biopipelines fallback). Inference runs
        # under container_prefix when a container is set; host-side helpers
        # below run under the activated env either way.
        script_content += self.activate_environment()
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _build_common_inference_args(self) -> List[str]:
        """Inference args shared by every input PDB (no contigs/pdb/prefix/inpaint).

        The per-PDB args (contigs, input_pdb, output_prefix, inpaint_seq/str)
        are added inside the bash loop in _generate_script_run_rfdiffusion;
        everything here is constant across input structures.
        """
        aa_args = []
        aa_args.append("inference.ckpt_path=RFDiffusionAA_paper_weights.pt")
        aa_args.append(f"diffuser.T={self.steps}")
        aa_args.append(f"inference.num_designs={self.num_designs}")
        aa_args.append(f"inference.design_startnum={self.design_startnum}")
        # Ligand residue code resolved from the compounds stream at runtime
        # into $LIGAND_CODE (one ligand, broadcast to every PDB).
        aa_args.append("inference.ligand=$LIGAND_CODE")

        if self.ppi_design:
            aa_args.append(f"inference.ppi_design={self.ppi_design}")
        if self.ppi_hotspot_residues:
            aa_args.append(f"ppi.hotspot_res=[\\\'{','.join(self.ppi_hotspot_residues)}\\\']")
        if self.ppi_binder_length is not None:
            aa_args.append(f"ppi.binderlen={self.ppi_binder_length}")

        if self.autogenerate_contigs:
            aa_args.append(f"inference.autogenerate_contigs={self.autogenerate_contigs}")
        if self.model_only_neighbors:
            aa_args.append(f"inference.model_only_neighbors={self.model_only_neighbors}")
        if self.num_recycles > 1:
            aa_args.append(f"inference.num_recycles={self.num_recycles}")
        if self.scaffold_guided:
            aa_args.append(f"inference.scaffold_guided={self.scaffold_guided}")
        if not self.align_motif:
            aa_args.append(f"inference.align_motif={self.align_motif}")
        if self.deterministic:
            aa_args.append(f"inference.deterministic={self.deterministic}")

        # inpaint_seq / inpaint_length are broadcast scalars (distinct from the
        # per-PDB `inpaint` arg, which maps to contigmap.inpaint_seq in-loop).
        if self.inpaint_seq:
            aa_args.append(f"contigmap.inpaint_seq={self.inpaint_seq}")
        if self.inpaint_length is not None:
            aa_args.append(f"contigmap.length={self.inpaint_length}")

        if self.partial_steps > 0:
            aa_args.append(f"diffuser.partial_T={self.partial_steps}")

        if self.guiding_potentials:
            aa_args.append(f"potentials.guiding_potentials={self.guiding_potentials}")
        if self.guide_scale is not None:
            aa_args.append(f"potentials.guide_scale={self.guide_scale}")
        if self.guide_decay is not None:
            aa_args.append(f"potentials.guide_decay={self.guide_decay}")
        # substrate code for substrate_contacts: explicit string, resolved stream
        # code, or default to the ligand code ($LIGAND_CODE, set by ligand_snippet).
        # Only emitted when a substrate_contacts potential is actually in use.
        # potentials.substrate is NOT in RFdiffusionAA's base hydra config struct, so
        # it must be APPENDED with '+' (a plain override raises "Key 'substrate' is
        # not in struct"). Emit '+potentials.substrate=...'.
        if self.guiding_potentials and "substrate_contacts" in str(self.guiding_potentials):
            if isinstance(self.substrate, str):
                aa_args.append(f"+potentials.substrate={self.substrate}")
            elif self.substrate_stream is not None:
                self.substrate_stream.save_json(self.substrate_json)
                sub_id = Resolve.stream_ids(self.substrate_json, index=0)
                sub_code = Resolve.stream_item(self.substrate_json, sub_id, column='code')
                aa_args.append(f"+potentials.substrate={sub_code}")
            else:
                aa_args.append("+potentials.substrate=$LIGAND_CODE")

        return aa_args

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion-AllAtom execution part of the script.

        With a PDB input, runs once per input structure in a bash loop over the
        input ids (resolved at runtime). Per-PDB contigs / inpaint / inpaint_str
        are pre-resolved into contig_options.json by pipe_rfdiffusion_contigs.py.
        The single ligand is broadcast to every PDB (its code resolved once).
        Without a PDB it is a single unconditional run.
        """
        # Resolve the ligand `code` from the compounds stream at runtime: pick
        # the first id at runtime (not ids[0] at config time), then read its
        # `code` value. One ligand is reused across the whole loop.
        self.ligand_stream.save_json(self.ligand_json)
        ligand_snippet = f"""LIGAND_ID={Resolve.stream_ids(self.ligand_json, index=0)}
LIGAND_CODE={Resolve.stream_item(self.ligand_json, '$LIGAND_ID', column='code')}
"""

        common_args = self._build_common_inference_args()
        common = ' '.join(common_args)
        # Bash-array elements: double-quote each so $VAR expands (e.g. $LIGAND_CODE)
        # while literal " inside the value (guiding_potentials) is preserved as \".
        common_array = ' '.join('"' + a.replace('"', '\\"') + '"' for a in common_args)
        repo_dir = self.folders["RFdiffusionAllAtom"]
        mode = "container" if self.uses_container() else "environment"

        if not self.pdb_stream:
            # Unconditional — single run, contigs is a literal.
            contigs = self._contigs_arg[1]
            aa = f"contigmap.contigs=[\\\'{contigs}\\\'] inference.input_pdb=null"
            aa += f" inference.output_prefix={self.stream_path('structures', self.pipeline_name)}"
            if self._inpaint_arg[1]:
                aa += f" contigmap.inpaint_seq=[\\\'{self._inpaint_arg[1]}\\\']"
            if self._inpaint_str_arg[1]:
                aa += f" contigmap.inpaint_str={self._inpaint_str_arg[1]}"
            return f"""{ligand_snippet}echo "Starting RFdiffusion-AllAtom ({mode} mode, unconditional)"
echo "Output folder: {self.output_folder}"

cd {repo_dir}
{self.container_prefix()}python {self.inference_py_file} {aa} {common}

"""

        # PDB-conditioned: pre-resolve per-PDB selections, then loop.
        contigs_arg = _selection_cli_arg(self._contigs_arg)
        inpaint_arg = _selection_cli_arg(self._inpaint_arg)
        inpaint_str_arg = _selection_cli_arg(self._inpaint_str_arg)
        structures_dir = self.stream_folder("structures")

        with open(self.contig_args_json, "w") as f:
            json.dump({
                "structures_json": str(self.pdb_ds_json),
                "contigs": contigs_arg,
                "inpaint": inpaint_arg,
                "inpaint_str": inpaint_str_arg,
                "output_json": str(self.contigs_json),
            }, f, indent=2)

        return f"""{ligand_snippet}echo "Resolving per-PDB contig options"
python {self.contigs_py} "{self.contig_args_json}"

echo "Starting RFdiffusion-AllAtom ({mode} mode)"
echo "Output folder: {self.output_folder}"
cd {repo_dir}

for STRUCT_ID in {Resolve.stream_ids(self.pdb_ds_json)}; do
    INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, '$STRUCT_ID')}
    OPTS=$(python {self.contigs_reader_py} "{self.contigs_json}" "$STRUCT_ID")
    CONTIGS=$(echo "$OPTS" | sed -n '1p')
    INPAINT_SEL=$(echo "$OPTS" | sed -n '2p')
    INPAINT_STR_SEL=$(echo "$OPTS" | sed -n '3p')
    AA_OPTIONS=(
        "contigmap.contigs=['$CONTIGS']"
        "inference.input_pdb=$INPUT_PDB"
        "inference.output_prefix={structures_dir}/$STRUCT_ID"
        {common_array}
    )
    if [ -n "$INPAINT_SEL" ]; then
        AA_OPTIONS+=("contigmap.inpaint_seq=['$INPAINT_SEL']")
    fi
    if [ -n "$INPAINT_STR_SEL" ]; then
        AA_OPTIONS+=("contigmap.inpaint_str=$INPAINT_STR_SEL")
    fi
    echo "Running RFdiffusion-AllAtom for $STRUCT_ID"
    {self.container_prefix()}python {self.inference_py_file} "${{AA_OPTIONS[@]}}"
done

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        # pipe_rfdiffusion_table.py scans the structures/ stream folder for every
        # produced .pdb/.trb pair (one prefix per input PDB, or one for the
        # unconditional case).
        structures_dir = self.stream_folder("structures")
        return f"""echo "Creating results table"
python {self.table_py_file} "{structures_dir}" "{self.main_table}"

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # Each design id is "<pdb_id>_<n>"; derive the `structures.id` parent
        # provenance from that suffix at runtime (one or many parent PDBs).
        prov_arg = ' --provenance-from-suffix "structures.id"' if self.pdb_stream else ""
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion-AllAtom execution."""
        start = self.design_startnum
        end = self.design_startnum + self.num_designs - 1
        suffix_pattern = f"<{start}..{end}>"
        if self.pdb_stream:
            structure_ids = generate_multiplied_ids_pattern(
                self.pdb_stream.ids, suffix_pattern,
                input_stream_name="structures"
            )
        else:
            structure_ids = [f"{self.pipeline_name}_{suffix_pattern}"]
        file_template = [self.stream_path("structures", "<id>.pdb")]

        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map(); here we only declare the
        # stream and its map_table path.
        structures_map = self.stream_map_path("structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=file_template,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "pdb", "fixed", "designed", "source_fixed", "plddt_mean", "status"],
                description="RFdiffusion-AllAtom structure generation results with fixed/designed regions"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all RFdiffusion-AllAtom parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfdaa_params": {
                "pdb_input_ids": self.pdb_stream.ids if self.pdb_stream else None,
                "contigs": str(self.contigs),
                "inpaint": str(self.inpaint),
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum,
                "ligand_ids": list(self.ligand_stream.ids),
                "ppi_design": self.ppi_design,
                "ppi_hotspot_residues": self.ppi_hotspot_residues,
                "ppi_binder_length": self.ppi_binder_length,
                "autogenerate_contigs": self.autogenerate_contigs,
                "model_only_neighbors": self.model_only_neighbors,
                "num_recycles": self.num_recycles,
                "scaffold_guided": self.scaffold_guided,
                "align_motif": self.align_motif,
                "deterministic": self.deterministic,
                "inpaint_str": str(self.inpaint_str) if self.inpaint_str is not None else None,
                "inpaint_seq": self.inpaint_seq,
                "inpaint_length": self.inpaint_length,
                "guiding_potentials": self.guiding_potentials,
                "guide_scale": self.guide_scale,
                "guide_decay": self.guide_decay,
                "substrate": self.substrate if isinstance(self.substrate, str) else None,
            }
        })
        return base_dict


class RFDAA_PrepareLigand(BaseConfig):
    """
    Preparation tool for RFdiffusion-AllAtom to add a dummy peptide to ligand-only PDB files.
    """

    TOOL_NAME = "RFDAA_PrepareLigand"
    TOOL_VERSION = "1.1"

    # Lazy path descriptors
    #   prepared_pdb    — the single output PDB, lives in structures/.
    #   structures_csv  — standalone TableInfo (tables/structures.csv).
    #   ligand_ds_json  — config-time input DataStream serialization.
    prepared_pdb = Path(lambda self: self.stream_path("structures", "prepared_ligand.pdb"))
    structures_csv = Path(lambda self: self.table_path("structures"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_rfdaa_prepare_ligand.py"))
    ligand_ds_json = Path(lambda self: self.configuration_path("input_ligand.json"))

    def __init__(self,
                 ligand: Union[DataStream, StandardizedOutput],
                 **kwargs):
        """
        Initialize RFDAA_PrepareLigand tool.

        Args:
            ligand: Ligand structure as DataStream or StandardizedOutput
            **kwargs: Additional parameters

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | file_path
        """
        # Resolve ligand input — store stream for runtime resolution
        if isinstance(ligand, StandardizedOutput):
            self.ligand_stream = ligand.streams.structures
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")
        # The ligand file is resolved at runtime via Resolve.stream_ids(index=0)
        # in generate_script — never index ids[0] at config time (lazy IDs).

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate parameters."""
        if not self.ligand_stream:
            raise ValueError("ligand input is required")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display."""
        config_lines = super().get_config_display()
        config_lines.append(f"LIGAND_SOURCE: {', '.join(self.ligand_stream.ids)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to combine ligand with dummy peptide."""
        # configuration/ + structures/ folders auto-created by the pipeline.
        self.ligand_stream.save_json(self.ligand_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFDAA_PrepareLigand execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""LIGAND_ID={Resolve.stream_ids(self.ligand_ds_json, index=0, valid_set=True)}
LIGAND_FILE={Resolve.stream_item(self.ligand_ds_json, '$LIGAND_ID')}

echo "Preparing ligand structure for RFdiffusion-AllAtom"
echo "Input ligand: $LIGAND_FILE"
echo "Output: {self.prepared_pdb}"

python "{self.helper_script}" \\
  --ligand_pdb "$LIGAND_FILE" \\
  --output_pdb "{self.prepared_pdb}" \\
  --output_csv "{self.structures_csv}" \\
  --structures_map "{self.structures_map}" \\
  --pdbs_folder "{self.folders['pdbs']}"

if [ $? -eq 0 ]; then
    echo "Successfully prepared ligand structure"
else
    echo "Error: Failed to prepare ligand structure"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        structure_ids = ["prepared_ligand"]
        structure_files = [self.prepared_pdb]

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "file_path"],
                description="Prepared ligand structure with dummy peptide"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "ligand_input_ids": list(self.ligand_stream.ids)
            }
        })
        return base_dict
