# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion2 configuration for enzyme active-site scaffolding and small-molecule binder design.

Wraps RosettaCommons/RFdiffusion2 (Apptainer container, Hydra entrypoint
``rf_diffusion/run_inference.py --config-name=aa``). Its defining capability is
atomic-level active-site specification: a motif is given as individual side-chain
atoms (``contigmap.contig_atoms``) that the model scaffolds either at fixed
sequence positions (indexed) or anywhere in the chain (unindexed guidepost).
Also supports RASA-conditioned small-molecule binder design.

Distinct from RFdiffusionAllAtom (the earlier all-atom diffusion model) and
RFdiffusion3 (the foundry/rfd3 successor).
"""

import json
import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
    from .combinatorics import generate_multiplied_ids_pattern
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve, TableReference
    from combinatorics import generate_multiplied_ids_pattern
    from input_standardization import resolve_basic_input
    from ligand import Ligand


def _normalize_contigs_arg(value):
    """Normalize the contigs argument to a (kind, token) pair.

    ``contigs`` accepts a plain string (broadcast to every input PDB) or a
    per-PDB column reference (a ``TableReference`` from ``tool.tables.X.col``, or
    a ``(TableInfo, "col")`` tuple). Returns ``("literal", str)`` or
    ``("table", "TABLE_REFERENCE:path:col")``.
    """
    if isinstance(value, TableReference):
        return ("table", str(value))
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return ("table", str(TableReference(value[0].info.path, value[1])))
    if isinstance(value, str):
        return ("literal", value)
    raise ValueError(
        f"contigs must be a string or a (table, column) reference, got {type(value)}"
    )


def _hydra_contig_atoms(contig_atoms: Dict[str, str]) -> str:
    """Render contig_atoms dict into the Hydra override RFdiffusion2 expects.

    Hydra cannot lex a bare dict literal with quoted keys ({'A106':...}); the
    value must be a *quoted string* that RFdiffusion2 parses itself. Verified
    against the container's Hydra OverridesParser: the double-quoted form keeps
    it as the string "{'A106':'NE,CD,CZ'}". Emitted into a double-quoted
    bash-array element where the inner " is backslash-escaped so it survives to
    Hydra. Hydra receives, e.g.:
        contigmap.contig_atoms="{'A106':'NE,CD,CZ','A166':'OD1,CG'}"
    """
    body = ",".join(f"'{res}':'{atoms}'" for res, atoms in contig_atoms.items())
    return f'contigmap.contig_atoms="{{{body}}}"'


def _hydra_partially_fixed_ligand(partially_fixed_ligand: Dict[str, List[str]]) -> str:
    """Render partially_fixed_ligand dict into an appended Hydra override.

    ``inference.partially_fixed_ligand`` is not in the base config struct, so it
    must be appended with ``++``. A dict of code -> atom-list; emitted into a
    double-quoted bash-array element so no shell escaping is needed. Hydra
    receives, e.g.:
        ++inference.partially_fixed_ligand={NAD:[O7N,C7N],OXM:[O3,C2]}
    """
    body = ",".join(
        f"{code}:[{','.join(atoms)}]" for code, atoms in partially_fixed_ligand.items()
    )
    return f"++inference.partially_fixed_ligand={{{body}}}"


class RFdiffusion2(BaseConfig):
    """
    RFdiffusion2 for atomic active-site scaffolding and small-molecule binder design.

    Two modes, selected by the inputs you give:

    - **Active-site scaffolding.** Provide ``pdb`` (the theozyme / active-site
      structure), ``ligand`` (one or more bound HETATM codes), ``contigs``, and
      ``contig_atoms`` (the per-residue side-chain atoms that define the motif).
      ``contig_as_guidepost=True`` scaffolds the motif unindexed (the model
      chooses where it lands); ``False`` keeps the indexed positions.

    - **Small-molecule binder design.** Provide ``pdb`` (ligand-only or
      ligand+origin), ``ligand``, ``contigs`` (a length, e.g. ``"150"``),
      ``length``, and ``relative_sasa`` to bury the design around the ligand.
    """

    TOOL_NAME = "RFdiffusion2"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("RFdiffusion2", "")
        parent_dir = os.path.dirname(repo_dir)
        sif = f"{repo_dir}/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif"
        # The install artifacts setup.py fetches: the container image AND the
        # inference weights. setup.py swallows every download error and still
        # prints success with exit 0, so its exit code is meaningless — verify
        # the real files. RFD_173.pt is the demo default ckpt; the .sif is what
        # inference runs in.
        weights = f"{repo_dir}/rf_diffusion/model_weights/RFD_173.pt"
        # A weight present and non-empty proves the proxy reached files.ipd.uw.edu;
        # a 0-byte / missing file is the "Connection refused" failure mode.
        verify = f'[ -f "{sif}" ] && [ -s "{weights}" ]'
        skip = "" if force_reinstall else f"""# Check if already installed (image + weights, both non-empty)
if [ -d "{repo_dir}" ] && {verify}; then
    echo "RFdiffusion2 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        # setup.py skips files already present, so a re-run only fetches the
        # weights/image still missing (e.g. after a proxy blip) without re-pulling
        # the multi-GB .sif. Use 'overwrite' manually to force a clean refetch.
        return f"""echo "=== Installing RFdiffusion2 ==="
{skip}mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/RosettaCommons/RFdiffusion2.git
fi
cd {repo_dir}
export PYTHONPATH="{repo_dir}"
# Weights download via urllib needs the outbound proxy. scheduler.init exports
# it on compute nodes; re-source here so a direct .install() shell also has it.
for f in /etc/profile.d/proxy.sh; do [ -f "$f" ] && source "$f"; done
# Downloads model weights + the Apptainer container (.sif). Large; can take 30+
# minutes. setup.py reports success even when downloads fail — we verify below.
python setup.py

# Verify the REAL artifacts (setup.py's own success message is unreliable).
if {verify}; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion2 installation complete ==="
    echo "Container mode: set containers.RFdiffusion2 to {sif} in config.yaml"
else
    echo "ERROR: RFdiffusion2 verification failed."
    [ -f "{sif}" ] || echo "       Container image missing: {sif}"
    [ -s "{weights}" ] || echo "       Model weights missing/empty: {weights} (download likely blocked — check the outbound proxy on this node)"
    echo "       Re-run install once the node can reach files.ipd.uw.edu."
    exit 1
fi
"""

    # Lazy path descriptors
    #   main_table  — standalone TableInfo CSV (tables/structures.csv).
    #   pdb_ds_json — config-time input DataStream serialization.
    main_table = Path(lambda self: self.table_path("structures"))
    inference_py_file = Path(lambda self: "rf_diffusion/run_inference.py")
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_table.py"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))
    # Per-PDB contig resolution (shared with the RFdiffusion family).
    contigs_json = Path(lambda self: self.configuration_path("contig_options.json"))
    contigs_py = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_contigs.py"))
    contigs_reader_py = Path(lambda self: self.pipe_script_path("resolve_rfdiffusion_contigs.py"))

    def __init__(self,
                 ligand: Union[str, DataStream, StandardizedOutput],
                 pdb: Union[DataStream, StandardizedOutput],
                 contigs: Union[str, TableReference, "tuple"],
                 contig_atoms: Optional[Dict[str, str]] = None,
                 contig_as_guidepost: bool = True,
                 only_guidepost_positions: Optional[str] = None,
                 partially_fixed_ligand: Optional[Dict[str, List[str]]] = None,
                 length: Optional[str] = None,
                 relative_sasa: Optional[float] = None,
                 num_designs: int = 1,
                 design_startnum: int = 1,
                 deterministic: bool = False,
                 seed_offset: Optional[int] = None,
                 **kwargs):
        """
        Initialize RFdiffusion2 configuration.

        Args:
            ligand: Ligand(s) as a compounds stream. The 3-letter residue ``code``
                of every id in the stream is read at runtime and joined into
                ``inference.ligand='CODE1,CODE2'``. For multiple bound ligands pass
                ONE stream carrying all codes, e.g. ``Ligand(code=["NAD", "OXM"])``
                (or any compounds-producing tool). A bare string is shorthand for an
                internal ``Ligand(code=...)``.
            pdb: Input PDB structure (the active-site/theozyme or the ligand-only
                binder target) as DataStream or StandardizedOutput. Iterated at
                execution time (one run per input structure).
            contigs: Contig specification (RFdiffusion2 uses ``,`` separators,
                e.g. ``"46,A106-106,59,A166-166,2,A169-169,23,A193-193,46"``, or a
                bare length like ``"150"`` for binder design). A string is
                broadcast to every input PDB; a ``(table, column)`` reference is
                resolved per input-PDB id at runtime.
            contig_atoms: Per-residue side-chain atoms defining the atomic motif,
                e.g. ``{"A106": "NE,CD,CZ", "A166": "OD1,CG"}`` -> ``contigmap.contig_atoms``.
                Broadcast to every input PDB.
            contig_as_guidepost: True (default) scaffolds the motif unindexed (the
                model places it); False keeps the indexed contig positions.
                Maps to ``inference.contig_as_guidepost``.
            only_guidepost_positions: Restrict guidepost treatment to these
                positions (e.g. ``"A106"``); the rest stay indexed. Maps to
                ``inference.only_guidepost_positions``.
            partially_fixed_ligand: Per-ligand subset of atoms to keep fixed,
                e.g. ``{"NAD": ["O7N", "C7N"], "OXM": ["O3", "C2"]}`` ->
                ``++inference.partially_fixed_ligand`` (appended; not in base struct).
            length: Total design length for binder mode (e.g. ``"150-150"``) ->
                ``contigmap.length``.
            relative_sasa: Relative-SASA target for RASA-conditioned binder design
                (e.g. ``0`` to fully bury). Enables
                ``inference.conditions.relative_sasa_v2`` (active + rasa).
            num_designs: Number of designs per input PDB -> ``inference.num_designs``.
            design_startnum: Starting design number -> ``inference.design_startnum``.
            deterministic: Fixed RNG for reproducible outputs -> ``inference.deterministic``.
            seed_offset: Seed offset when deterministic -> ``inference.seed_offset``.
            **kwargs: Additional parameters passed to BaseConfig.

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | pdb | fixed | designed | source_fixed | plddt_mean | status
        """
        # Input PDB — required; iterated at execution time (never index ids[0] at
        # config time, which breaks under lazy IDs).
        if isinstance(pdb, StandardizedOutput):
            self.pdb_stream = pdb.streams.structures
        elif isinstance(pdb, DataStream):
            self.pdb_stream = pdb
        else:
            raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Ligand — compounds stream; every id's `code` is joined at runtime. A
        # bare string is shorthand for an internal Ligand(code=...).
        self.ligand_stream: DataStream = resolve_basic_input(
            ligand, Ligand, "compounds", "code", allow_none=False)

        self.contigs = contigs
        self.contig_atoms = contig_atoms
        self.contig_as_guidepost = contig_as_guidepost
        self.only_guidepost_positions = only_guidepost_positions
        self.partially_fixed_ligand = partially_fixed_ligand
        self.length = length
        self.relative_sasa = relative_sasa
        self.num_designs = num_designs
        self.design_startnum = design_startnum
        self.deterministic = deterministic
        self.seed_offset = seed_offset

        self._contigs_arg = _normalize_contigs_arg(contigs)

        if isinstance(contigs, str) and '/' in contigs:
            print("Warning: Character '/' found in contigs. RFdiffusion2 uses ','.")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion2 parameters."""
        contigs_kind, contigs_token = self._contigs_arg
        if contigs_kind == "literal" and not contigs_token:
            raise ValueError("contigs parameter is required for RFdiffusion2")
        if contigs_kind == "literal":
            _validate_freeform_string("contigs", contigs_token)

        if not self.pdb_stream or len(self.pdb_stream) == 0:
            raise ValueError("pdb input is required and must not be empty")

        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand (a compounds stream, e.g. Ligand(code=...)) is required and must not be empty")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.contig_atoms is not None:
            if not isinstance(self.contig_atoms, dict):
                raise ValueError("contig_atoms must be a dict of {residue: 'atom,atom,...'}")
            for res, atoms in self.contig_atoms.items():
                _validate_freeform_string("contig_atoms key", res)
                _validate_freeform_string("contig_atoms value", atoms)

        if self.partially_fixed_ligand is not None:
            if not isinstance(self.partially_fixed_ligand, dict):
                raise ValueError("partially_fixed_ligand must be a dict of {code: [atom, ...]}")
            for code, atoms in self.partially_fixed_ligand.items():
                _validate_freeform_string("partially_fixed_ligand key", code)
                if not isinstance(atoms, (list, tuple)):
                    raise ValueError("partially_fixed_ligand values must be lists of atom names")
                for atom in atoms:
                    _validate_freeform_string("partially_fixed_ligand atom", atom)

        if self.only_guidepost_positions is not None:
            _validate_freeform_string("only_guidepost_positions", self.only_guidepost_positions)
        if self.length is not None:
            _validate_freeform_string("length", self.length)
        if self.relative_sasa is not None and not isinstance(self.relative_sasa, (int, float)):
            raise ValueError("relative_sasa must be a number")
        if self.seed_offset is not None and not isinstance(self.seed_offset, int):
            raise ValueError("seed_offset must be an integer")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion2 configuration display lines."""
        config_lines = super().get_config_display()

        mode = "container" if self.uses_container() else "environment (SE3nv)"
        config_lines.extend([
            f"MODE: {mode}",
            f"PDB: {', '.join(self.pdb_stream.ids)}",
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            "LIGAND: (codes resolved from compounds stream at runtime)",
            f"CONTIG AS GUIDEPOST: {self.contig_as_guidepost}",
        ])
        if self.contig_atoms:
            config_lines.append(f"CONTIG ATOMS: {self.contig_atoms}")
        if self.only_guidepost_positions:
            config_lines.append(f"ONLY GUIDEPOST POSITIONS: {self.only_guidepost_positions}")
        if self.partially_fixed_ligand:
            config_lines.append(f"PARTIALLY FIXED LIGAND: {self.partially_fixed_ligand}")
        if self.length:
            config_lines.append(f"LENGTH: {self.length}")
        if self.relative_sasa is not None:
            config_lines.append(f"RELATIVE SASA: {self.relative_sasa}")
        if self.deterministic:
            config_lines.append(f"DETERMINISTIC: {self.deterministic}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RFdiffusion2 execution script."""
        self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion2 execution script\n"
        script_content += self.generate_completion_check_header()
        # e3nn / torch.load weights_only default flipped in torch 2.6+; the shared
        # RFdiffusion-AA stack relies on the legacy behavior.
        script_content += "export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1\n"
        script_content += self.activate_environment()
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _build_common_inference_args(self) -> List[str]:
        """Inference args shared by every input PDB (no contigs/pdb/prefix).

        Per-PDB args (input_pdb, output_prefix, contigs) are added in the loop;
        everything here is constant across input structures. The ligand codes are
        broadcast to every PDB via the runtime-resolved $LIGAND_CODES.
        """
        args = []
        # Single-quote the code list: unquoted 'NAD,OXM' is read by Hydra as a
        # ChoiceSweep (multirun), not the comma-string RFdiffusion2 expects. The
        # quotes are literal inside the double-quoted bash-array element.
        args.append("inference.ligand='$LIGAND_CODES'")
        args.append(f"inference.num_designs={self.num_designs}")
        args.append(f"inference.design_startnum={self.design_startnum}")
        args.append(f"inference.contig_as_guidepost={self.contig_as_guidepost}")

        if self.contig_atoms:
            args.append(_hydra_contig_atoms(self.contig_atoms))
        if self.only_guidepost_positions:
            args.append(f"inference.only_guidepost_positions={self.only_guidepost_positions}")
        if self.partially_fixed_ligand:
            args.append(_hydra_partially_fixed_ligand(self.partially_fixed_ligand))
        if self.length:
            args.append(f"contigmap.length={self.length}")
        if self.relative_sasa is not None:
            args.append("inference.conditions.relative_sasa_v2.active=True")
            args.append(f"inference.conditions.relative_sasa_v2.rasa={self.relative_sasa}")
        if self.deterministic:
            args.append(f"inference.deterministic={self.deterministic}")
        if self.seed_offset is not None:
            args.append(f"inference.seed_offset={self.seed_offset}")

        return args

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion2 execution part of the script.

        Runs once per input PDB in a bash loop over the input ids (resolved at
        runtime). Per-PDB contigs are pre-resolved into contig_options.json by
        pipe_rfdiffusion_contigs.py. The ligand codes are broadcast to every PDB.
        """
        # Resolve every ligand `code` from the compounds stream at runtime and
        # join with ',' into $LIGAND_CODES (RFdiffusion2 accepts multiple codes).
        self.ligand_stream.save_json(self.ligand_json)
        ligand_snippet = f"""LIGAND_CODES=""
for LIG_ID in {Resolve.stream_ids(self.ligand_json)}; do
    LIG_CODE={Resolve.stream_item(self.ligand_json, '$LIG_ID', column='code')}
    if [ -z "$LIGAND_CODES" ]; then
        LIGAND_CODES="$LIG_CODE"
    else
        LIGAND_CODES="$LIGAND_CODES,$LIG_CODE"
    fi
done
"""

        common_args = self._build_common_inference_args()
        # Bash-array elements: double-quote each so $VAR expands (e.g. $LIGAND_CODES)
        # while literal " inside the value (contig_atoms) is preserved as \".
        common_array = ' '.join('"' + a.replace('"', '\\"') + '"' for a in common_args)
        repo_dir = self.folders["RFdiffusion2"]
        mode = "container" if self.uses_container() else "environment"

        contigs_arg = self._contigs_arg[1] if self._contigs_arg[0] == "literal" else str(self._contigs_arg[1])
        # The shared contigs resolver takes contigs/inpaint/inpaint_str; RFd2 uses
        # only contigs, so pass '-' for the other two.
        contigs_cli = contigs_arg if contigs_arg else "-"
        structures_dir = self.stream_folder("structures")

        return f"""{ligand_snippet}echo "Resolving per-PDB contig options"
python {self.contigs_py} "{self.pdb_ds_json}" "{contigs_cli}" "-" "-" "{self.contigs_json}"

echo "Starting RFdiffusion2 ({mode} mode)"
echo "Output folder: {self.output_folder}"
cd {repo_dir}
export PYTHONPATH="{repo_dir}"

for STRUCT_ID in {Resolve.stream_ids(self.pdb_ds_json)}; do
    INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, '$STRUCT_ID')}
    OPTS=$(python {self.contigs_reader_py} "{self.contigs_json}" "$STRUCT_ID")
    CONTIGS=$(echo "$OPTS" | sed -n '1p')
    RF2_OPTIONS=(
        "--config-name=aa"
        "contigmap.contigs=['$CONTIGS']"
        "inference.input_pdb=$INPUT_PDB"
        "inference.output_prefix={structures_dir}/$STRUCT_ID"
        {common_array}
    )
    echo "Running RFdiffusion2 for $STRUCT_ID"
    {self.container_prefix()}python {self.inference_py_file} "${{RF2_OPTIONS[@]}}"
done

# RFdiffusion2 names outputs "<id>_<n>-atomized-bb-<bool>.pdb/.trb"; strip the
# "-atomized-bb-<bool>" suffix so files match the canonical "<id>_<n>" stream
# ids. Top-level only — the unidealized/ subdir copies are left as extras.
echo "Normalizing RFdiffusion2 output filenames"
for f in "{structures_dir}"/*-atomized-bb-*.pdb "{structures_dir}"/*-atomized-bb-*.trb; do
    [ -e "$f" ] || continue
    base=$(basename "$f")
    ext="${{base##*.}}"
    stem="${{base%.*}}"
    clean="${{stem%%-atomized-bb-*}}"
    mv -f "$f" "{structures_dir}/$clean.$ext"
done

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script (reuses RFdiffusion family helper)."""
        structures_dir = self.stream_folder("structures")
        return f"""echo "Creating results table"
python {self.table_py_file} "{structures_dir}" "{self.main_table}"

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # Each design id is "<pdb_id>_<n>"; derive the `structures.id` parent
        # provenance from that suffix at runtime.
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}" --provenance-from-suffix "structures.id"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion2 execution."""
        start = self.design_startnum
        end = self.design_startnum + self.num_designs - 1
        suffix_pattern = f"<{start}..{end}>"
        structure_ids = generate_multiplied_ids_pattern(
            self.pdb_stream.ids, suffix_pattern,
            input_stream_name="structures"
        )
        file_template = [self.stream_path("structures", "<id>.pdb")]

        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map(); here we only declare it.
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
                description="RFdiffusion2 structure generation results with fixed/designed regions"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }
