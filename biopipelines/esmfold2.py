# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold2 tool for all-atom biomolecular complex prediction.

ESMFold2 (EvolutionaryScale, MIT-licensed) predicts the all-atom structure of a
biomolecular complex from sequence, building on ESMC representations through a
recurrent folding stack and a diffusion module. It folds protein / DNA / RNA /
ligand complexes (including antibody scFvs as multi-chain proteins) and reports
per-atom pLDDT, pTM, interface pTM (ipTM) and, optionally, PAE. Without an MSA
in the single-sequence regime, or recycling precomputed MSAs per protein chain.

Inference-time scaling is a first-class control: increasing the number of
recurrent folding loops, diffusion samples, or independent seeds raises
prediction quality, with best-of-N selection by ipTM (complexes) or pLDDT
(monomers).

Reference:
    Candido, Hayes, Derry, Rao, Lin et al. "Language Modeling Materializes a
    World Model of Protein Biology" (bioRxiv 2026).
    Code/weights: https://github.com/evolutionaryscale/esm (biohub/ESMFold2).
"""

import json
import os
import re
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import (
        generate_combinatorics_config, get_mode,
        predict_output_ids_with_provenance, Bundle, Each,
    )
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import (
        generate_combinatorics_config, get_mode,
        predict_output_ids_with_provenance, Bundle, Each,
    )
    from datastream_resolver import resolve_input_to_datastream


_MODIFICATION_KEYS = {"chain", "position", "ccd"}
_BOND_KEYS = {"atom1", "atom2"}
# Same shape as ligand._CCD_CODE_RE; kept local to avoid importing that module.
_CCD_CODE_RE = re.compile(r'^[A-Za-z0-9]{1,5}$')


def _validate_chain_position(label, chain, position):
    if not isinstance(chain, str) or not chain:
        raise ValueError(f"{label}: chain must be a non-empty chain id (e.g. 'A')")
    # bool is an int subclass; True would otherwise pass as position 1.
    if isinstance(position, bool) or not isinstance(position, int):
        raise ValueError(f"{label}: position must be an int, got {type(position).__name__}")
    if position < 1:
        raise ValueError(f"{label}: positions are 1-indexed, got {position}")


def _check_keys(label, obj, expected):
    if not isinstance(obj, dict):
        raise ValueError(f"{label} must be a dict with keys {sorted(expected)}")
    unknown = set(obj) - expected
    if unknown:
        raise ValueError(f"{label}: unknown key(s) {sorted(unknown)}; expected {sorted(expected)}")
    missing = expected - set(obj)
    if missing:
        raise ValueError(f"{label}: missing key(s) {sorted(missing)}")


def _validate_modifications(modifications):
    if not isinstance(modifications, (list, tuple)) or not modifications:
        raise ValueError(
            "modifications must be a non-empty list of {'chain', 'position', 'ccd'} dicts")
    for i, mod in enumerate(modifications):
        label = f"modifications[{i}]"
        if isinstance(mod, dict) and "smiles" in mod:
            raise ValueError(
                f"{label}: modified residues are specified by CCD code only "
                "(upstream Modification.smiles is not implemented); use e.g. ccd='SEP'")
        _check_keys(label, mod, _MODIFICATION_KEYS)
        _validate_chain_position(label, mod["chain"], mod["position"])
        if not isinstance(mod["ccd"], str) or not mod["ccd"]:
            raise ValueError(f"{label}: ccd must be a non-empty CCD code (e.g. 'SEP')")


def _validate_glycans(glycans):
    """Each glycan is one ligand chain given as an ordered list of CCD codes.

    Upstream puts every code in its own residue of that chain but creates no
    bonds between them -- ligand_bonds is populated only for SMILES ligands.
    So the glycosidic linkages, and the bond to the protein, must be declared
    through covalent_bonds. Declaring any bond on the chain is also what makes
    upstream strip the leaving atoms (the anomeric O1), so an undeclared
    linkage is not merely unbonded, it is chemically wrong.
    """
    if not isinstance(glycans, (list, tuple)) or not glycans:
        raise ValueError(
            "glycans must be a non-empty list of CCD-code lists, "
            "e.g. [['NAG', 'NAG', 'BMA']]")
    for i, chain in enumerate(glycans):
        label = f"glycans[{i}]"
        if not isinstance(chain, (list, tuple)) or not chain:
            raise ValueError(f"{label} must be a non-empty list of CCD codes")
        for j, code in enumerate(chain):
            if not isinstance(code, str) or not code:
                raise ValueError(f"{label}[{j}] must be a non-empty CCD code")
            if not _CCD_CODE_RE.match(code):
                raise ValueError(
                    f"{label}[{j}]: invalid CCD code {code!r}; "
                    "expected 1-5 alphanumeric characters")


def _validate_covalent_bonds(bonds):
    if not isinstance(bonds, (list, tuple)) or not bonds:
        raise ValueError(
            "covalent_bonds must be a non-empty list of {'atom1', 'atom2'} dicts")
    for i, bond in enumerate(bonds):
        label = f"covalent_bonds[{i}]"
        _check_keys(label, bond, _BOND_KEYS)
        for key in ("atom1", "atom2"):
            atom = bond[key]
            if not isinstance(atom, (list, tuple)) or len(atom) != 3:
                raise ValueError(f"{label}.{key} must be [chain, position, atom_name]")
            _validate_chain_position(f"{label}.{key}", atom[0], atom[1])
            if not isinstance(atom[2], str) or not atom[2]:
                raise ValueError(f"{label}.{key}: atom_name must be a non-empty string")


class ESMFold2(BaseConfig):
    """
    ESMFold2: all-atom biomolecular complex prediction from sequence.

    Folds a complex assembled from one or more input axes (proteins, ssDNA,
    dsDNA, ssRNA, dsRNA, ligands), with Bundle/Each combinatorics matching the
    rest of the framework. One mmCIF per output complex, plus a confidence
    table (pLDDT, pTM, ipTM, optional max PAE).

    Inputs:
        proteins:   amino acid sequences (StandardizedOutput or DataStream)
        ligands:    a compounds stream (Ligand / CompoundLibrary); folded by CCD code
        ssDNA/dsDNA/ssRNA/dsRNA: nucleic-acid sequence streams (double-stranded
                    axes add the reverse-complement chain automatically)
        msas:       precomputed MSAs to recycle per protein chain (StandardizedOutput
                    with an msas table, e.g. from MMseqs2 or a previous run). When
                    omitted, ESMFold2 runs single-sequence (no MSA).
        modifications: per-residue CCD substitutions on a polymer chain, as
                    [{"chain": "A", "position": 122, "ccd": "SEP"}]. The residue is
                    replaced by the CCD entry and atom-tokenized, so SEP/TPO/PTR give
                    real phospho-residues rather than their unmodified parents. Also
                    applies to DNA/RNA chains (e.g. PSU). CCD codes only.
        covalent_bonds: inter-entity covalent links, as
                    [{"atom1": ["A", 122, "SG"], "atom2": ["C", 1, "C1"]}]. Atom names
                    are resolved against the residue's CCD atoms at runtime.
        num_loops:  recurrent folding loops (inference-time scaling)
        num_sampling_steps: diffusion sampling steps
        num_diffusion_samples: diffusion samples per seed
        num_seeds:  independent seeds; best by ipTM (complex) / pLDDT (monomer) is kept
        msa_max_depth: MSA subsampling depth per loop
        include_pae: also predict and tabulate PAE (max PAE in the confidence table)
        top_only:   keep only the best sample per complex (True) or surface every
                    diffusion sample as <id>_1..N (False)
        model_name: HuggingFace checkpoint id (default biohub/ESMFold2)

    Chain ids for modifications and covalent_bonds are assigned in axis order —
    proteins, then ssDNA/dsDNA/ssRNA/dsRNA, then ligands — starting at A, with a
    double-stranded axis taking two consecutive letters. They address a position
    in the complex, not a sequence id, so when an axis iterates the same
    constraint applies to every item of that axis.

    There is deliberately no `pocket` parameter. Upstream defines
    PocketConditioning and stores it on StructurePredictionInput, but
    build_feature_tensors sets pocket_feature to zeros unconditionally and never
    reads it (esm 3.3.0), so the constraint would be accepted and discarded. Do
    not add one back without checking that path first.

    Outputs:
        Streams:
            structures: predicted mmCIF, one per complex (B-factor = per-atom pLDDT)
        Tables:
            structures: id | file
            confidence: id | file | plddt | ptm | iptm [| max_pae]
            compounds:  ligand chemistry passthrough (code = the code ESMFold2 used)

    Usage::

        with Pipeline(project="Examples", job="ESMFold2-demo"):
            Resources(gpu="A100", memory="64GB", time="6:00:00")
            rec = Sequence(["MIEIK..."], ids=["EGFR"])
            binders = ProteinMPNN(structures=rfd, num_sequences=8)
            cplx = ESMFold2(proteins=Each(rec, binders), num_seeds=4)

        # phosphoserine at S122 of chain A, folded against a phosphatase
        ESMFold2(proteins=Each(substrate, phosphatase),
                 modifications=[{"chain": "A", "position": 122, "ccd": "SEP"}])
    """

    TOOL_NAME = "ESMFold2"
    TOOL_VERSION = "2.1"
    ENV_NAME = "esmfold2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the esmfold2 env (picks esmfold2.<variant>.yaml automatically) and
        verify by importing the transformers ESMFold2 model class."""
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "{env} environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing ESMFold2 ==="
{skip}{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create {env} environment."
    exit 1
fi

# Verify installation
if {cls._env_run(env, env_manager)}python -c "from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ESMFold2 installation complete ==="
else
    echo "ERROR: ESMFold2 verification failed (cannot import ESMFold2Model)"
    exit 1
fi
"""

    # ---------------------------------------------------------------------------
    # Lazy path descriptors
    # ---------------------------------------------------------------------------

    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    constraints_config_file   = Path(lambda self: self.configuration_path("constraints_config.json"))
    msas_table_file           = Path(lambda self: self.configuration_path(".input_msas.csv"))
    predictions_folder        = Path(lambda self: self.stream_folder("structures"))
    structures_map_csv        = Path(lambda self: self.stream_map_path("structures"))
    compounds_map_csv         = Path(lambda self: self.stream_map_path("compounds"))
    confidence_csv            = Path(lambda self: self.table_path("confidence"))
    missing_csv               = Path(lambda self: self.table_path("missing"))

    inference_py = Path(lambda self: self.pipe_script_path("pipe_esmfold2_inference.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_esmfold2_postprocessing.py"))

    # ---------------------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------------------

    def __init__(
        self,
        proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
        ssDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
        dsDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
        ssRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
        dsRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
        ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
        msas: Optional[StandardizedOutput] = None,
        modifications: Optional[List[Dict[str, Any]]] = None,
        covalent_bonds: Optional[List[Dict[str, Any]]] = None,
        glycans: Optional[List[List[str]]] = None,
        num_loops: int = 10,
        num_sampling_steps: int = 100,
        num_diffusion_samples: int = 1,
        num_seeds: int = 1,
        msa_max_depth: int = 1024,
        include_pae: bool = False,
        top_only: bool = True,
        model_name: str = "biohub/ESMFold2",
        **kwargs
    ):
        self.msas = msas

        # Store raw inputs for combinatorics config generation
        self.proteins = proteins
        self.ssDNA = ssDNA
        self.dsDNA = dsDNA
        self.ssRNA = ssRNA
        self.dsRNA = dsRNA
        self.ligands = ligands

        # Resolve inputs to DataStreams
        self.proteins_stream = resolve_input_to_datastream(proteins, fallback_stream="sequences")
        self.ssDNA_stream = resolve_input_to_datastream(ssDNA, fallback_stream="sequences")
        self.dsDNA_stream = resolve_input_to_datastream(dsDNA, fallback_stream="sequences")
        self.ssRNA_stream = resolve_input_to_datastream(ssRNA, fallback_stream="sequences")
        self.dsRNA_stream = resolve_input_to_datastream(dsRNA, fallback_stream="sequences")
        self.ligands_stream = None
        if ligands is not None:
            self.ligands_stream = resolve_input_to_datastream(ligands, fallback_stream="compounds")

        # BaseConfig reads kwargs with .get() and ignores the rest, so an
        # unsupported `pocket=` would be swallowed in silence — which is the
        # failure mode the parameter was removed to avoid.
        if "pocket" in kwargs:
            raise ValueError(
                "ESMFold2 has no pocket parameter. Upstream accepts "
                "PocketConditioning and then discards it: as of esm 3.3.0, "
                "build_feature_tensors zeroes pocket_feature unconditionally and "
                "never reads input.pocket, so the constraint would change nothing.")

        self.modifications = modifications
        self.covalent_bonds = covalent_bonds
        self.glycans = glycans

        self.num_loops = num_loops
        self.num_sampling_steps = num_sampling_steps
        self.num_diffusion_samples = num_diffusion_samples
        self.num_seeds = num_seeds
        self.msa_max_depth = msa_max_depth
        self.include_pae = include_pae
        self.top_only = top_only
        self.model_name = model_name

        super().__init__(**kwargs)

    # ---------------------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------------------

    def validate_params(self):
        has_input = any([
            self.proteins_stream is not None,
            self.ssDNA_stream is not None,
            self.dsDNA_stream is not None,
            self.ssRNA_stream is not None,
            self.dsRNA_stream is not None,
        ])
        if not has_input:
            raise ValueError("at least one sequence axis (proteins/ssDNA/dsDNA/ssRNA/dsRNA) is required")

        for name, val in (
            ("num_loops", self.num_loops),
            ("num_sampling_steps", self.num_sampling_steps),
            ("num_diffusion_samples", self.num_diffusion_samples),
            ("num_seeds", self.num_seeds),
            ("msa_max_depth", self.msa_max_depth),
        ):
            if not isinstance(val, int) or val < 1:
                raise ValueError(f"{name} must be a positive integer")

        if not isinstance(self.top_only, bool):
            raise ValueError("top_only must be a bool")
        if not isinstance(self.include_pae, bool):
            raise ValueError("include_pae must be a bool")

        if self.modifications is not None:
            _validate_modifications(self.modifications)
        if self.covalent_bonds is not None:
            _validate_covalent_bonds(self.covalent_bonds)
        if self.glycans is not None:
            _validate_glycans(self.glycans)

        _validate_freeform_string("model_name", self.model_name)

    # ---------------------------------------------------------------------------
    # Configure inputs
    # ---------------------------------------------------------------------------

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def _build_combinatorics_kwargs(self) -> Dict:
        kwargs = {}
        if self.proteins is not None:
            kwargs['proteins'] = (self.proteins, "sequences", "protein")
        if self.ssDNA is not None:
            kwargs['ssDNA'] = (self.ssDNA, "sequences", "ssdna")
        if self.dsDNA is not None:
            kwargs['dsDNA'] = (self.dsDNA, "sequences", "dsdna")
        if self.ssRNA is not None:
            kwargs['ssRNA'] = (self.ssRNA, "sequences", "ssrna")
        if self.dsRNA is not None:
            kwargs['dsRNA'] = (self.dsRNA, "sequences", "dsrna")
        if self.ligands is not None:
            kwargs['ligands'] = (self.ligands, "compounds", "ligand")
        return kwargs

    def _constraints_payload(self) -> Optional[Dict[str, Any]]:
        """Modification and covalent-bond constraints, or None if unused.

        Positions stay 1-indexed here and are converted to upstream's 0-indexed
        form in the inference script, so the emitted JSON reads the way the user
        wrote it.
        """
        payload: Dict[str, Any] = {}
        if self.modifications is not None:
            payload["modifications"] = self.modifications
        if self.covalent_bonds is not None:
            payload["covalent_bonds"] = self.covalent_bonds
        if self.glycans is not None:
            payload["glycans"] = self.glycans
        return payload or None

    def _missing_input_sources(self):
        return (
            self.proteins, self.proteins_stream,
            self.ssDNA, self.ssDNA_stream,
            self.dsDNA, self.dsDNA_stream,
            self.ssRNA, self.ssRNA_stream,
            self.dsRNA, self.dsRNA_stream,
            self.ligands, self.ligands_stream,
        )

    # ---------------------------------------------------------------------------
    # Config display
    # ---------------------------------------------------------------------------

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.extend([
            f"MODEL: {self.model_name}",
            f"LOOPS: {self.num_loops} | SAMPLING STEPS: {self.num_sampling_steps}",
            f"DIFFUSION SAMPLES: {self.num_diffusion_samples} | SEEDS: {self.num_seeds}",
            f"MSA RECYCLING: {'on' if self.msas is not None else 'off (single-sequence)'}",
            f"INCLUDE PAE: {self.include_pae}",
            f"TOP ONLY: {self.top_only}",
        ])
        if self.modifications is not None:
            mods = ", ".join(
                f"{m['chain']}{m['position']}->{m['ccd']}" for m in self.modifications)
            lines.append(f"MODIFICATIONS: {mods}")
        if self.covalent_bonds is not None:
            bonds = ", ".join(
                f"{b['atom1'][0]}{b['atom1'][1]}.{b['atom1'][2]}-"
                f"{b['atom2'][0]}{b['atom2'][1]}.{b['atom2'][2]}"
                for b in self.covalent_bonds)
            lines.append(f"COVALENT BONDS: {bonds}")
        return lines

    # ---------------------------------------------------------------------------
    # Script generation
    # ---------------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        generate_combinatorics_config(
            self.combinatorics_config_file, **self._build_combinatorics_kwargs()
        )

        constraints_flag = ""
        constraints = self._constraints_payload()
        if constraints is not None:
            os.makedirs(os.path.dirname(self.constraints_config_file), exist_ok=True)
            with open(self.constraints_config_file, "w", encoding="utf-8") as f:
                json.dump(constraints, f, indent=2)
            constraints_flag = f' --constraints "{self.constraints_config_file}"'

        msas_table = ""
        if self.msas is not None:
            msa_stream = getattr(self.msas.streams, "msas", None)
            if msa_stream is None or not msa_stream.map_table:
                raise ValueError("msas must be a tool output carrying an msas table")
            msas_table = msa_stream.map_table

        script = "#!/bin/bash\n"
        script += "# ESMFold2 all-atom complex prediction\n"
        script += self.generate_completion_check_header()

        # Phase 1: inference under the esmfold2 env
        script += self.activate_environment()  # esmfold2 (primary env)

        pae_flag = " --include-pae" if self.include_pae else ""
        top_flag = "" if self.top_only else " --all-samples"
        msa_flag = f' --msas-table "{msas_table}"' if msas_table else ""

        script += f'echo "Running ESMFold2 ({self.model_name})"\n\n'
        script += f"""python "{self.inference_py}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --output-dir "{self.predictions_folder}" \\
    --model-name "{self.model_name}" \\
    --num-loops {self.num_loops} \\
    --num-sampling-steps {self.num_sampling_steps} \\
    --num-diffusion-samples {self.num_diffusion_samples} \\
    --num-seeds {self.num_seeds} \\
    --msa-max-depth {self.msa_max_depth}{pae_flag}{top_flag}{msa_flag}{constraints_flag}

if [ $? -ne 0 ]; then
    echo "Error: ESMFold2 inference failed"
    exit 1
fi

"""

        # Phase 2: post-processing under the biopipelines env
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")

        compounds_flag = ""
        if self.ligands_stream is not None:
            compounds_flag = (
                f' --ligands-table "{self.ligands_stream.map_table}"'
                f' --compounds-map "{self.compounds_map_csv}"'
            )

        script += f"""python "{self.postprocess_py}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --predictions-dir "{self.predictions_folder}" \\
    --structures-map "{self.structures_map_csv}" \\
    --confidence-csv "{self.confidence_csv}"{compounds_flag}

if [ $? -eq 0 ]; then
    echo "ESMFold2 post-processing completed successfully"
else
    echo "Error: ESMFold2 post-processing failed"
    exit 1
fi

"""
        # Condition mirrors get_output_files() so declared <=> produced.
        if self._collect_upstream_missing_paths(*self._missing_input_sources()):
            script += self._generate_missing_table_propagation()
        script += self.generate_completion_check_footer()
        return script

    # ---------------------------------------------------------------------------
    # Output files
    # ---------------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        predicted_ids, _provenance = predict_output_ids_with_provenance(
            **self._build_combinatorics_kwargs()
        )

        structure_files = [self.stream_path("structures", "<id>.cif")]
        structure_ids = predicted_ids
        if not self.top_only:
            k = self.num_diffusion_samples
            structure_ids = [f"{pid}_<1..{k}>" for pid in predicted_ids]

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map_csv,
            format="cif",
        )

        confidence_columns = ["id", "file", "plddt", "ptm", "iptm"]
        if self.include_pae:
            confidence_columns.append("max_pae")

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map_csv,
                columns=["id", "file"],
                description="ESMFold2 predicted structures",
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=confidence_columns,
                description="ESMFold2 confidence scores (pLDDT, pTM, ipTM)",
            ),
        }

        compounds = None
        if self.ligands_stream is not None:
            compounds = DataStream(
                name="compounds",
                ids=self.ligands_stream.ids,
                files=[],
                map_table=self.compounds_map_csv,
                format="csv",
            )
            tables["compounds"] = TableInfo(
                name="compounds",
                path=self.compounds_map_csv,
                columns=["id", "format", "code", "smiles", "ccd"],
                description="Ligand compounds with the residue code ESMFold2 used",
            )

        if self._collect_upstream_missing_paths(*self._missing_input_sources()):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": compounds if compounds is not None else DataStream.empty("compounds", "csv"),
            "tables": tables,
            "output_folder": self.output_folder,
            "rendering_parameters": {
                "structures": {
                    "color_by": "plddt",
                    "plddt_upper": 100,
                }
            }
        }

    def _generate_missing_table_propagation(self) -> str:
        return self.generate_missing_propagation(
            *self._missing_input_sources(), missing_csv=self.missing_csv
        )
