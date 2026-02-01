# Tool Refactoring Guide

## Refactoring Checklist

### Structure Generation
- [x] RFdiffusion
- [x] RFdiffusionAllAtom
- [x] RFdiffusion3
- [x] BoltzGen (+ BoltzGenMerge, BoltzGenImport)

### Sequence Design
- [x] ProteinMPNN
- [x] LigandMPNN
- [x] MutationComposer
- [x] SDM (SiteDirectedMutagenesis)
- [x] Fuse
- [x] StitchSequences
- [x] SplitChains
- [x] DNAEncoder

### Structure Prediction
- [x] AlphaFold
- [x] GhostFold (removed)
- [x] ESMFold
- [x] Boltz2
- [x] RF3 (removed)
- [x] OnionNet

### Analysis
- [x] DynamicBind
- [x] ResidueAtomDistance
- [x] PLIP
- [x] DistanceSelector
- [x] ConformationalChange
- [x] MutationProfiler
- [x] SequenceMetricCorrelation
- [x] BayesianAdjuster
- [x] ProteinLigandContacts
- [x] PoseDistance
- [x] SASA
- [x] SelectionEditor

### Data Management
- [x] Filter
- [x] Rank
- [x] SelectBest
- [x] RemoveDuplicates
- [X] MergeTables
- [x] ConcatenateTables
- [x] SliceTable
- [x] ExtractMetrics
- [x] AverageByTable

### Utilities
- [x] LoadOutput / LoadOutputs (legacy support retained)
- [x] MMseqs2
- [x] MMseqs2Server
- [x] MMseqs2LCF (+ MMseqs2ServerLCF)
- [x] CompoundLibrary
- [x] PDB
- [x] Ligand
- [x] PyMOL
- [x] Plot

---

## New Classes

### `DataStream` (datastream.py)
Unified container for tool I/O. Replaces separate `structures`/`structure_ids`, `sequences`/`sequence_ids`, etc.

```python
DataStream(name="structures", ids=["a", "b"], files=["/path/a.pdb", "/path/b.pdb"], map_table="/path/map.csv", format="pdb")
DataStream.empty("compounds", "sdf")  # Empty stream
```

### `Path` (file_paths.py)
Lazy path descriptor. Replaces `_initialize_file_paths()` and `_setup_file_paths()`.

```python
class MyTool(BaseConfig):
    main_table = Path(lambda self: os.path.join(self.output_folder, "results.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_script.py"))
```

## BaseConfig Provides

- `self.pipeline_name` - extracted from output folder using `TOOL_NAME`
- `self.log_file` - computed from folder naming pattern

## Input Rules

2. **No fallbacks** - code crashes if input is wrong
3. **No legacy support** - remove all deprecated parameters
4. Structure and sequence inputs are `DataStream` or `StandardizedOutput` only

## Tool Migration Steps

1. **Remove** `_initialize_file_paths()`, `_setup_file_paths()`, `_extract_pipeline_name()`
2. **Remove** all calls to these methods
3. **Remove** `ToolOutput` handling - only accept `DataStream` or `StandardizedOutput`
4. **Remove** string/list/folder input handling for pdb/structure arguments, except the PDB tool
5. **Remove** all fallback logic and legacy parameters, including 'main' table
6. **Add** `Path` descriptors at class level for all file paths
7. **Simplify** `configure_inputs` to just `self.folders = pipeline_folders`
8. **Update** `get_output_files()` to return `DataStream` objects

## Example: Refactored Tool

```python
class MyTool(BaseConfig):
    TOOL_NAME = "MyTool"

    # Path descriptors
    main_table = Path(lambda self: os.path.join(self.output_folder, "results.csv"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mytool.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_outputs: int = 1,
                 **kwargs):
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream = structures.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput")

        self.num_outputs = num_outputs
        super().__init__(**kwargs)

    def validate_params(self):
        if len(self.structures_stream) == 0:
            raise ValueError("structures required")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_output_files(self) -> Dict[str, Any]:
        structures = DataStream(
            name="structures",
            ids=[...],
            files=[...],
            map_table=self.main_table,
            format="pdb"
        )
        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": {...},
            "output_folder": self.output_folder
        }
```
