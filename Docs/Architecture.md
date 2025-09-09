# Pipeline System Architecture Overview

## System Design Philosophy

The ProteinNotebooks pipeline system is designed around several key principles:

1. **Clean Separation of Concerns**: Each tool manages only its own outputs
2. **No File Copying**: Tools work with files in their original locations
3. **Dependency Management**: Automatic resolution of tool dependencies
4. **Environment Isolation**: Proper conda environment switching
5. **Consistent Interface**: Unified API across all tools

## Core Components

### 1. Pipeline Orchestrator (`Pipeline` class)

The central coordinator that:
- Manages tool execution order
- Handles environment switching
- Generates unified execution scripts
- Provides SLURM integration
- Tracks dependencies

```python
pipeline = Pipeline("MyPipeline", "JobName", "Description")
tool1 = pipeline.add(ToolA(...))
tool2 = pipeline.add(ToolB(tool1))  # tool2 depends on tool1
pipeline.save()  # Generate execution script
```

### 2. Base Configuration (`BaseConfig` class)

Abstract base class providing:
- **Parameter validation**: `validate_params()`
- **Input configuration**: `configure_inputs()`
- **Script generation**: `generate_script()`
- **Output definition**: `get_output_files()`
- **Resource management**: SLURM resource specification
- **Environment compatibility**: Conda environment validation

### 3. Tool Output Wrapper (`ToolOutput` class)

Encapsulates tool outputs with:
- **Type-safe access**: `output.output_pdbs`, `output.output_sequences`
- **Flexible querying**: `output.get_output_files("structures")`
- **Metadata tracking**: execution order, environment, folder paths
- **Dependency linking**: automatic dependency resolution

### 4. Folder Management (`FolderManager` class)

Handles directory structure:
```
Pipeline_001/
├── RunTime/              # Execution scripts only
│   ├── pipeline.sh
│   ├── slurm.sh
│   └── tool_scripts/
├── 1_ToolA/             # ToolA outputs only  
│   ├── results.txt
│   └── data/
└── 2_ToolB/             # ToolB outputs only
    ├── final.pdb
    └── analysis.csv
```

## Data Flow Architecture

### Input/Output Chain

```
PDB File → RFdiffusion → ProteinMPNN → AlphaFold → Results
   ↓           ↓            ↓           ↓
 PDBs/    1_RFdiffusion  2_ProteinMPNN  3_AlphaFold
          structures.pdb   sequences.fa   ranked_0.pdb
```

### Dependency Resolution

Tools can accept inputs from:
1. **Previous Tools**: `ToolB(input_sequences=tool_a_output)`
2. **File Paths**: `ToolA(input_structure="protein.pdb")`
3. **Direct Data**: `ToolC(ligand_smiles="CCO")`

The pipeline automatically:
- Resolves file paths from previous tools
- Maintains execution order based on dependencies  
- Passes correct folder references
- Handles environment switching

### Environment Management

```python
# Automatic environment switching in generated script
source activate ProteinEnv
./1_rfdiffusion.sh | tee _1_rfdiffusion.log

source activate ligandmpnn_env  
./2_ligandmpnn.sh | tee _2_ligandmpnn.log

source activate Boltz2Env
./3_boltz2.sh | tee _3_boltz2.log
```

## Tool Implementation Pattern

### 1. Class Definition

```python
class NewTool(BaseConfig):
    TOOL_NAME = "NewTool"
    DEFAULT_ENV = "tool_env"
    COMPATIBLE_ENVS = ["tool_env", "alt_env"]
    DEFAULT_RESOURCES = {"gpu": "V100", "memory": "16GB", "time": "24:00:00"}
```

### 2. Initialization

```python
def __init__(self, input_param: Union[str, ToolOutput], **kwargs):
    self.input_param = input_param
    self.input_is_tool_output = isinstance(input_param, ToolOutput)
    super().__init__(**kwargs)
```

### 3. Input Configuration

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    if self.input_is_tool_output:
        # Get files from previous tool
        files = self.input_param.get_output_files("relevant_type")
        self.input_file = files[0]
        self.dependencies.append(self.input_param.config)
    else:
        # Handle direct file input
        self.input_file = os.path.join(pipeline_folders["PDBs"], self.input_param)
```

### 4. Script Generation

```python
def generate_script(self, script_path: str) -> str:
    # Use input_file where it is (no copying!)
    # Output to self.output_folder only
    return f"""#!/bin/bash
tool_command --input "{self.input_file}" --output "{self.output_folder}"
"""
```

### 5. Output Definition

```python
def get_output_files(self) -> Dict[str, List[str]]:
    return {
        "structures": [os.path.join(self.output_folder, "result.pdb")],
        "datasheets": [os.path.join(self.output_folder, "summary.csv")]
    }
```

## Key Design Patterns

### 1. Clean Output Organization

Each tool outputs **only** to its own numbered folder, example:
- `1_RFdiffusion/` - only RFdiffusion outputs
- `2_ProteinMPNN/` - only ProteinMPNN outputs  
- `3_AlphaFold/` - only AlphaFold outputs
- `RunTime/` - only execution scripts

### 2. Flexible Input Handling

Tools accept multiple input types transparently:
```python
# From previous tool
tool_b = pipeline.add(ToolB(tool_a_output))

# From file
tool_a = pipeline.add(ToolA("input.pdb"))

# From data
tool_c = pipeline.add(ToolC(["SMILES1", "SMILES2"]))
```

### 4. Consistent Output Types

All tools can provide standard output types:
- `structures`: PDB/mmCIF files
- `sequences`: csv file containing columns 'id', 'sequence'
- `datasheets`: CSV/JSON analysis files, in form of dictionary. the one called 'main' contains a mapping of ids to output
- `output_folder`: Tool's output directory

## Script Generation Flow

### 1. Pipeline Assembly

```python
pipeline = Pipeline("MyPipeline", "JobName", "Description")
tool1 = pipeline.add(ToolA(...))  # Creates 1_ToolA/ folder
tool2 = pipeline.add(ToolB(...))  # Creates 2_ToolB/ folder
tool3 = pipeline.add(ToolC(...))  # Creates 3_ToolC/ folder
```

### 2. Dependency Resolution

```python
pipeline.save()  # Triggers:
# 1. Validate all tools
# 2. Configure inputs based on dependencies
# 3. Generate individual tool scripts
# 4. Create unified pipeline script
# 5. Set up environment switching
```

### 3. Generated Script Structure

```bash
#!/bin/bash
source activate ProteinEnv
echo Configuration
./config.sh | tee pipeline_config.txt

echo ToolA
./1_toola.sh | tee _1_toola.log

source activate different_env
echo ToolB  
./2_toolb.sh | tee _2_toolb.log

echo Job done
```

## Performance Considerations

### Resource Management

- GPU allocation per tool
- Memory requirements tracked
- Time limits enforced
- Environment switching minimized

### Parallelization

While tools run sequentially within a pipeline, multiple pipelines can run in parallel:
```python
# Multiple independent jobs
for protein in protein_list:
    pipeline = Pipeline(f"Design_{protein}", protein, "Parallel design")
    # ... configure pipeline ...
    pipeline.save()  # Each gets unique folders
```

### Storage Efficiency

- No file duplication between tools
- Intermediate files in tool-specific folders
- Automatic cleanup opportunities
- Compressed result archives
