# Developer Guide - Pipeline Architecture

Technical documentation for extending and customizing the ProteinNotebooks pipeline system.

## Architecture Overview

The pipeline system is built around three core modules:

```
ConfigScripts/    # Tool configuration classes
├── base_config.py        # Abstract base class
├── rfdiffusion.py       # RFdiffusion implementation  
├── protein_mpnn.py      # ProteinMPNN implementation
├── alphafold.py         # AlphaFold implementation
└── [future tools]       # Boltz2, LigandMPNN, etc.

Utilities/        # Core pipeline management
├── pipeline.py          # Main Pipeline orchestrator
└── folders.py           # Folder structure management

IOScripts/        # Format conversion utilities
└── converters.py        # PDB, FASTA, JSON conversions
```

## Core Design Principles

### 1. Rich Return Objects

When tools are added to pipelines, they return `ToolOutput` objects with rich metadata:

```python
rfd = pipeline.add(RFdiffusion.RFdiffusion(...))

# Rich object with multiple access patterns
rfd.output_pdbs              # List of PDB file paths
rfd.get_output_files('pdbs') # Same as above  
rfd.tool_type               # "RFdiffusion"
rfd.environment             # "ProteinEnv"
rfd.dependencies            # List of upstream tools
```

### 2. Automatic I/O Resolution

Tools automatically find their inputs from upstream tools:

```python
# ProteinMPNN automatically finds RFdiffusion PDB outputs
mpnn = pipeline.add(ProteinMPNN.ProteinMPNN(
    input_structures=rfd  # Resolves to rfd.output_pdbs
))
```

### 3. Environment Management

Pipeline automatically groups tools by environment and minimizes switches:

```python
# Generated script will be:
# source activate ProteinEnv
#   - Run RFdiffusion  
#   - Run ProteinMPNN
# source activate Boltz2Env
#   - Run Boltz2
```

## Implementing New Tools

### Step 1: Create Tool Configuration Class

Extend `BaseConfig` for each new tool:

```python
from .base_config import BaseConfig

class NewTool(BaseConfig):
    TOOL_NAME = "NewTool"
    DEFAULT_ENV = "ProteinEnv" 
    COMPATIBLE_ENVS = ["ProteinEnv", "OtherEnv"]
    
    def __init__(self, param1: str, param2: int = 10, **kwargs):
        self.param1 = param1
        self.param2 = param2
        super().__init__(**kwargs)
    
    def validate_params(self):
        if not self.param1:
            raise ValueError("param1 is required")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        # Set up input sources from pipeline context
        pass
    
    def generate_script(self, script_path: str) -> str:
        # Generate bash execution script
        return "#!/bin/bash\necho 'Running NewTool'"
    
    def get_output_files(self) -> Dict[str, List[str]]:
        # Return expected output files
        return {"results": ["/path/to/output.txt"]}
```

### Step 2: Handle Input Resolution

Tools must handle various input types:

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    if isinstance(self.input_data, ToolOutput):
        # Input from previous tool
        source_files = self.input_data.get_output_files("pdbs")
        self.dependencies.append(self.input_data.config)
        
    elif isinstance(self.input_data, str):
        # File or folder path
        if self.input_data.endswith('.pdb'):
            source_files = [self.input_data]
        else:
            # Folder of files
            folder_files = glob.glob(f"{self.input_data}/*.pdb")
    
    # Copy to runtime folder for processing
    runtime_folder = os.path.join(self.output_folder, "RunTime")
    for file in source_files:
        runtime_path = os.path.join(runtime_folder, os.path.basename(file))
        self.input_sources[os.path.basename(file)] = file
```

### Step 3: Script Generation

Generate bash scripts that integrate with existing helper scripts:

```python
def generate_script(self, script_path: str) -> str:
    runtime_folder = os.path.dirname(script_path)
    
    # Copy input files
    for filename, source in self.input_sources.items():
        shutil.copy(source, os.path.join(runtime_folder, filename))
    
    # Tool-specific paths
    tool_executable = os.path.join(self.folders["NewTool"], "bin", "newtool")
    output_folder = os.path.join(self.output_folder, "NewToolOutput")
    
    return f"""#!/bin/bash
echo "Starting NewTool: {self.job_name}"
mkdir -p {output_folder}

# Run tool
{tool_executable} --input {runtime_folder} --output {output_folder}

echo "NewTool completed"
"""
```

### Step 4: Output Specification

Define expected outputs for downstream tools:

```python
def get_output_files(self) -> Dict[str, List[str]]:
    output_folder = os.path.join(self.output_folder, "NewToolOutput")
    
    # Generate expected file paths
    result_files = []
    for i in range(self.num_outputs):
        result_files.append(f"{output_folder}/result_{i}.txt")
    
    return {
        "results": result_files,
        "summary": [f"{output_folder}/summary.json"],
        "output_folder": [output_folder]
    }
```

## Extending Pipeline Functionality

### Custom Analysis Tools

Create analysis components that can be added to pipelines:

```python
class CustomAnalysis(BaseConfig):
    TOOL_NAME = "CustomAnalysis"
    
    def __init__(self, structures: List[ToolOutput], metric: str, **kwargs):
        self.structures = structures
        self.metric = metric
        super().__init__(**kwargs)
    
    def configure_inputs(self, pipeline_folders):
        # Add all upstream tools as dependencies
        for struct in self.structures:
            self.dependencies.append(struct.config)
```

### Conditional Execution

Implement conditional pipeline steps:

```python
# In pipeline construction
rfd = pipeline.add(RFdiffusion.RFdiffusion(...))

# Only run expensive analysis if many designs were generated
if len(rfd.output_pdbs) > 10:
    analysis = pipeline.add(ExtensiveAnalysis.ExtensiveAnalysis(
        structures=rfd
    ))
else:
    analysis = pipeline.add(BasicAnalysis.BasicAnalysis(
        structures=rfd  
    ))
```

### Branching Workflows

Create parallel processing branches:

```python
# Single input, multiple processing paths
mpnn = pipeline.add(ProteinMPNN.ProteinMPNN(input_structures=rfd))

# Parallel folding with different methods
af = pipeline.add(AlphaFold.AlphaFold(sequences=mpnn))
omega = pipeline.add(OmegaFold.OmegaFold(sequences=mpnn))

# Compare results
comparison = pipeline.add(Analysis.Compare(
    structures=[af, omega],
    reference=rfd
))
```

## Advanced Features

### Resource Management

Tools can specify detailed resource requirements:

```python
class GPUIntensiveTool(BaseConfig):
    DEFAULT_RESOURCES = {
        "gpu": "A100",           # Specific GPU requirement
        "memory": "32GB",        # Memory requirement  
        "time": "2-00:00:00",   # 2 days
        "nodes": 1,             # Multi-node if needed
        "tasks_per_node": 1
    }
```

### Batch Processing

Implement tools that process multiple inputs efficiently:

```python
class BatchProcessor(BaseConfig):
    def generate_script(self, script_path: str) -> str:
        # Generate script that processes all inputs in parallel
        script = "#!/bin/bash\n"
        
        for i, input_file in enumerate(self.input_files):
            script += f"process_input {input_file} &\n"
        
        script += "wait  # Wait for all background jobs\n"
        return script
```

### Environment Validation

Add environment compatibility checking:

```python
class SpecialTool(BaseConfig):
    @classmethod
    def validate_environment_setup(cls, env_name: str) -> bool:
        """Check if environment has required dependencies."""
        try:
            # Check for required executables, packages, etc.
            subprocess.run(["which", "special_executable"], 
                         check=True, capture_output=True)
            return True
        except subprocess.CalledProcessError:
            return False
```

## Testing and Validation

### Unit Testing Tools

Test individual tool configurations:

```python
import unittest
from ConfigScripts import NewTool
from Utilities import Pipeline

class TestNewTool(unittest.TestCase):
    def setUp(self):
        self.pipeline = Pipeline.Pipeline("test")
    
    def test_parameter_validation(self):
        with self.assertRaises(ValueError):
            NewTool.NewTool(param1="")  # Should fail
    
    def test_output_generation(self):
        tool = NewTool.NewTool(param1="test")
        outputs = tool.get_output_files()
        self.assertIn("results", outputs)
```

### Integration Testing

Test complete workflows:

```python
def test_complete_workflow():
    pipeline = Pipeline.Pipeline("integration_test")
    
    # Build pipeline
    tool1 = pipeline.add(Tool1.Tool1(...))
    tool2 = pipeline.add(Tool2.Tool2(input=tool1))
    
    # Validate pipeline
    assert pipeline.validate_pipeline()
    
    # Generate scripts
    script_path = pipeline.save()
    assert os.path.exists(script_path)
```

## Performance Optimization

### Efficient Resource Usage

- **Group by environment**: Minimize environment switches
- **Parallel execution**: Use job arrays for independent tasks  
- **Memory management**: Process large datasets in chunks
- **Caching**: Store expensive computations (MSAs, models)

### Script Optimization

```python
def generate_optimized_script(self, script_path: str) -> str:
    # Use efficient file operations
    script = "#!/bin/bash\nset -e\n"
    
    # Parallel processing where possible
    if len(self.input_files) > 4:
        script += "# Process inputs in parallel\n"
        for file in self.input_files:
            script += f"process_file {file} &\n"
        script += "wait\n"
    
    # Cleanup intermediate files
    script += "# Cleanup\nrm -rf temp_files/\n"
    
    return script
```

## Debugging and Troubleshooting

### Logging Integration

Tools should provide detailed logging:

```python
def generate_script(self, script_path: str) -> str:
    log_file = os.path.join(self.output_folder, f"_{self.TOOL_NAME.lower()}.log")
    
    return f"""#!/bin/bash
exec > >(tee {log_file}) 2>&1  # Log everything

echo "Starting {self.TOOL_NAME} at $(date)"
echo "Parameters: {json.dumps(self.to_dict(), indent=2)}"

# Tool execution
run_tool_command

echo "Completed at $(date)"
"""
```

### Error Handling

```python
def validate_outputs(self):
    """Validate that expected outputs were generated."""
    expected = self.get_output_files()
    
    for output_type, files in expected.items():
        for file_path in files:
            if not os.path.exists(file_path):
                raise RuntimeError(f"Expected output not found: {file_path}")
            
            if os.path.getsize(file_path) == 0:
                raise RuntimeError(f"Output file is empty: {file_path}")
```

## Migration and Compatibility

### Legacy Script Integration

Integrate existing helper scripts:

```python
def generate_script(self, script_path: str) -> str:
    # Use existing helper scripts
    helper_script = os.path.join(self.folders["HelpScripts"], "legacy_tool.py")
    
    return f"""#!/bin/bash
# Use legacy helper script
python {helper_script} {self.output_folder} {' '.join(self.parameters)}
"""
```

### Backward Compatibility

Maintain compatibility with existing workflows while adding new features:

```python
class ModernTool(BaseConfig):
    def __init__(self, legacy_mode: bool = False, **kwargs):
        self.legacy_mode = legacy_mode
        super().__init__(**kwargs)
    
    def generate_script(self, script_path: str) -> str:
        if self.legacy_mode:
            return self._generate_legacy_script(script_path)
        else:
            return self._generate_modern_script(script_path)
```

This architecture provides a robust foundation for extending the pipeline system while maintaining simplicity for end users.