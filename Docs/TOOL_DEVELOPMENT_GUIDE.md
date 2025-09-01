# Tool Development Guide

## Overview

This guide explains how to create and modify tool classes in the PipelineScripts system. The system provides a unified interface for protein modeling tools with automatic dependency management, environment switching, and clean output organization. The tools write scripts that will give some output and provide the predicted structure of the output, however the actual output files will only be there at runtime of the generated scripts.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Base Classes](#base-classes)
3. [Creating New Tools](#creating-new-tools)
4. [Input/Output Management](#inputoutput-management)
5. [Script Generation](#script-generation)
6. [Best Practices](#best-practices)
7. [Examples](#examples)

## Architecture Overview

The pipeline system follows a clean separation of concerns:

```
Pipeline
├── Tool 1 (BaseConfig)
│   ├── Input Configuration
│   ├── Script Generation
│   └── Output Definition
├── Tool 2 (BaseConfig)
│   └── ...
└── Tool N (BaseConfig)
```

### Key Principles

1. **No File Copying**: Tools work with files in their original locations.
2. **Clean Separation**: Each tool outputs only to its own folder
3. **Dependency Chain**: Tools can use outputs from previous tools as inputs
4. **Environment Management**: Automatic conda environment switching
5. **Consistent Interface**: All tools follow the same patterns

## Base Classes

### BaseConfig

The `BaseConfig` class is the foundation for all tools. It provides:

#### Core Attributes

```python
class BaseConfig(ABC):
    TOOL_NAME = "base"                    # Override in subclasses
    DEFAULT_ENV = "ProteinEnv"            # Default conda environment
    COMPATIBLE_ENVS = ["ProteinEnv"]      # List of compatible environments
    DEFAULT_RESOURCES = {                 # Default SLURM resources
        "gpu": "V100", 
        "memory": "15GB", 
        "time": "24:00:00"
    }
```

#### Instance Variables

- `self.tool_name`: Tool identifier
- `self.job_name`: User-provided job name
- `self.environment`: Conda environment to use
- `self.resources`: SLURM resource requirements
- `self.dependencies`: List of dependent tools
- `self.pipeline_ref`: Reference to parent pipeline
- `self.execution_order`: Order in pipeline (1, 2, 3...)
- `self.input_sources`: Dictionary of input file mappings
- `self.output_files`: Dictionary of expected outputs
- `self.output_folder`: Tool's dedicated output folder
- `self.folders`: Dictionary of all pipeline folders

#### Abstract Methods (Must Implement)

```python
@abstractmethod
def validate_params(self):
    """Validate tool-specific parameters."""
    pass

@abstractmethod
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    """Configure input sources from pipeline context."""
    pass

@abstractmethod
def generate_script(self, script_path: str) -> str:
    """Generate bash script for tool execution."""
    pass

@abstractmethod
def get_output_files(self) -> Dict[str, List[str]]:
    """Get expected output files after execution."""
    pass
```

#### Optional Override Methods

```python
def get_config_display(self) -> List[str]:
    """Configuration display lines for pipeline output."""
    return ["Job: " + self.job_name] if self.job_name else []

def get_expected_output_paths(self) -> Dict[str, List[str]]:
    """Get expected output paths without validation."""
    return self.get_output_files()

def to_dict(self) -> Dict[str, Any]:
    """Serialize configuration to dictionary."""
    return {
        'tool_name': self.tool_name,
        'environment': self.environment,
        'resources': self.resources,
        'parameters': self.params,
        'execution_order': self.execution_order
    }
```

### ToolOutput

The `ToolOutput` class wraps tool outputs and provides convenient access methods:

#### Key Attributes

- `self.config`: Reference to the tool configuration
- `self.tool_type`: Tool name (e.g., "ProteinMPNN")
- `self.environment`: Tool's environment
- `self.output_folder`: Tool's output folder
- `self.execution_order`: Tool's execution order
- `self._output_files`: Dictionary of actual output files

#### Convenience Properties

```python
@property
def output_pdbs(self) -> List[str]:
    """Get PDB files."""
    return self.get_output_files('pdbs')

@property
def output_sequences(self) -> List[str]:
    """Get sequence files."""
    return self.get_output_files('sequences')

@property
def output_structures(self) -> List[str]:
    """Get structure files (alias for PDbs)."""
    return self.output_pdbs

@property
def output_datasheets(self) -> List[str]:
    """Get datasheet files (CSV, JSON, etc.)."""
    return self.get_output_files('datasheets')
```

## Creating New Tools

### Step 1: Define Tool Class

```python
from typing import Dict, List, Any, Optional, Union
from .base_config import BaseConfig, ToolOutput

class MyTool(BaseConfig):
    # Tool identification
    TOOL_NAME = "MyTool"
    DEFAULT_ENV = "my_env"
    COMPATIBLE_ENVS = ["my_env", "alternative_env"]
    DEFAULT_RESOURCES = {"gpu": "V100", "memory": "16GB", "time": "12:00:00"}
    
    def __init__(self, 
                 input_param: Union[str, ToolOutput],
                 specific_param: int = 10,
                 **kwargs):
        """Initialize tool configuration."""
        # Store tool-specific parameters
        self.input_param = input_param
        self.specific_param = specific_param
        
        # Track input types
        self.input_is_tool_output = isinstance(input_param, ToolOutput)
        
        # Initialize base class
        super().__init__(**kwargs)
```

### Step 2: Implement Required Methods

```python
def validate_params(self):
    """Validate tool-specific parameters."""
    if not self.input_param:
        raise ValueError("input_param is required")
    
    if self.specific_param <= 0:
        raise ValueError("specific_param must be positive")

def configure_inputs(self, pipeline_folders: Dict[str, str]):
    """Configure input sources."""
    self.folders = pipeline_folders
    
    if self.input_is_tool_output:
        # Input from previous tool
        tool_output: ToolOutput = self.input_param
        source_files = tool_output.get_output_files("relevant_type")
        
        if not source_files:
            raise ValueError(f"No relevant files from {tool_output.tool_type}")
        
        # Store input file paths (don't copy files!)
        self.input_file_path = source_files[0]
        
        # Add dependency
        self.dependencies.append(tool_output.config)
        
    elif isinstance(self.input_param, str):
        # String input - file path
        if self.input_param.endswith('.ext'):
            input_path = os.path.join(pipeline_folders["PDBs"], self.input_param)
            if os.path.exists(input_path):
                self.input_file_path = input_path
            else:
                raise ValueError(f"File not found: {input_path}")
        else:
            raise ValueError("Invalid input format")
    
    else:
        raise ValueError(f"Invalid input type: {type(self.input_param)}")

def generate_script(self, script_path: str) -> str:
    """Generate bash script for execution."""
    # Tool executable paths
    tool_exe = os.path.join(self.folders["data"], "my_tool", "run_tool.py")
    
    # Build command-line options
    options = f"--input '{self.input_file_path}'"
    options += f" --output '{self.output_folder}'"
    options += f" --param {self.specific_param}"
    
    # Generate script content
    script_content = f"""#!/bin/bash
# MyTool execution script
# Generated by ProteinNotebooks pipeline system

echo "Running MyTool"
python {tool_exe} {options}

echo "MyTool completed"
"""
    return script_content

def get_output_files(self) -> Dict[str, List[str]]:
    """Get expected output files."""
    # Expected outputs in tool's output folder
    main_output = os.path.join(self.output_folder, f"{self.job_name}_result.txt")
    summary_csv = os.path.join(self.output_folder, f"{self.job_name}_summary.csv")
    
    return {
        "results": [main_output],
        "datasheets": [summary_csv],
        "output_folder": [self.output_folder]
    }
```

## Input/Output Management

### Input Handling Patterns

#### 1. Tool Output Inputs (Most Common)

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    if self.input_is_tool_output:
        tool_output: ToolOutput = self.input_param
        
        # Get specific output type from previous tool
        structures = tool_output.get_output_files("structures")
        sequences = tool_output.get_output_files("sequences")
        
        # Use the files where they are (no copying!)
        if structures:
            self.input_structure_path = structures[0]
        
        # Add dependency for proper execution order
        self.dependencies.append(tool_output.config)
```

#### 2. File Path Inputs

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    elif isinstance(self.input_param, str):
        # Handle different file types
        if self.input_param.endswith('.pdb'):
            file_path = os.path.join(pipeline_folders["PDBs"], self.input_param)
        elif self.input_param.endswith('.fasta'):
            file_path = os.path.join(pipeline_folders["data"], self.input_param)
        else:
            raise ValueError(f"Unsupported file type: {self.input_param}")
        
        if not os.path.exists(file_path):
            raise ValueError(f"File not found: {file_path}")
        
        self.input_file_path = file_path
```

#### 3. Direct Data Inputs

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    elif isinstance(self.input_param, list):
        # Handle list of SMILES, sequences, etc.
        self.input_data = self.input_param
    elif isinstance(self.input_param, dict):
        # Handle structured data
        self.input_config = self.input_param
```

### Output Definition Patterns

#### 1. Standard Output Types

Always provide these standard output types where applicable:

```python
def get_output_files(self) -> Dict[str, List[str]]:
    outputs = {}
    
    # Structures (PDB, mmCIF files)
    if self.produces_structures:
        structure_files = [os.path.join(self.output_folder, "structures", "*.pdb")]
        outputs["structures"] = structure_files
        outputs["pdbs"] = structure_files  # Alias
    
    # Sequences (FASTA files)
    if self.produces_sequences:
        seq_files = [os.path.join(self.output_folder, "sequences.fasta")]
        outputs["sequences"] = seq_files
        outputs["fasta_files"] = seq_files  # Alias
    
    # Datasheets (CSV, JSON analysis files)
    datasheets = [
        os.path.join(self.output_folder, f"{self.job_name}_results.csv"),
        os.path.join(self.output_folder, f"{self.job_name}_summary.json")
    ]
    outputs["datasheets"] = datasheets
    
    # Always include output folder for tools that need it
    outputs["output_folder"] = [self.output_folder]
    
    return outputs
```

#### 2. Tool-Specific Outputs

```python
def get_output_files(self) -> Dict[str, List[str]]:
    outputs = {
        # Standard outputs
        "structures": [...],
        "datasheets": [...],
        
        # Tool-specific outputs
        "alignments": [os.path.join(self.output_folder, "alignments.a3m")],
        "energies": [os.path.join(self.output_folder, "energies.txt")],
        "pymol_session": [os.path.join(self.output_folder, "visualization.pse")],
        
        # Intermediate files (if other tools need them)
        "config_files": [os.path.join(self.output_folder, "configs")],
        "temp_folder": [os.path.join(self.output_folder, "temp")]
    }
    return outputs
```

### Folder Organization

Each tool gets its own numbered output folder:

```
Pipeline_001/
├── RunTime/                    # Execution scripts only
│   ├── pipeline.sh
│   ├── slurm.sh
│   ├── 1_mytool.sh
│   └── 2_othertool.sh
├── 1_MyTool/                   # Tool outputs only
│   ├── results.txt
│   ├── summary.csv
│   └── subfolder/
└── 2_OtherTool/               # Next tool outputs
    ├── final_results.pdb
    └── analysis.json
```

## Script Generation

### Template Structure

```python
def generate_script(self, script_path: str) -> str:
    """Generate bash script following standard pattern."""
    
    # 1. Get tool installation paths
    tool_folder = os.path.join(self.folders["data"], "tool_name")
    tool_exe = os.path.join(tool_folder, "tool_executable")
    
    # 2. Build command-line options
    options = self._build_command_options()
    
    # 3. Helper script paths
    helper_script = os.path.join(self.folders["HelpScripts"], "helper.py")
    
    # 4. Generate script content with standard pattern
    script_content = f"""#!/bin/bash
# {self.TOOL_NAME} execution script
# Generated by ProteinNotebooks pipeline system

echo "Running {self.TOOL_NAME}"
cd {tool_folder}
{tool_exe} {options}

echo "Processing results"
python {helper_script} {self.output_folder} {self.job_name}

echo "{self.TOOL_NAME} completed"
"""
    
    return script_content

def _build_command_options(self) -> str:
    """Build command-line options string."""
    options = []
    
    # Input files
    options.append(f"--input '{self.input_file_path}'")
    
    # Output directory
    options.append(f"--output '{self.output_folder}'")
    
    # Tool-specific parameters
    if hasattr(self, 'num_sequences'):
        options.append(f"--num-seq {self.num_sequences}")
    
    # Boolean flags
    if getattr(self, 'use_gpu', True):
        options.append("--gpu")
    
    return " ".join(options)
```

### Environment-Specific Handling

```python
def generate_script(self, script_path: str) -> str:
    # Different commands for different environments
    if self.environment == "Boltz2Env":
        main_command = "boltz predict"
    elif self.environment == "ligandmpnn_env":
        main_command = "cd /path/to/LigandMPNN && python run.py"
    else:
        main_command = "python tool_script.py"
    
    script_content = f"""#!/bin/bash
echo "Environment: {self.environment}"
{main_command} {options}
"""
    return script_content
```

## Best Practices

### 1. Input Validation

```python
def validate_params(self):
    """Always validate all parameters thoroughly."""
    # Required parameters
    if not self.required_param:
        raise ValueError("required_param is mandatory")
    
    # Type validation
    if not isinstance(self.numeric_param, (int, float)):
        raise ValueError("numeric_param must be a number")
    
    # Range validation
    if self.numeric_param <= 0:
        raise ValueError("numeric_param must be positive")
    
    # Enum validation
    valid_options = ["option1", "option2", "option3"]
    if self.choice_param not in valid_options:
        raise ValueError(f"choice_param must be one of: {valid_options}")
    
    # File existence (for direct file inputs)
    if isinstance(self.file_param, str) and self.file_param.startswith('/'):
        if not os.path.exists(self.file_param):
            raise ValueError(f"File not found: {self.file_param}")
```

### 2. Robust Input Configuration

```python
def configure_inputs(self, pipeline_folders: Dict[str, str]):
    """Handle multiple input types gracefully."""
    self.folders = pipeline_folders
    
    # Store original input for reference
    self.original_input = self.input_param
    
    # Type-specific handling with clear error messages
    if self.input_is_tool_output:
        self._configure_tool_output_input()
    elif isinstance(self.input_param, str):
        self._configure_string_input()
    elif isinstance(self.input_param, list):
        self._configure_list_input()
    else:
        raise ValueError(
            f"Unsupported input type: {type(self.input_param)}. "
            f"Expected str, list, or ToolOutput"
        )

def _configure_tool_output_input(self):
    """Configure input from previous tool output."""
    tool_output: ToolOutput = self.input_param
    
    # Try different output types in order of preference
    for output_type in ["structures", "pdbs", "sequences"]:
        files = tool_output.get_output_files(output_type)
        if files:
            self.input_files = files
            self.input_type = output_type
            break
    else:
        available = list(tool_output._output_files.keys())
        raise ValueError(
            f"No compatible outputs from {tool_output.tool_type}. "
            f"Available: {available}"
        )
    
    # Add dependency
    self.dependencies.append(tool_output.config)
```

### 3. Flexible Output Definition

```python
def get_output_files(self) -> Dict[str, List[str]]:
    """Define outputs flexibly based on configuration."""
    outputs = {}
    
    # Base outputs (always present)
    outputs["output_folder"] = [self.output_folder]
    
    # Conditional outputs based on tool configuration
    if self.produces_structures():
        outputs["structures"] = self._get_structure_outputs()
        outputs["pdbs"] = outputs["structures"]  # Alias
    
    if self.produces_sequences():
        outputs["sequences"] = self._get_sequence_outputs()
    
    if self.produces_analysis():
        outputs["datasheets"] = self._get_analysis_outputs()
    
    # Tool-specific outputs
    outputs.update(self._get_tool_specific_outputs())
    
    return outputs

def produces_structures(self) -> bool:
    """Check if this tool configuration produces structures."""
    return getattr(self, 'output_format', 'none') in ['pdb', 'mmcif']

def _get_structure_outputs(self) -> List[str]:
    """Get structure output paths based on expected patterns."""
    if self.single_structure:
        return [os.path.join(self.output_folder, f"{self.job_name}.pdb")]
    else:
        return [os.path.join(self.output_folder, "structures")]  # Folder
```

### 4. Comprehensive Configuration Display

```python
def get_config_display(self) -> List[str]:
    """Provide detailed configuration info for users."""
    config_lines = super().get_config_display()
    
    # Input information
    if self.input_is_tool_output:
        config_lines.append(f"Input: {self.input_param.tool_type} output")
    else:
        config_lines.append(f"Input: {self.input_param}")
    
    # Key parameters
    config_lines.extend([
        f"Parameter 1: {self.param1}",
        f"Parameter 2: {self.param2}",
        f"Output format: {self.output_format}"
    ])
    
    # Conditional information
    if hasattr(self, 'optional_param') and self.optional_param:
        config_lines.append(f"Optional: {self.optional_param}")
    
    return config_lines
```

## Examples

### Simple File Processing Tool

```python
class FileProcessor(BaseConfig):
    TOOL_NAME = "FileProcessor"
    DEFAULT_ENV = "ProteinEnv"
    
    def __init__(self, input_file: str, format_type: str = "csv", **kwargs):
        self.input_file = input_file
        self.format_type = format_type
        super().__init__(**kwargs)
    
    def validate_params(self):
        if self.format_type not in ["csv", "json", "xml"]:
            raise ValueError("format_type must be csv, json, or xml")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders
        self.input_path = os.path.join(pipeline_folders["data"], self.input_file)
        if not os.path.exists(self.input_path):
            raise ValueError(f"File not found: {self.input_path}")
    
    def generate_script(self, script_path: str) -> str:
        processor_py = os.path.join(self.folders["HelpScripts"], "process_file.py")
        output_file = os.path.join(self.output_folder, f"processed.{self.format_type}")
        
        return f"""#!/bin/bash
echo "Processing file: {self.input_file}"
python {processor_py} "{self.input_path}" "{output_file}" --format {self.format_type}
echo "File processing completed"
"""
    
    def get_output_files(self) -> Dict[str, List[str]]:
        return {
            "processed_file": [os.path.join(self.output_folder, f"processed.{self.format_type}")],
            "datasheets": [os.path.join(self.output_folder, f"processed.{self.format_type}")]
        }
```

### Complex Multi-Input Tool

```python
class ComplexAnalyzer(BaseConfig):
    TOOL_NAME = "ComplexAnalyzer"
    DEFAULT_ENV = "analysis_env"
    
    def __init__(self, 
                 structures: Union[ToolOutput, str], 
                 sequences: Union[ToolOutput, str],
                 parameters: Dict[str, Any] = None,
                 **kwargs):
        self.structures = structures
        self.sequences = sequences
        self.parameters = parameters or {}
        
        self.structures_is_tool_output = isinstance(structures, ToolOutput)
        self.sequences_is_tool_output = isinstance(sequences, ToolOutput)
        
        super().__init__(**kwargs)
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders
        
        # Configure structure inputs
        if self.structures_is_tool_output:
            struct_files = self.structures.get_output_files("structures")
            if not struct_files:
                raise ValueError("No structures found in tool output")
            self.structure_path = struct_files[0]
            self.dependencies.append(self.structures.config)
        else:
            self.structure_path = os.path.join(pipeline_folders["PDBs"], self.structures)
        
        # Configure sequence inputs
        if self.sequences_is_tool_output:
            seq_files = self.sequences.get_output_files("sequences")
            if not seq_files:
                raise ValueError("No sequences found in tool output") 
            self.sequence_path = seq_files[0]
            if self.sequences.config not in self.dependencies:
                self.dependencies.append(self.sequences.config)
        else:
            self.sequence_path = os.path.join(pipeline_folders["data"], self.sequences)
    
    def generate_script(self, script_path: str) -> str:
        analyzer_py = os.path.join(self.folders["HelpScripts"], "analyze_complex.py")
        config_file = os.path.join(self.output_folder, "analysis_config.json")
        
        # Create configuration file
        config_content = f"""
import json
config = {{
    "structures": "{self.structure_path}",
    "sequences": "{self.sequence_path}",
    "parameters": {json.dumps(self.parameters)},
    "output_dir": "{self.output_folder}"
}}
with open("{config_file}", "w") as f:
    json.dump(config, f, indent=2)
"""
        
        return f"""#!/bin/bash
echo "Creating analysis configuration"
python -c '{config_content}'

echo "Running complex analysis"
python {analyzer_py} --config "{config_file}"

echo "Analysis completed"
"""
    
    def get_output_files(self) -> Dict[str, List[str]]:
        return {
            "analysis_results": [os.path.join(self.output_folder, "results.csv")],
            "plots": [os.path.join(self.output_folder, "plots")],
            "datasheets": [
                os.path.join(self.output_folder, "results.csv"),
                os.path.join(self.output_folder, "summary.json")
            ]
        }
```

## Testing Your Tools

### Quick Test Script

```python
#!/usr/bin/env python3
from PipelineScripts import Pipeline, MyTool

# Test tool in isolation
pipeline = Pipeline("TestPipeline", "TestJob", "Testing MyTool")

# Add your tool
tool_output = pipeline.add(MyTool(
    input_param="test_input.txt",
    specific_param=5
))

# Check configuration
print("Tool configured successfully!")
print("Expected outputs:", tool_output.get_output_files())

# Generate pipeline
try:
    script_path = pipeline.save()
    print(f"Pipeline script generated: {script_path}")
except Exception as e:
    print(f"Error: {e}")
```

This comprehensive guide should help you create robust, well-integrated tools for the pipeline system. Remember to always follow the established patterns for consistency and maintainability!