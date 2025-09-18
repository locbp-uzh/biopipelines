"""
PyMOL visualization configuration for creating molecular sessions.

Creates PyMOL session files (.pse) from structure outputs with alignment,
coloring based on datasheet metrics, and multi-tool structure combination.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class PyMOL(BaseConfig):
    """
    PyMOL visualization configuration for creating molecular sessions.
    
    Creates PyMOL session files from structure outputs with proper alignment,
    coloring based on datasheet columns, and support for multi-tool combinations.
    """
    
    # Tool identification
    TOOL_NAME = "Pymol"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "none", "memory": "8GB", "time": "2:00:00"}
    
    def __init__(self,
                 structures: Union[ToolOutput, List[Tuple[ToolOutput, str]]] = None,
                 color_by: Optional[str] = None,
                 reference_structure: Optional[str] = None,
                 alignment: str = "align",
                 session_name: Optional[str] = None,
                 **kwargs):
        """
        Initialize PyMOL visualization configuration.
        
        Args:
            structures: Structure source(s). Either:
                - Single ToolOutput for simple visualization
                - List of (ToolOutput, prefix) tuples for multi-tool combination
            color_by: Datasheet column reference for coloring (e.g., "tool.datasheets.analysis.plddt")
            reference_structure: Structure prefix to use as alignment reference (default: first structure)
            alignment: Alignment method - "align", "cealign", "super", "alignto", or "none"
            session_name: Name for output PyMOL session file (default: auto-generated)
            **kwargs: Additional configuration parameters
        """
        
        # Handle structures input format
        self.structure_sources = self._normalize_structure_inputs(structures)
        self.color_by = color_by
        self.reference_structure = reference_structure
        self.alignment = alignment
        self.session_name = session_name or "pymol_session"
        
        # Validate parameters
        self._validate_alignment_method()
        
        super().__init__(**kwargs)
    
    def _normalize_structure_inputs(self, structures):
        """Normalize structure inputs to consistent format."""
        if structures is None:
            raise ValueError("structures parameter is required")
        
        if isinstance(structures, ToolOutput):
            return [(structures, "structures")]
        elif isinstance(structures, list):
            normalized = []
            for item in structures:
                if isinstance(item, tuple) and len(item) == 2:
                    tool_output, prefix = item
                    if not isinstance(tool_output, ToolOutput):
                        raise ValueError(f"Expected ToolOutput in tuple, got {type(tool_output)}")
                    normalized.append((tool_output, prefix))
                else:
                    raise ValueError("List items must be (ToolOutput, prefix) tuples")
            return normalized
        else:
            raise ValueError(f"structures must be ToolOutput or list of (ToolOutput, prefix) tuples, got {type(structures)}")
    
    def _validate_alignment_method(self):
        """Validate alignment method parameter."""
        valid_methods = ["align", "cealign", "super", "alignto", "none"]
        if self.alignment not in valid_methods:
            raise ValueError(f"alignment must be one of {valid_methods}, got '{self.alignment}'")
    
    def validate_params(self):
        """Validate PyMOL-specific parameters."""
        if not self.structure_sources:
            raise ValueError("At least one structure source must be provided")
        
        # Validate color_by format if provided
        if self.color_by:
            parts = self.color_by.split('.')
            if len(parts) < 4 or parts[1] != "output" or parts[2] != "datasheets":
                raise ValueError(
                    "color_by must follow format 'tool.datasheets.sheet.column', "
                    f"got '{self.color_by}'"
                )
        
        # Validate reference_structure exists in prefixes if specified
        if self.reference_structure:
            available_prefixes = [prefix for _, prefix in self.structure_sources]
            if self.reference_structure not in available_prefixes:
                raise ValueError(
                    f"reference_structure '{self.reference_structure}' not found in structure prefixes: "
                    f"{available_prefixes}"
                )
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from pipeline context."""
        self.input_sources = {}
        
        # Configure structure inputs
        for i, (tool_output, prefix) in enumerate(self.structure_sources):
            key = f"structures_{i}"
            if hasattr(tool_output, 'output_folder'):
                self.input_sources[key] = {
                    'type': 'tool_output',
                    'path': tool_output.output_folder,
                    'prefix': prefix,
                    'tool_output': tool_output
                }
            else:
                raise ValueError(f"Structure source {i} does not have output_folder")
        
        # Set output folder
        self.output_folder = os.path.join(
            pipeline_folders['job_folder'],
            f"{self.execution_order:03d}_{self.TOOL_NAME}"
        )
    
    def predict_output_structure(self) -> StandardizedOutput:
        """Predict the output structure for PyMOL visualization."""
        
        # PyMOL session file
        session_file = os.path.join(self.output_folder, f"{self.session_name}.pse")
        
        # Count total structures across all sources
        total_structures = 0
        for tool_output, _ in self.structure_sources:
            if hasattr(tool_output, 'structure_count'):
                total_structures += tool_output.structure_count
            else:
                total_structures += 1  # Assume at least one structure
        
        return StandardizedOutput(
            structures_folder=self.output_folder,
            structures_count=1,  # One PSE session file
            datasheets=[],  # PyMOL doesn't generate datasheets
            output_info={
                'session_file': session_file,
                'structure_sources': len(self.structure_sources),
                'total_input_structures': total_structures,
                'alignment_method': self.alignment,
                'color_scheme': self.color_by or 'default',
                'reference_structure': self.reference_structure
            }
        )
    
    def generate_script(self, script_path: str) -> str:
        """Generate bash script for PyMOL session creation."""
        
        # Create configuration for helper script
        config = {
            'structure_sources': [],
            'color_by': self.color_by,
            'reference_structure': self.reference_structure,
            'alignment': self.alignment,
            'session_name': self.session_name,
            'output_folder': self.output_folder
        }
        
        # Add structure source configurations
        for i, (tool_output, prefix) in enumerate(self.structure_sources):
            source_config = {
                'input_folder': self.input_sources[f'structures_{i}']['path'],
                'prefix': prefix,
                'tool_name': getattr(tool_output, 'tool_name', 'unknown')
            }
            config['structure_sources'].append(source_config)
        
        config_file = os.path.join(self.output_folder, 'pymol_config.json')
        session_file = os.path.join(self.output_folder, f"{self.session_name}.pse")
        
        script_content = f"""#!/bin/bash
set -e

# Create output directory
mkdir -p {self.output_folder}

# Write configuration
cat > {config_file} << 'EOF'
{json.dumps(config, indent=2)}
EOF

# Run PyMOL session creation
python /Users/gianluca/Desktop/PhD/Coding/Git/biopipelines/HelpScripts/pipe_pymol.py \\
    --config {config_file} \\
    --output {session_file}

echo "PyMOL session created: {session_file}"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def get_parameter_summary(self) -> List[str]:
        """Get human-readable summary of PyMOL parameters."""
        summary = [
            f"Structure sources: {len(self.structure_sources)}",
            f"Alignment: {self.alignment}"
        ]
        
        if self.color_by:
            summary.append(f"Color by: {self.color_by}")
        
        if self.reference_structure:
            summary.append(f"Reference: {self.reference_structure}")
        else:
            summary.append("Reference: first structure")
        
        summary.append(f"Session: {self.session_name}.pse")
        
        return summary
    
    def __str__(self) -> str:
        """String representation of the PyMOL configuration."""
        sources_str = ", ".join([f"{prefix}" for _, prefix in self.structure_sources])
        return f"Pymol(sources=[{sources_str}], alignment={self.alignment})"