"""
Core Pipeline class for automated workflow construction.

Manages tool sequencing, environment switching, dependency resolution,
and script generation for protein modeling pipelines.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union
from collections import defaultdict
from datetime import datetime

from .folders import FolderManager
try:
    from .base_config import BaseConfig, ToolOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput


class Pipeline:
    """
    Main pipeline orchestration class.
    
    Manages tools, dependencies, environment switching, and script generation
    for automated protein modeling workflows.
    """
    
    def __init__(self, pipeline_name: str, job_name: str, job_description: str="Description missing", debug: bool=False):
        if ' ' in job_name: job_name=job_name.replace(' ','_') #It will create issues at runtime otherwise
        
        self.pipeline_name = pipeline_name
        self.job_name = job_name
        self.job_description = job_description
        self.debug = debug

        self.folder_manager = FolderManager(pipeline_name,job_name,debug)
        self.folders = self.folder_manager.get_folders()
        
        # Tool management
        self.tools = []
        self.tool_outputs = []
        self.execution_order = 0
        
        # Environment tracking
        self.environments_used = set()
        self.environment_groups = defaultdict(list)
        
        # Script generation
        self.pipeline_script = ""
        self.slurm_script = ""
        self.scripts_generated = False

        # Resource management
        self.global_resources = {
            "gpu": "V100",
            "memory": "15GB", 
            "time": "24:00:00"
        }
    
    def add(self, tool_config: BaseConfig, env: str = None, **kwargs) -> ToolOutput:
        """
        Add a tool to the pipeline.
        
        Args:
            tool_config: Configured tool instance (e.g., RFdiffusion())
            env: Environment override
            **kwargs: Additional pipeline-specific options
            
        Returns:
            ToolOutput object with metadata and output information
        """
        if env is not None: tool_config.environment = env
        tool_config.resources.update(self.global_resources) #will warn if resources are insufficient
        
        # Set execution order and create step-numbered folder immediately
        self.execution_order += 1
        step_folder_name = f"{self.execution_order}_{tool_config.TOOL_NAME}"
        tool_output_folder = os.path.join(self.folders["output"], step_folder_name)
        os.makedirs(tool_output_folder, exist_ok=True)
        
        # Set pipeline context with unified folder structure
        tool_config.set_pipeline_context(
            self, self.execution_order, tool_output_folder
        )
        
        # Configure inputs immediately so tool outputs are properly set
        tool_config.configure_inputs(self.folders)
        
        # Track environment usage
        self.environments_used.add(tool_config.environment)
        self.environment_groups[tool_config.environment].append(tool_config)
        
        # Add to pipeline
        self.tools.append(tool_config)
        
        # Create and configure output object
        tool_output = ToolOutput(tool_config)
        
        # Immediately populate expected outputs using pure path construction
        try:
            expected_outputs = tool_config.get_output_files()
            tool_output.update_outputs(expected_outputs)
        except Exception as e:
            # If get_output_files fails, tool may need dependencies resolved first
            # We'll populate it later during script generation
            print(f"Warning: Could not immediately populate outputs for {tool_config.TOOL_NAME}: {e}")
        
        self.tool_outputs.append(tool_output)
        
        return tool_output
    
    def validate_pipeline(self) -> bool:
        """
        Validate entire pipeline for consistency and compatibility.
        
        Returns:
            True if pipeline is valid
            
        Raises:
            ValueError: If pipeline has issues
        """
        if not self.tools:
            raise ValueError("Pipeline is empty")
        
        # Check tool dependencies
        for i, tool in enumerate(self.tools):
            for dep in tool.dependencies:
                if isinstance(dep, BaseConfig):
                    if dep not in self.tools[:i]:
                        raise ValueError(
                            f"Tool {tool.TOOL_NAME} depends on {dep.TOOL_NAME} "
                            f"which hasn't been added yet"
                        )
        
        # Check environment compatibility
        for env in self.environments_used:
            env_tools = [t.TOOL_NAME for t in self.environment_groups[env]]
            print(f"Environment {env}: {', '.join(env_tools)}")
        
        return True
    
    def save(self) -> str:
        """
        Generate and save pipeline execution script.
        
        Returns:
            Path to generated pipeline script
        """
        if not self.tools:
            raise ValueError("Cannot save empty pipeline")
        
        # Print tool outputs with execution order
        for i, tool_output in enumerate(self.tool_outputs, 1):
            print("="*30+f"{i}.{tool_output.config.TOOL_NAME}"+"="*30)
            print(tool_output.output)
        
        print("="*30+"Pipeline"+"="*30)
        self.validate_pipeline()
        
        # Generate pipeline script
        script_path = os.path.join(
            self.folders["runtime"], 
            "pipeline.sh"
        )
        
        script_content = self._generate_pipeline_script()
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        self.pipeline_script = script_path
        self.scripts_generated = True
        
        # Export tool outputs for potential reuse with LoadOutput
        self._export_tool_outputs()
        
        print(f"Pipeline saved to: {script_path}")
        return script_path
    
    def _generate_pipeline_script(self) -> str:
        """Generate unified pipeline script following notebook pattern."""
        # Pipeline folders already created during initialization
        # Just need to configure tools and generate scripts
        
        # Process tools in order, setting up outputs for dependencies first
        for i, tool in enumerate(self.tools, 1):
            # Ensure all previous tools have their outputs set up
            for j in range(i):  # 0 to i-1 (all previous tools)
                prev_tool = self.tools[j]
                prev_tool_output = self.tool_outputs[j]
                
                # Update previous tool outputs if not already done
                if not prev_tool_output._output_files:
                    try:
                        expected_outputs = prev_tool.get_output_files()
                        prev_tool_output.update_outputs(expected_outputs)
                    except Exception as e:
                        # If we can't get outputs, create empty ones so dependencies don't fail
                        print(f"Warning: Could not get outputs from {prev_tool.TOOL_NAME}: {e}")
                        prev_tool_output.update_outputs({})
            
            # Now configure this tool's inputs (dependencies should have outputs now)
            tool.configure_inputs(self.folders)
            
            # Set up this tool's outputs
            try:
                expected_outputs = tool.get_output_files()
                self.tool_outputs[i-1].update_outputs(expected_outputs)
            except Exception as e:
                # If we can't get outputs yet, leave empty for now
                print(f"Warning: Could not get outputs from {tool.TOOL_NAME} during pipeline setup: {e}")
        
        # Generate configuration display
        config_lines = [
            f'echo "Pipeline: {self.pipeline_name}"',
            f'echo "Job name: {self.job_name}"',
            f'echo "Description: {self.job_description}"',
            'echo'
        ]
        
        # Add tool configurations
        for i, tool in enumerate(self.tools, 1):
            config_lines.extend([
                f'echo "Step {i}: {tool.TOOL_NAME}"',
                f'echo "  Environment: {tool.environment}"'
            ])
            
            # Add tool-specific parameters (will be implemented in each tool)
            tool_config = tool.get_config_display()
            for config_line in tool_config:
                config_lines.append(f'echo "  {config_line}"')
            
            config_lines.append('echo')
        
        # Create config script
        config_script = os.path.join(self.folders["runtime"], "config.sh")
        with open(config_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(config_lines))
        os.chmod(config_script, 0o755)
        
        # Build main pipeline script
        script_lines = [
            "#!/bin/bash",
            f"# Generated pipeline: {self.pipeline_name}",
            f"# Tools: {', '.join([t.TOOL_NAME for t in self.tools])}",
            "",
            "set -e  # Exit on any error",
            ""
        ]
        
        # Determine primary environment (most common)
        env_counts = {}
        for tool in self.tools:
            env_counts[tool.environment] = env_counts.get(tool.environment, 0) + 1
        primary_env = max(env_counts, key=env_counts.get)
        
        # Activate primary environment
        script_lines.extend([
            f"source activate {primary_env}",
            "echo",
            "echo Configuration",
            f"{config_script} | tee {os.path.join(self.folders['output'], f'{self.pipeline_name}_config.txt')}",
            "echo"
        ])
        
        # Add each tool execution
        current_env = primary_env
        
        for i, tool in enumerate(self.tools, 1):
            # Environment switch if needed
            if tool.environment != current_env:
                script_lines.extend([
                    f"source activate {tool.environment}",
                    "echo"
                ])
                current_env = tool.environment
            
            # Generate tool-specific script
            tool_script_path = os.path.join(self.folders["runtime"], f"{i}_{tool.TOOL_NAME.lower()}.sh")
            tool_script_content = tool.generate_script(tool_script_path)
            
            # Write tool script
            with open(tool_script_path, 'w') as f:
                f.write(tool_script_content)
            os.chmod(tool_script_path, 0o755)
            
            # Update tool output files now that paths are correct
            try:
                expected_outputs = tool.get_output_files()
                # Find the corresponding ToolOutput and update it
                for tool_output in self.tool_outputs:
                    if tool_output.config == tool:
                        tool_output.update_outputs(expected_outputs)
                        break
            except Exception as e:
                # Continue if output files can't be determined yet
                # This is expected - files don't exist until execution
                print(f"Debug: Could not determine output files for {tool.TOOL_NAME} during script generation: {e}")
            
            # Add tool execution following notebook pattern
            log_file = os.path.join(self.folders["output"], f"_{i}_{tool.TOOL_NAME.lower()}.log")
            script_lines.extend([
                f"echo {tool.TOOL_NAME}",
                f"{tool_script_path} | tee {log_file}",
                "echo"
            ])
        
        # Final steps
        script_lines.extend([
            "echo",
            "echo Job done",
            f"echo Results in:",
            f"echo {self.folders['output']}"
        ])
        
        return "\n".join(script_lines)
    
    def slurm(self, email: str = "", gpu: str = None, memory: str = None, 
              time: str = None, **slurm_options):
        """
        Generate SLURM job submission script.
        
        Args:
            user: Username for job submission
            gpu: GPU type requirement
            memory: Memory requirement
            time: Time limit
            **slurm_options: Additional SLURM options
            
        Returns:
            SLURM script content
        """
        if not self.scripts_generated:
            raise ValueError("Must call save() before generating SLURM script")
        
        email_line = "" if email == "" else f"""
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={email}"""
        
        # Use provided resources or defaults
        resources = self.global_resources.copy()
        if gpu: resources["gpu"] = gpu
        if memory: resources["memory"] = memory  
        if time: resources["time"] = time
        
        # Determine GPU constraint
        if resources["gpu"] == "high-memory":
            gpu_line = "#SBATCH --gpus=1\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB\""
        elif resources["gpu"] == "gpu":
            gpu_line = "#SBATCH --gpus=1"
        else:
            gpu_line = f"#SBATCH --gpus={resources['gpu']}:1"
        
        # Convert memory to MB if needed
        memory_mb = resources["memory"]
        if memory_mb.endswith("GB"):
            memory_mb = str(int(float(memory_mb[:-2]) * 1000))
        
        # Generate SLURM script
        slurm_content = f"""#!/usr/bin/bash
{gpu_line}
#SBATCH --mem={memory_mb}
#SBATCH --time={resources["time"]}
#SBATCH --output=job.out
#SBATCH --begin=now+0hour
{email_line}

# Check if nvidia-smi is available
if ! command -v nvidia-smi &> /dev/null
then
    echo "Could not load GPU correctly: nvidia-smi could not be found"
    exit 1
fi

gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader)
echo "GPU Type: $gpu_type"

module load mamba

# Execute pipeline
{self.pipeline_script}
"""
        
        # Save SLURM script
        slurm_path = os.path.join(
            self.folders["runtime"],
            f"slurm.sh"
        )
        
        with open(slurm_path, 'w') as f:
            f.write(slurm_content)
        os.chmod(slurm_path, 0o755)
        
        self.slurm_script = slurm_path
        
        print("="*30+"Job"+"="*30)
        print(f"{self.pipeline_name}: {self.job_name}")
        print("="*30+"Slurm Script"+"="*30)
        # Print line by line to ensure proper formatting
        for line in slurm_content.split('\n'):
            if line != "": 
                print(line)
    
    def resources(self, gpu: str = None, memory: str = None, time: str = None):
        if gpu: self.global_resources["gpu"] = gpu #T4, V100, A100, gpu, high-memory
        if memory: self.global_resources["memory"] = memory
        if time: self.global_resources["time"] = time
    
    def get_tool_outputs(self, tool_type: str = None) -> List[ToolOutput]:
        """
        Get outputs from tools in pipeline.
        
        Args:
            tool_type: Filter by tool type (optional)
            
        Returns:
            List of ToolOutput objects
        """
        if tool_type is None:
            return self.tool_outputs
        return [out for out in self.tool_outputs if out.tool_type == tool_type]
    
    def summary(self) -> str:
        """Get pipeline summary."""
        summary = [
            f"Pipeline: {self.pipeline_name}",
            f"Tools: {len(self.tools)}",
            f"Environments: {len(self.environments_used)}",
            ""
        ]
        
        for i, tool in enumerate(self.tools, 1):
            summary.append(
                f"{i}. {tool.TOOL_NAME} ({tool.environment}) -> {tool.job_name}"
            )
        
        return "\n".join(summary)
    
    def _export_tool_outputs(self):
        """
        Export tool output metadata to JSON files for potential reuse with LoadOutput.
        
        Creates a tool_outputs/ directory with JSON files for each tool containing
        complete output structure, configuration, and execution metadata.
        """
        from datetime import datetime
        
        # Create tool_outputs directory
        tool_outputs_dir = os.path.join(self.folders["output"], "tool_outputs")
        os.makedirs(tool_outputs_dir, exist_ok=True)
        
        exported_count = 0
        
        for i, (tool, tool_output) in enumerate(zip(self.tools, self.tool_outputs), 1):
            try:
                # Get current output structure
                output_structure = tool.get_output_files()
                
                # Convert DatasheetInfo objects to dictionaries for JSON serialization
                if 'datasheets' in output_structure:
                    from .base_config import DatasheetInfo
                    serializable_datasheets = {}
                    for name, datasheet_info in output_structure['datasheets'].items():
                        if isinstance(datasheet_info, DatasheetInfo):
                            serializable_datasheets[name] = datasheet_info.to_dict()
                        else:
                            serializable_datasheets[name] = datasheet_info
                    output_structure['datasheets'] = serializable_datasheets
                
                # Build complete tool metadata
                tool_metadata = {
                    "tool_name": tool.TOOL_NAME,
                    "tool_class": tool.__class__.__name__,
                    "job_name": tool.job_name,
                    "execution_order": i,
                    "environment": tool.environment,
                    "output_structure": output_structure,
                    "configuration": {
                        "tool_parameters": tool.to_dict(),
                        "resources": tool.resources.copy(),
                        "pipeline_context": {
                            "pipeline_name": self.pipeline_name,
                            "pipeline_job_name": self.job_name,
                            "pipeline_description": self.job_description
                        }
                    },
                    "export_metadata": {
                        "export_time": datetime.now().isoformat(),
                        "pipeline_version": "1.0",
                        "exported_by": "pipeline.save()"
                    }
                }
                
                # Add execution metadata placeholder (will be filled at runtime)
                tool_metadata["execution_metadata"] = {
                    "completion_time": None,
                    "execution_success": None,
                    "log_files": [],
                    "runtime_notes": "Execution metadata will be available after pipeline execution"
                }
                
                # Generate filename
                output_filename = f"{self.job_name}_{i}_{tool.TOOL_NAME}_output.json"
                output_path = os.path.join(tool_outputs_dir, output_filename)
                
                # Save to JSON file
                with open(output_path, 'w') as f:
                    json.dump(tool_metadata, f, indent=2)
                
                exported_count += 1
                
            except Exception as e:
                print(f"Warning: Could not export output metadata for {tool.TOOL_NAME}: {e}")
                continue
        
        if exported_count > 0:
            print(f"Exported {exported_count} tool output metadata files to: {tool_outputs_dir}")
        else:
            print("Warning: No tool output metadata could be exported")