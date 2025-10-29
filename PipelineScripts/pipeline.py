"""
Core Pipeline class for automated workflow construction.

Manages tool sequencing, environment switching, dependency resolution,
and script generation for protein modeling pipelines.
"""

import os
import json
import getpass
import inspect
import shutil
import contextvars
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

# Module-level context variable to track active pipeline for auto-registration
_active_pipeline: contextvars.ContextVar[Optional['Pipeline']] = contextvars.ContextVar('_active_pipeline', default=None)

# Email mapping for auto email functionality
EMAILS = {
    'rebeca': 'rebecca.andrews@uzh.ch',
    'dkossm': 'dorothea.kossmann@chem.uzh.ch',
    'hlaemm': 'henriette.laemmermann@chem.uzh.ch',
    'gquarg': 'gianluca.quargnali@chem.uzh.ch',
    'pabriv': 'pablo.riverafuentes@uzh.ch',
    'clsonn': 'clarissa.sonnemann@chem.uzh.ch',
    'crstev': 'craig.steven@chem.uzh.ch'
}


class Pipeline:
    """
    Main pipeline orchestration class.
    
    Manages tools, dependencies, environment switching, and script generation
    for automated protein modeling workflows.
    """
    
    def __init__(self, pipeline_name: str, job_name: str, job_description: str="Description missing", debug: bool=False, shared: bool=False):
        """
        Initialize a new pipeline instance.

        Args:
            pipeline_name: Name of the pipeline (used for output folders)
            job_name: Name of the specific job (used for output folders)
            job_description: Optional description of what this job does
            debug: If True, runs in debug mode without creating folders
            shared: If True, uses shared/public user for folder structure (default: False)
        """
        if ' ' in job_name: job_name=job_name.replace(' ','_') #It will create issues at runtime otherwise
        
        self.pipeline_name = pipeline_name
        self.job_name = job_name
        self.job_description = job_description
        self.debug = debug
        self.shared = shared

        self.folder_manager = FolderManager(pipeline_name, job_name, debug, shared)
        self.folders = self.folder_manager.get_folders()
        
        # Tool management
        self.tools = []
        self.tool_outputs = []
        self.execution_order = 0
        
        # Suffix for tool folder naming
        self.current_suffix = ""
        
        # Environment tracking
        self.environments_used = set()
        self.environment_groups = defaultdict(list)
        
        # Script generation
        self.pipeline_script = ""
        self.slurm_script = ""
        self.scripts_generated = False

        # Resource management
        self.global_resources = {
            "gpu": None,  # No GPU by default
            "memory": "15GB",
            "time": "24:00:00"
        }

        # Context manager state
        self._explicit_save_called = False  # Track if Save() was explicitly called

    def __enter__(self):
        """
        Enter context manager - set this pipeline as the active context.

        Returns:
            self (Pipeline instance)

        Raises:
            RuntimeError: If another pipeline context is already active (nested contexts not allowed)
        """
        current_pipeline = _active_pipeline.get()
        if current_pipeline is not None:
            raise RuntimeError(
                f"Cannot nest Pipeline contexts. "
                f"Pipeline '{current_pipeline.pipeline_name}' is already active."
            )
        _active_pipeline.set(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context manager - clear the active pipeline context and auto-submit to SLURM.

        Automatically calls self.slurm() on exit unless:
        - An exception occurred during execution
        - Save() was explicitly called (which saves without submitting)

        Args:
            exc_type: Exception type if an exception occurred
            exc_val: Exception value if an exception occurred
            exc_tb: Exception traceback if an exception occurred

        Returns:
            False (don't suppress exceptions)
        """
        try:
            # Only auto-submit if no exception and Save() wasn't explicitly called
            if exc_type is None and not self._explicit_save_called:
                self.slurm()
        finally:
            # Always clear the active pipeline context
            _active_pipeline.set(None)

        return False  # Don't suppress exceptions

    @staticmethod
    def get_active_pipeline() -> Optional['Pipeline']:
        """
        Get the currently active pipeline context, if any.

        Returns:
            Pipeline instance if within a context manager, None otherwise
        """
        return _active_pipeline.get()

    def set_suffix(self, suffix: str = ""):
        """
        Set suffix for subsequent tool folder names.
        
        Args:
            suffix: Suffix to append to tool folder names (e.g., "single_point")
                   Empty string clears the suffix.
        
        Example:
            pipeline.set_suffix("single_point")
            # Next tools will create folders like: 004_MutationComposer_single_point/
        """
        self.current_suffix = suffix
    
    def add(self, tool_config: BaseConfig, env: str = None, **kwargs):
        """
        Add a tool to the pipeline.

        Args:
            tool_config: Configured tool instance (e.g., RFdiffusion())
            env: Environment override
            **kwargs: Additional pipeline-specific options

        Returns:
            ToolAwareOutput object with standardized output interface and .tool access
        """
        if env is not None: tool_config.environment = env

        # Merge pipeline resources with tool resources (pipeline settings take precedence if not None)
        for key, value in self.global_resources.items():
            if value is not None:  # Only override if pipeline has a non-None value
                tool_config.resources[key] = value
        
        # Set execution order and create step-numbered folder immediately
        self.execution_order += 1
        if self.current_suffix:
            step_folder_name = f"{self.execution_order:03d}_{tool_config.TOOL_NAME}_{self.current_suffix}"
        else:
            step_folder_name = f"{self.execution_order:03d}_{tool_config.TOOL_NAME}"
        tool_output_folder = os.path.join(self.folders["output"], step_folder_name)
        os.makedirs(tool_output_folder, exist_ok=True)

        # Set pipeline context with unified folder structure
        tool_config.set_pipeline_context(
            self, self.execution_order, tool_output_folder, self.current_suffix
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

        return tool_output.output

    def _auto_register(self, tool_config: BaseConfig) -> ToolOutput:
        """
        Internal method for auto-registering tools created within a Pipeline context.

        This is called automatically when a tool is instantiated within a
        'with Pipeline() as pipeline:' context. It performs the same operations
        as add() but is designed to be called from BaseConfig.__new__.

        Args:
            tool_config: Configured tool instance

        Returns:
            ToolOutput object for the registered tool
        """
        # Merge pipeline resources with tool resources
        for key, value in self.global_resources.items():
            if value is not None:
                tool_config.resources[key] = value

        # Set execution order and create step-numbered folder
        self.execution_order += 1
        if self.current_suffix:
            step_folder_name = f"{self.execution_order:03d}_{tool_config.TOOL_NAME}_{self.current_suffix}"
        else:
            step_folder_name = f"{self.execution_order:03d}_{tool_config.TOOL_NAME}"
        tool_output_folder = os.path.join(self.folders["output"], step_folder_name)
        os.makedirs(tool_output_folder, exist_ok=True)

        # Set pipeline context
        tool_config.set_pipeline_context(
            self, self.execution_order, tool_output_folder, self.current_suffix
        )

        # Configure inputs immediately
        tool_config.configure_inputs(self.folders)

        # Track environment usage
        self.environments_used.add(tool_config.environment)
        self.environment_groups[tool_config.environment].append(tool_config)

        # Add to pipeline
        self.tools.append(tool_config)

        # Create and configure output object
        tool_output = ToolOutput(tool_config)

        # Populate expected outputs
        try:
            expected_outputs = tool_config.get_output_files()
            tool_output.update_outputs(expected_outputs)
        except Exception as e:
            print(f"Warning: Could not immediately populate outputs for {tool_config.TOOL_NAME}: {e}")

        self.tool_outputs.append(tool_output)

        return tool_output

    def validate_pipeline(self, debug) -> bool:
        """
        Validate entire pipeline for consistency and compatibility.
        
        Returns:
            True if pipeline is valid
            
        Raises:
            ValueError: If pipeline has issues
        """
        if not self.tools:
            raise ValueError("Pipeline is empty")
        if not debug: return
        
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
    
    def save(self, debug=False) -> str:
        """
        Generate and save pipeline execution script.
        
        Returns:
            Path to generated pipeline script
        """
        if not self.tools:
            raise ValueError("Cannot save empty pipeline")
        
        if debug:
            # Print tool outputs with execution order
            for i, tool_output in enumerate(self.tool_outputs, 1):
                print("="*30+f"{i}.{tool_output.config.TOOL_NAME}"+"="*30)
                print(tool_output.output)
            
            print("="*30+"Pipeline"+"="*30)
        self.validate_pipeline(debug)
        
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

        # Save the original pipeline Python script to runtime folder
        self._save_original_pipeline_script()

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
            "umask 002  # Make all files group-writable by default",
            ""
        ]
        
        # Determine primary environment (most common)
        env_counts = {}
        for tool in self.tools:
            env_counts[tool.environment] = env_counts.get(tool.environment, 0) + 1
        primary_env = max(env_counts, key=env_counts.get)
        
        # Activate primary environment (if one is needed)
        if primary_env is not None:
            script_lines.extend([
                f"source activate {primary_env}",
                "echo"
            ])

        script_lines.extend([
            "echo Configuration",
            f"{config_script} | tee {os.path.join(self.folders['output'], f'{self.pipeline_name}_config.txt')}",
            "echo"
        ])
        
        # Add each tool execution
        current_env = primary_env
        
        for i, tool in enumerate(self.tools, 1):
            # Environment switch if needed
            if tool.environment != current_env and tool.environment is not None:
                script_lines.extend([
                    f"source activate {tool.environment}",
                    "echo"
                ])
                current_env = tool.environment
            
            # Generate tool-specific script
            tool_script_path = os.path.join(self.folders["runtime"], f"{i:03d}_{tool.TOOL_NAME.lower()}.sh")
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
            log_file = os.path.join(self.folders["logs"], f"{i:03d}_{tool.TOOL_NAME.lower()}.log")
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
    
    def slurm(self, email: str = "auto", gpu: str = None, memory: str = None,
              time: str = None, **slurm_options):
        """
        Generate SLURM job submission script.

        Args:
            email: Email for job notifications. Options:
                - "auto": Automatically detect email from current username (default)
                - "": No email notifications
                - specific email: Use provided email address
            gpu: GPU type requirement
            memory: Memory requirement
            time: Time limit
            **slurm_options: Additional SLURM options passed as param=value pairs.
                Each parameter will be converted to #SBATCH --param=value.
                Examples:
                - nodes=1 → #SBATCH --nodes=1
                - ntasks_per_node=1 → #SBATCH --ntasks-per-node=1
                - cpus_per_task=32 → #SBATCH --cpus-per-task=32
                - partition="gpu" → #SBATCH --partition=gpu

        Returns:
            SLURM script content
        """
        if not self.scripts_generated:
            self.save()

        # Handle auto email detection
        if email == "auto":
            current_user = getpass.getuser()
            if current_user in EMAILS:
                email = EMAILS[current_user]
                print(f"Auto-detected email: {email} (user: {current_user})")
            else:
                print(f"Warning: No email mapping found for user '{current_user}'. Disabling email notifications.")
                email = ""

        email_line = "" if email == "" else f"""
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={email}"""
        
        # Use provided resources or defaults
        resources = self.global_resources.copy()
        if gpu: resources["gpu"] = gpu
        if memory: resources["memory"] = memory  
        if time: resources["time"] = time
        
        # Determine GPU constraint
        gpu_spec = resources["gpu"]

        if gpu_spec is None or gpu_spec == "none":
            # No GPU requested
            gpu_line = ""
        elif gpu_spec == "high-memory":
            gpu_line = "#SBATCH --gpus=1\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB\""
        elif gpu_spec == "gpu":
            gpu_line = "#SBATCH --gpus=1"
        elif gpu_spec.startswith("!"):
            # Exclusion syntax (e.g., "!T4")
            # Use multiple memory constraints to exclude T4 (16GB) by requesting 32GB or 80GB
            excluded_model = gpu_spec[1:]  # Remove the "!" prefix
            if excluded_model.upper() == "T4":
                # T4 has 16GB, so request 32GB or 80GB to exclude it
                gpu_line = "#SBATCH --gpus=1\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB\""
            else:
                # For other exclusions, try the tilde constraint syntax
                gpu_line = f"#SBATCH --gpus=1\n#SBATCH --constraint=\"~GPU{excluded_model}\""
        elif gpu_spec in ["32GB", "80GB"] or "|" in gpu_spec:
            # Memory-based GPU specification
            if "|" in gpu_spec:
                # Multiple memory options (e.g., "32GB|80GB")
                memory_options = gpu_spec.split("|")
                constraint_parts = [f"GPUMEM{mem}" for mem in memory_options]
                constraint = "|".join(constraint_parts)
            else:
                # Single memory option (e.g., "32GB")
                constraint = f"GPUMEM{gpu_spec}"
            gpu_line = f"#SBATCH --gpus=1\n#SBATCH --constraint=\"{constraint}\""
        else:
            # Specific GPU model (T4, V100, A100, H100)
            gpu_line = f"#SBATCH --gpus={gpu_spec}:1"
        
        # Use memory specification directly (SLURM accepts human-readable formats)
        memory_spec = resources["memory"]

        # Process additional SLURM options
        additional_sbatch_lines = ""
        all_slurm_options = {}

        # Include SLURM options from resources() method
        if "slurm_options" in resources:
            all_slurm_options.update(resources["slurm_options"])

        # Include SLURM options from slurm() method call (these override resources() options)
        if slurm_options:
            all_slurm_options.update(slurm_options)

        if all_slurm_options:
            sbatch_lines = []
            for param, value in all_slurm_options.items():
                # Handle special case for cpus -> cpus-per-task
                if param == "cpus":
                    slurm_param = "cpus-per-task"
                else:
                    # Convert underscores to hyphens for SLURM parameter names
                    slurm_param = param.replace('_', '-')
                sbatch_lines.append(f"#SBATCH --{slurm_param}={value}")
            additional_sbatch_lines = "\n" + "\n".join(sbatch_lines)

        # Conditional GPU setup
        if gpu_spec is None or gpu_spec == "none":
            gpu_setup = ""
        else:
            gpu_setup = """
# Check if nvidia-smi is available
if ! command -v nvidia-smi &> /dev/null
then
    echo "Could not load GPU correctly: nvidia-smi could not be found"
    exit 1
fi

gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader)
echo "GPU Type: $gpu_type"
"""
        module_load = "module load mamba singularityce"

        # Generate SLURM script
        slurm_content = f"""#!/usr/bin/bash
{gpu_line}
#SBATCH --mem={memory_spec}
#SBATCH --time={resources["time"]}
#SBATCH --output=job.out
#SBATCH --begin=now+0hour{additional_sbatch_lines}
{email_line}

# Make all files group-writable by default
umask 002
{gpu_setup}
{module_load}

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
        
        print(f"Slurm saved to: {slurm_path}")
        print("="*30+"Job"+"="*30)
        print(f"{self.pipeline_name}: {self.job_name} ({self.folder_manager.job_id})")
        print("="*30+"Slurm Script"+"="*30)
        # Print line by line to ensure proper formatting
        for line in slurm_content.split('\n'):
            if line != "":
                print(line)
        print("="*30+"SBATCH"+"="*30)
        output_path = os.path.join(
            self.folders["runtime"],
            f"slurm.out"
        )
        # Format job name as "<pipeline>: <job> (<N>)"
        formatted_job_name = f"{self.pipeline_name}: {self.job_name} ({self.folder_manager.job_id})"
        print(f"sbatch --job-name=\"{formatted_job_name}\" --output {output_path} {slurm_path}")
    
    def resources(self, gpu: str = None, memory: str = None, time: str = None, **slurm_options):
        """
        Configure computational resources for the pipeline.

        Args:
            gpu: GPU specification. Options:
                - Specific models: "T4", "V100", "A100", "H100"
                - Memory-based: "32GB", "80GB", "32GB|80GB" (uses --constraint)
                - Generic: "gpu" (any available GPU)
                - High-memory: "high-memory" (equivalent to "32GB|80GB")
                - Exclusion: "!T4", "!V100", etc. (exclude specific model)
                - "none": No GPU required
            memory: RAM allocation (e.g., "16GB", "32GB", "64GB")
            time: Wall time limit in HH:MM:SS format (e.g., "24:00:00", "48:00:00")
            **slurm_options: Additional SLURM parameters. Examples:
                - cpus=32 → #SBATCH --cpus-per-task=32
                - nodes=1 → #SBATCH --nodes=1
                - ntasks_per_node=1 → #SBATCH --ntasks-per-node=1
                - partition="gpu" → #SBATCH --partition=gpu

        Examples:
            # Specific GPU model
            pipeline.resources(gpu="A100", memory="32GB", time="24:00:00")

            # Memory-based GPU selection (uses SLURM constraints)
            pipeline.resources(gpu="80GB", memory="16GB", time="12:00:00")

            # Multiple memory options
            pipeline.resources(gpu="32GB|80GB", memory="32GB", time="48:00:00")

            # High-memory GPUs
            pipeline.resources(gpu="high-memory", memory="64GB", time="24:00:00")

            # Exclude T4 GPUs (get V100, A100, or H100)
            pipeline.resources(gpu="!T4", memory="32GB", time="24:00:00")

            # CPU-only with multiple cores and nodes
            pipeline.resources(gpu="none", memory="128GB", time="24:00:00", cpus=32, nodes=1)

            # GPU with specific CPU allocation
            pipeline.resources(gpu="V100", memory="32GB", time="12:00:00", cpus=16)
        """
        if gpu:
            self.global_resources["gpu"] = gpu
        if memory:
            self.global_resources["memory"] = memory
        if time:
            self.global_resources["time"] = time

        # Store additional SLURM options for later use in slurm() method
        if slurm_options:
            if "slurm_options" not in self.global_resources:
                self.global_resources["slurm_options"] = {}
            self.global_resources["slurm_options"].update(slurm_options)
    
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
        tool_outputs_dir = os.path.join(self.folders["output"], "ToolOutputs")
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
                
                # Generate filename with 3-digit execution order
                if hasattr(tool, 'suffix') and tool.suffix:
                    output_filename = f"{i:03d}_{tool.TOOL_NAME}_{tool.suffix}.json"
                else:
                    output_filename = f"{i:03d}_{tool.TOOL_NAME}.json"
                output_path = os.path.join(tool_outputs_dir, output_filename)
                
                # Save to JSON file with custom serialization
                with open(output_path, 'w') as f:
                    json.dump(tool_metadata, f, indent=2, default=self._json_serializer)
                
                exported_count += 1
                
            except Exception as e:
                print(f"Warning: Could not export output metadata for {tool.TOOL_NAME}: {e}")
                continue
        
        if exported_count > 0:
            print(f"Exported {exported_count} tool output metadata files to: {tool_outputs_dir}")
        else:
            print("Warning: No tool output metadata could be exported")
    
    def _json_serializer(self, obj):
        """Custom JSON serializer for ToolOutput and other non-serializable objects."""
        if hasattr(obj, 'to_dict'):
            return obj.to_dict()
        elif hasattr(obj, '__dict__'):
            return obj.__dict__
        else:
            return str(obj)

    def _save_original_pipeline_script(self):
        """
        Save the original Python pipeline script to the runtime folder.

        This copies the calling script (e.g., example_pipeline.py) to the runtime
        folder for reproducibility and debugging purposes.
        """
        try:
            # Get the calling script path using inspect
            # We need to go up the call stack to find the original script
            current_frame = inspect.currentframe()
            calling_frame = current_frame

            # Walk up the stack to find the first frame outside this module
            while calling_frame:
                calling_frame = calling_frame.f_back
                if calling_frame is None:
                    break

                frame_filename = calling_frame.f_code.co_filename

                # Skip frames from this module and built-in modules
                if (not frame_filename.endswith('pipeline.py') and
                    not frame_filename.startswith('<') and
                    frame_filename.endswith('.py')):

                    # Found the original calling script
                    original_script_path = os.path.abspath(frame_filename)

                    if os.path.exists(original_script_path):
                        # Get just the filename
                        script_filename = os.path.basename(original_script_path)

                        # Copy to runtime folder
                        runtime_script_path = os.path.join(self.folders["runtime"], script_filename)
                        shutil.copy2(original_script_path, runtime_script_path)

                        print(f"Original pipeline script copied to: {runtime_script_path}")
                        return runtime_script_path

            # If we couldn't find the calling script automatically, try to find it from sys.argv
            import sys
            if len(sys.argv) > 0 and sys.argv[0].endswith('.py'):
                original_script_path = os.path.abspath(sys.argv[0])
                if os.path.exists(original_script_path):
                    script_filename = os.path.basename(original_script_path)
                    runtime_script_path = os.path.join(self.folders["runtime"], script_filename)
                    shutil.copy2(original_script_path, runtime_script_path)
                    print(f"Original pipeline script copied to: {runtime_script_path}")
                    return runtime_script_path

            print("Warning: Could not automatically detect original pipeline script")

        except Exception as e:
            print(f"Warning: Could not save original pipeline script: {e}")

        return None


# Module-level convenience functions for use within Pipeline context
def Resources(**kwargs):
    """
    Set global resources for all tools in the active pipeline.

    Must be called within a Pipeline context manager.

    Args:
        gpu: GPU memory (e.g., "32GB", "16GB", or None for no GPU)
        memory: System RAM (e.g., "16GB", "32GB")
        time: Wall time limit (e.g., "24:00:00", "2:00:00")

    Example:
        with Pipeline("Test", "Job", "Description"):
            Resources(gpu="32GB", memory="16GB", time="24:00:00")
            tool1 = SomeTool(...)

    Raises:
        RuntimeError: If called outside a Pipeline context
    """
    pipeline = Pipeline.get_active_pipeline()
    if pipeline is None:
        raise RuntimeError(
            "Resources() must be called within a Pipeline context. "
            "Use: with Pipeline(...): Resources(...)"
        )
    pipeline.resources(**kwargs)


def Suffix(suffix: str = ""):
    """
    Set suffix for subsequent tool folder names in the active pipeline.

    Must be called within a Pipeline context manager.

    Args:
        suffix: Suffix to append to tool folder names (e.g., "001", "batch1")
               Empty string clears the suffix.

    Example:
        with Pipeline("Test", "Job", "Description"):
            Suffix("001")
            tool1 = SomeTool(...)  # Creates folder like: 001_SomeTool_001

            Suffix("002")
            tool2 = OtherTool(...)  # Creates folder like: 002_OtherTool_002

    Raises:
        RuntimeError: If called outside a Pipeline context
    """
    pipeline = Pipeline.get_active_pipeline()
    if pipeline is None:
        raise RuntimeError(
            "Suffix() must be called within a Pipeline context. "
            "Use: with Pipeline(...): Suffix('...')"
        )
    pipeline.set_suffix(suffix)


def Save():
    """
    Save pipeline configuration without submitting to SLURM.

    Must be called within a Pipeline context manager. Useful if you want
    to save the pipeline but not submit it immediately.

    Example:
        with Pipeline("Test", "Job", "Description"):
            tool1 = SomeTool(...)
            Save()  # Saves but doesn't submit

    Note: Calling Save() prevents auto-submission on context exit.

    Raises:
        RuntimeError: If called outside a Pipeline context
    """
    pipeline = Pipeline.get_active_pipeline()
    if pipeline is None:
        raise RuntimeError(
            "Save() must be called within a Pipeline context. "
            "Use: with Pipeline(...): Save()"
        )
    pipeline.save()
    # Mark that save was explicitly called to prevent auto-submit
    pipeline._explicit_save_called = True