# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Core Pipeline class for automated workflow construction.

Manages tool sequencing, environment switching, dependency resolution,
and script generation for protein modeling pipelines.
"""

import os
import sys
import json
import getpass
import inspect
import shutil
import subprocess
import contextvars
from typing import Dict, List, Any, Optional, Union
from collections import defaultdict
from datetime import datetime

from .folders import FolderManager
from .config_manager import ConfigManager
try:
    from .base_config import BaseConfig, ToolOutput
    from .combinatorics import Bundle, Each
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput
    from combinatorics import Bundle, Each

# Module-level context variable to track active pipeline for auto-registration
_active_pipeline: contextvars.ContextVar[Optional['Pipeline']] = contextvars.ContextVar('_active_pipeline', default=None)

class Pipeline:
    """
    Main pipeline orchestration class.
    
    Manages tools, dependencies, environment switching, and script generation
    for automated protein modeling workflows.
    """
    
    def __init__(self, project: str, job: str, description: str="Description missing", on_the_fly: Optional[bool]=None, local_output: Optional[bool]=None):
        """
        Initialize a new pipeline instance.

        Args:
            project: Name of the folder (used for output folders)
            job: Name of the specific job (a unique numeric id NNN will be appended to it) (used for output folders)
            description: Optional description of what this job does
            on_the_fly: If True, each tool's script is executed immediately when added.
                        Useful for interactive use in Jupyter notebooks or running locally
                        with plain Python. Skips SLURM submission on exit.
                        If None (default), auto-detects: True when running in a Jupyter
                        notebook (.ipynb), False otherwise.
            local_output: If True, write output to ./tests/ (current working
                        directory) instead of the config-defined path.
                        If None (default), follows on_the_fly.
        """
        if ' ' in job: job=job.replace(' ','_') #It will create issues at runtime otherwise

        self.project = project
        self.job = job
        self.description = description

        if on_the_fly is None:
            on_the_fly = self._detect_notebook()
        self.on_the_fly = on_the_fly

        if local_output is None:
            local_output = on_the_fly
        self.local_output = local_output

        self.folder_manager = FolderManager(project, job, local_output=local_output)
        self.folders = self.folder_manager.get_folders()
        # Include resolved container paths as container:<ToolName> keys
        for tool_name, container_path in self.folder_manager.get_containers().items():
            self.folders[f"container:{tool_name}"] = container_path
        
        # Tool management
        self.tools = []
        self.tool_outputs = []
        self.execution_order = 0

        # Suffix for tool folder naming
        self.current_suffix = ""

        # Script generation
        self.pipeline_script = ""
        self.slurm_script = ""
        self.scripts_generated = False

        # Batch and resource management
        self.batch_resources = []  # List of resource dicts, one per batch
        self.batch_start_indices = []  # Tool indices where each batch starts
        self.current_batch = -1  # Current batch index (-1 means no Resources() called yet)

        # External job dependencies
        self.external_dependencies = []  # List of external SLURM job IDs to depend on

        # Context manager state
        self._explicit_save_called = False  # Track if Save() was explicitly called

        # In on_the_fly mode, activate the pipeline context immediately
        # so tools can auto-register across notebook cells without requiring
        # the `with` statement. Replace any previously active pipeline.
        if self.on_the_fly:
            _active_pipeline.set(self)

    @staticmethod
    def _detect_notebook() -> bool:
        """
        Detect whether the pipeline is being run from a Jupyter notebook.

        Checks the call stack for IPython's interactiveshell module,
        which is present when code executes inside a Jupyter notebook.

        Returns:
            True if running inside a Jupyter notebook, False otherwise.
        """
        try:
            from IPython import get_ipython
            shell = get_ipython()
            if shell is not None:
                if (shell.__class__.__name__ == 'ZMQInteractiveShell'
                        or shell.__class__.__module__ == 'google.colab._shell'):
                    return True
        except (ImportError, NameError):
            pass
        return False

    def __enter__(self):
        """
        Enter context manager - set this pipeline as the active context.

        Returns:
            self (Pipeline instance)

        Raises:
            RuntimeError: If a different pipeline context is already active
                         (nested contexts not allowed, unless on_the_fly already
                         set this pipeline as active)
        """
        current_pipeline = _active_pipeline.get()
        if current_pipeline is not None and current_pipeline is not self:
            raise RuntimeError(
                f"Cannot nest Pipeline contexts. "
                f"Pipeline '{current_pipeline.project}' is already active."
            )
        _active_pipeline.set(self)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context manager - auto-submit to SLURM or keep context alive.

        In normal mode: clears the active pipeline and auto-submits to SLURM.
        In on_the_fly mode: keeps the pipeline active so subsequent notebook
        cells can continue adding tools.

        Args:
            exc_type: Exception type if an exception occurred
            exc_val: Exception value if an exception occurred
            exc_tb: Exception traceback if an exception occurred

        Returns:
            False (don't suppress exceptions)
        """
        if self.on_the_fly:
            # Keep the pipeline active for subsequent notebook cells
            return False

        try:
            # Only auto-submit if no exception and Save() wasn't explicitly called
            if exc_type is None and not self._explicit_save_called:
                self.slurm()
        finally:
            # Clear the active pipeline context
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
    
    def add(self, tool_config: BaseConfig, **kwargs):
        """
        Add a tool to the pipeline.

        Args:
            tool_config: Configured tool instance (e.g., RFdiffusion())
            **kwargs: Additional pipeline-specific options

        Returns:
            ToolAwareOutput object with standardized output interface and .tool access
        """
        # Ensure Resources() has been called before adding tools
        if self.current_batch == -1:
            raise RuntimeError(
                "Resources() must be called before adding tools. "
                "Use: with Pipeline(...): Resources(...); tool = Tool(...)"
            )

        # Merge current batch resources with tool resources
        current_resources = self.batch_resources[self.current_batch]
        for key, value in current_resources.items():
            if value is not None:
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
        # Ensure Resources() has been called before adding tools
        # In on_the_fly mode, auto-initialize resources if not set
        if self.current_batch == -1:
            if self.on_the_fly:
                self.resources()
            else:
                raise RuntimeError(
                    "Resources() must be called before adding tools. "
                    "Use: with Pipeline(...): Resources(...); tool = Tool(...)"
                )

        # Merge current batch resources with tool resources
        current_resources = self.batch_resources[self.current_batch]
        for key, value in current_resources.items():
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

        # On-the-fly execution: generate and run the tool's script immediately
        if self.on_the_fly:
            self._execute_tool_on_the_fly(tool_config)

        return tool_output

    def _execute_tool_on_the_fly(self, tool_config: BaseConfig):
        """
        Generate and immediately execute a tool's bash script.

        Called during on_the_fly mode after a tool is registered. Streams
        stdout/stderr in real-time so output is visible in notebooks and terminals.

        Args:
            tool_config: The tool to execute

        Raises:
            RuntimeError: If the tool's script exits with a non-zero code
        """
        step_idx = self.execution_order
        if hasattr(tool_config, 'suffix') and tool_config.suffix:
            tool_script_name = f"{step_idx:03d}_{tool_config.TOOL_NAME}_{tool_config.suffix}.sh"
            log_name = f"{step_idx:03d}_{tool_config.TOOL_NAME}_{tool_config.suffix}.log"
        else:
            tool_script_name = f"{step_idx:03d}_{tool_config.TOOL_NAME}.sh"
            log_name = f"{step_idx:03d}_{tool_config.TOOL_NAME}.log"

        tool_script_path = os.path.join(self.folders["runtime"], tool_script_name)
        log_path = os.path.join(self.folders["logs"], log_name)

        # Generate the tool script
        script_content = tool_config.generate_script(tool_script_path)
        with open(tool_script_path, 'w') as f:
            f.write(script_content)
        os.chmod(tool_script_path, 0o755)

        print(f"\n{'='*60}")
        print(f"Running {tool_config.TOOL_NAME} (step {step_idx})")
        print(f"{'='*60}")

        # Execute the script, streaming output in real-time
        with open(log_path, 'w') as log_file:
            process = subprocess.Popen(
                ['bash', '-e', tool_script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1  # Line-buffered
            )
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log_file.write(line)
                log_file.flush()
            process.wait()

        if process.returncode != 0:
            raise RuntimeError(
                f"{tool_config.TOOL_NAME} failed with exit code {process.returncode}. "
                f"See log: {log_path}"
            )

        print(f"{'='*60}")
        print(f"{tool_config.TOOL_NAME} completed successfully")
        print(f"{'='*60}\n")

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
            f'echo "Pipeline: {self.project}"',
            f'echo "Job name: {self.job}"',
            f'echo "Description: {self.description}"',
            'echo'
        ]
        
        # Add tool configurations
        for i, tool in enumerate(self.tools, 1):
            config_lines.extend([
                f'echo "Step {i}: {tool.TOOL_NAME}"',
                f'echo "  Environments: {tool.environments}"'
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
            f"# Generated pipeline: {self.project}",
            f"# Tools: {', '.join([t.TOOL_NAME for t in self.tools])}",
            "",
            "set -e  # Exit on any error",
            "umask 002  # Make all files group-writable by default",
            "",
            "# Initialize environment manager",
            ConfigManager().get_shell_hook_command(),
            "",
            "echo Configuration",
            f"{config_script} | tee {os.path.join(self.folders['output'], f'{self.project}_config.txt')}",
            "echo"
        ]

        # Add each tool execution (each tool handles its own environment activation)
        for i, tool in enumerate(self.tools, 1):
            # Generate tool-specific script
            # Include suffix in script name if present
            if hasattr(tool, 'suffix') and tool.suffix:
                tool_script_path = os.path.join(self.folders["runtime"], f"{i:03d}_{tool.TOOL_NAME}_{tool.suffix}.sh")
            else:
                tool_script_path = os.path.join(self.folders["runtime"], f"{i:03d}_{tool.TOOL_NAME}.sh")
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
            # Include suffix in log file name if present
            if hasattr(tool, 'suffix') and tool.suffix:
                log_file = os.path.join(self.folders["logs"], f"{i:03d}_{tool.TOOL_NAME}_{tool.suffix}.log")
            else:
                log_file = os.path.join(self.folders["logs"], f"{i:03d}_{tool.TOOL_NAME}.log")
            script_lines.extend([
                f"echo {tool.TOOL_NAME}",
                f"{tool_script_path} 2>&1 | tee {log_file}",
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
    
    def slurm(self, email: str = "auto"):
        """
        Generate SLURM job submission script(s).

        If multiple batches exist (multiple Resources() calls), generates:
        - Separate SLURM script for each batch
        - Submission chain script that submits jobs with dependencies

        Args:
            email: Email for job notifications. Options:
                - "auto": Automatically detect email from current username (default)
                - "": No email notifications
                - specific email: Use provided email address

        Returns:
            None (prints submission instructions)
        """
        if not self.scripts_generated:
            self.save()

        # Validate that Resources() was called
        if self.current_batch == -1:
            raise RuntimeError("Resources() must be called before generating SLURM scripts")

        # Handle auto email detection
        if email == "auto":
            current_user = getpass.getuser()
            emails = ConfigManager().get_emails()
            if current_user in emails:
                email = emails[current_user]
                print(f"Auto-detected email: {email} (user: {current_user})")
            else:
                print(f"Warning: No email mapping found for user '{current_user}'. Disabling email notifications.")
                print(f"Add your username to the 'cluster.emails' section in config.yaml.")
                email = ""

        email_line = "" if email == "" else f"""
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={email}"""

        num_batches = len(self.batch_resources)

        # Determine if single batch or multiple batches
        if num_batches == 1:
            # Single batch - generate single SLURM script
            self._generate_single_batch_slurm(email_line)
        else:
            # Multiple batches - generate chained SLURM scripts
            self._generate_multi_batch_slurm(email_line, num_batches)

    def _generate_gpu_line(self, gpu_spec):
        """Generate SBATCH GPU constraint line based on GPU specification."""
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return ""
        elif gpu_spec == "high-memory":
            return "#SBATCH --gpus=1\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB|GPUMEM96GB\""
        elif gpu_spec == "gpu" or gpu_spec == "any":
            return "#SBATCH --gpus=1"
        elif gpu_spec.startswith("!"):
            excluded_model = gpu_spec[1:]
            if excluded_model.upper() == "L4":
                return "#SBATCH --gpus=1\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB|GPUMEM96GB\""
            else:
                return f"#SBATCH --gpus=1\n#SBATCH --constraint=\"~GPU{excluded_model}\""
        elif gpu_spec in ["24GB", "32GB", "80GB", "96GB"] or "|" in gpu_spec:
            if "|" in gpu_spec:
                memory_options = gpu_spec.split("|")
                constraint_parts = [f"GPUMEM{mem}" for mem in memory_options]
                constraint = "|".join(constraint_parts)
            else:
                constraint = f"GPUMEM{gpu_spec}"
            return f"#SBATCH --gpus=1\n#SBATCH --constraint=\"{constraint}\""
        else:
            return f"#SBATCH --gpus={gpu_spec}:1"

    def _generate_gpu_setup(self, gpu_spec):
        """Generate GPU setup script based on GPU specification."""
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return ""
        return """
# Display GPU information
gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader 2>/dev/null || echo "Unknown")
echo "GPU Type: $gpu_type"
"""

    def _generate_additional_sbatch_lines(self, slurm_options):
        """Generate additional SBATCH lines from slurm_options dict."""
        if not slurm_options:
            return ""
        sbatch_lines = []
        for param, value in slurm_options.items():
            if param == "cpus":
                slurm_param = "cpus-per-task"
            else:
                slurm_param = param.replace('_', '-')
            sbatch_lines.append(f"#SBATCH --{slurm_param}={value}")
        return "\n" + "\n".join(sbatch_lines)

    def _generate_dependency_line(self, job_ids: List[str]) -> str:
        """Generate SBATCH dependency line from list of job IDs."""
        if not job_ids:
            return ""
        return f"\n#SBATCH --dependency=afterok:{':'.join(job_ids)}"

    def _generate_single_batch_slurm(self, email_line):
        """Generate single SLURM script for single batch pipeline."""
        resources = self.batch_resources[0]
        gpu_line = self._generate_gpu_line(resources["gpu"])
        gpu_setup = self._generate_gpu_setup(resources["gpu"])
        additional_sbatch_lines = self._generate_additional_sbatch_lines(resources.get("slurm_options", {}))
        dependency_line = self._generate_dependency_line(self.external_dependencies)
        module_load = ConfigManager().get_module_load_line()

        slurm_content = f"""#!/usr/bin/bash
{gpu_line}
#SBATCH --mem={resources["memory"]}
#SBATCH --time={resources["time"]}
#SBATCH --output=job.out
#SBATCH --begin=now+0hour{additional_sbatch_lines}{dependency_line}
{email_line}

# Make all files group-writable by default
umask 002
{gpu_setup}
{module_load}

# Execute pipeline
{self.pipeline_script}
"""
        slurm_path = os.path.join(self.folders["runtime"], "slurm.sh")
        with open(slurm_path, 'w') as f:
            f.write(slurm_content)
        os.chmod(slurm_path, 0o755)
        self.slurm_script = slurm_path

        print(f"Slurm saved to: {slurm_path}")
        print("="*30+"Job"+"="*30)
        print(f"{self.project}: {self.job} ({self.folder_manager.job_id})")
        print("="*30+"Slurm Script"+"="*30)
        for line in slurm_content.split('\n'):
            if line != "":
                print(line)
        print("="*30+"SBATCH"+"="*30)
        output_path = os.path.join(self.folders["runtime"], "slurm.out")
        formatted_job_name = f"{self.project}: {self.job} ({self.folder_manager.job_id})"
        print(f"sbatch --job-name=\"{formatted_job_name}\" --output {output_path} {slurm_path}")

    def _generate_multi_batch_slurm(self, email_line, num_batches):
        """Generate multiple SLURM scripts for multi-batch pipeline with <JOBID> placeholders."""
        module_load = ConfigManager().get_module_load_line()

        # Determine batch ranges
        batch_ranges = []
        for batch_idx in range(num_batches):
            start_idx = self.batch_start_indices[batch_idx]
            end_idx = self.batch_start_indices[batch_idx + 1] if batch_idx + 1 < num_batches else len(self.tools)
            batch_ranges.append((start_idx, end_idx))

        # Generate pipeline scripts for each batch
        for batch_idx, (start_idx, end_idx) in enumerate(batch_ranges):
            batch_script_content = self._generate_batch_script(batch_idx, start_idx, end_idx)
            batch_script_path = os.path.join(self.folders["runtime"], f"pipeline_batch{batch_idx + 1}.sh")
            with open(batch_script_path, 'w') as f:
                f.write(batch_script_content)
            os.chmod(batch_script_path, 0o755)

        # Generate SLURM scripts for each batch
        for batch_idx in range(num_batches):
            resources = self.batch_resources[batch_idx]
            gpu_line = self._generate_gpu_line(resources["gpu"])
            gpu_setup = self._generate_gpu_setup(resources["gpu"])
            additional_sbatch_lines = self._generate_additional_sbatch_lines(resources.get("slurm_options", {}))

            # Add dependency line for batch 2+ (internal) and batch 1 (external)
            if batch_idx > 0:
                dependency_line = "\n#SBATCH --dependency=afterok:<JOBID>"
            elif self.external_dependencies:
                dependency_line = self._generate_dependency_line(self.external_dependencies)
            else:
                dependency_line = ""

            batch_slurm_content = f"""#!/usr/bin/bash
{gpu_line}
#SBATCH --mem={resources["memory"]}
#SBATCH --time={resources["time"]}
#SBATCH --output=job_batch{batch_idx + 1}.out
#SBATCH --begin=now+0hour{additional_sbatch_lines}{dependency_line}
{email_line}

# Make all files group-writable by default
umask 002
{gpu_setup}
{module_load}

# Execute pipeline batch {batch_idx + 1}
{os.path.join(self.folders["runtime"], f"pipeline_batch{batch_idx + 1}.sh")}
"""
            batch_slurm_path = os.path.join(self.folders["runtime"], f"slurm_batch{batch_idx + 1}.sh")
            with open(batch_slurm_path, 'w') as f:
                f.write(batch_slurm_content)
            os.chmod(batch_slurm_path, 0o755)

        self.slurm_script = os.path.join(self.folders["runtime"], "slurm_batch1.sh")

        # Print summary
        print(f"Generated {num_batches} batch job(s) with dependencies")
        print("="*30+"Job Batches"+"="*30)
        print(f"{self.project}: {self.job} ({self.folder_manager.job_id})")
        for batch_idx in range(num_batches):
            start_idx, end_idx = batch_ranges[batch_idx]
            batch_tools = self.tools[start_idx:end_idx]
            res = self.batch_resources[batch_idx]
            print(f"Batch {batch_idx + 1}: Tools {start_idx + 1}-{end_idx} ({', '.join(t.TOOL_NAME for t in batch_tools)}) | GPU: {res['gpu']}, Mem: {res['memory']}, Time: {res['time']}")
        print("="*30+"Manual Submission"+"="*30)
        print("# Submit batches manually (replace <JOBID> with previous job ID):")
        for batch_idx in range(num_batches):
            print(f"sbatch {os.path.join(self.folders['runtime'], f'slurm_batch{batch_idx + 1}.sh')}")
        print("="*30+"Automatic Submission"+"="*30)
        print(f"# Use submit script to automatically chain submissions")

    def _generate_batch_script(self, batch_idx: int, start_tool_idx: int, end_tool_idx: int) -> str:
        """Generate pipeline script for a specific batch of tools."""
        batch_tools = self.tools[start_tool_idx:end_tool_idx]

        # Build configuration display
        config_lines = [
            "#!/bin/bash",
            f'echo "Pipeline: {self.project}"',
            f'echo "Job: {self.job}"',
            f'echo "Batch: {batch_idx + 1}"',
            f'echo "Description: {self.description}"',
            'echo',
            f'echo "Tools in this batch: {len(batch_tools)}"',
            'echo'
        ]

        for i, tool in enumerate(batch_tools, start=start_tool_idx + 1):
            config_lines.extend([
                f'echo "Step {i}: {tool.TOOL_NAME}"',
                f'echo "  Environments: {tool.environments}"'
            ])
            tool_config = tool.get_config_display()
            for config_line in tool_config:
                config_lines.append(f'echo "  {config_line}"')
            config_lines.append('echo')

        # Create batch-specific config script
        config_script = os.path.join(self.folders["runtime"], f"config_batch{batch_idx + 1}.sh")
        with open(config_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(config_lines))
        os.chmod(config_script, 0o755)

        # Build batch pipeline script
        script_lines = [
            "#!/bin/bash",
            f"# Generated pipeline batch {batch_idx + 1}: {self.project}",
            f"# Tools: {', '.join([t.TOOL_NAME for t in batch_tools])}",
            "",
            "set -e  # Exit on any error",
            "umask 002  # Make all files group-writable by default",
            "",
            "# Initialize environment manager",
            ConfigManager().get_shell_hook_command(),
            "",
            "echo Configuration",
            f"{config_script} | tee -a {os.path.join(self.folders['output'], f'{self.project}_config.txt')}",
            "echo"
        ]

        # Add each tool execution (each tool handles its own environment activation)
        for i, tool in enumerate(batch_tools, start=start_tool_idx + 1):
            # Tool script path
            if hasattr(tool, 'suffix') and tool.suffix:
                tool_script_path = os.path.join(self.folders["runtime"], f"{i:03d}_{tool.TOOL_NAME}_{tool.suffix}.sh")
                log_file = os.path.join(self.folders["logs"], f"{i:03d}_{tool.TOOL_NAME}_{tool.suffix}.log")
            else:
                tool_script_path = os.path.join(self.folders["runtime"], f"{i:03d}_{tool.TOOL_NAME}.sh")
                log_file = os.path.join(self.folders["logs"], f"{i:03d}_{tool.TOOL_NAME}.log")

            script_lines.extend([f"echo {tool.TOOL_NAME}", f"{tool_script_path} 2>&1 | tee {log_file}", "echo"])

        # Final steps
        script_lines.extend(["echo", f"echo Batch {batch_idx + 1} done"])
        if batch_idx == len(self.batch_resources) - 1:  # Last batch
            script_lines.extend(["echo", "echo All jobs complete", f"echo Results in:", f"echo {self.folders['output']}"])

        return "\n".join(script_lines)

    def set_dependencies(self, job_ids: Union[str, List[str]]):
        """
        Set external SLURM job dependencies.

        The first batch of this pipeline will wait for these jobs to complete
        successfully before starting.

        Args:
            job_ids: Single job ID string or list of job ID strings
                    e.g., "12345678" or ["12345678", "12345679"]
        """
        if isinstance(job_ids, str):
            job_ids = [job_ids]
        self.external_dependencies.extend(job_ids)

    def resources(self, gpu: str = None, memory: str = None, time: str = None, **slurm_options):
        """
        Configure computational resources and start a new batch.

        First call: Initializes batch 0 with specified resources.
        Subsequent calls: Start new batch, inheriting unspecified resources from previous batch.

        Args:
            gpu: GPU specification. Options:
                - Specific models: "L4", "V100", "A100", "H100", "H200"
                - Memory-based: "24GB", "32GB", "80GB", "96GB", "32GB|80GB|96GB" (uses --constraint)
                - Generic: "gpu" (any available GPU)
                - High-memory: "high-memory" (equivalent to "32GB|80GB|96GB")
                - Exclusion: "!L4", "!V100", etc. (exclude specific model)
                - None: Inherits from previous batch (or no GPU if first batch)
            memory: RAM allocation (e.g., "16GB", "32GB", "64GB")
                - None: Inherits from previous batch (or "15GB" if first batch)
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
        # Start new batch
        self.current_batch += 1
        self.batch_start_indices.append(len(self.tools))

        # Determine resources for this batch
        if self.current_batch == 0:
            # First batch: use provided values or defaults
            batch_res = {
                "gpu": gpu if gpu is not None else None,
                "memory": memory if memory is not None else "15GB",
                "time": time if time is not None else "24:00:00",
                "slurm_options": slurm_options if slurm_options else {}
            }
        else:
            # Subsequent batches: inherit from previous if not specified
            prev_res = self.batch_resources[self.current_batch - 1]
            batch_res = {
                "gpu": gpu if gpu is not None else prev_res["gpu"],
                "memory": memory if memory is not None else prev_res["memory"],
                "time": time if time is not None else prev_res["time"],
                "slurm_options": slurm_options if slurm_options else prev_res.get("slurm_options", {})
            }

        self.batch_resources.append(batch_res)

        if self.current_batch > 0:
            print(f"Starting batch {self.current_batch + 1} with resources: gpu={batch_res['gpu']}, memory={batch_res['memory']}, time={batch_res['time']}")
    
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
            f"Pipeline: {self.project}",
            f"Tools: {len(self.tools)}",
            ""
        ]

        for i, tool in enumerate(self.tools, 1):
            summary.append(
                f"{i}. {tool.TOOL_NAME} ({tool.environments}) -> {tool.job_name}"
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

                # Convert DataStream objects to dictionaries for JSON serialization
                from .datastream import DataStream
                from .base_config import TableInfo

                for key in ['structures', 'sequences', 'compounds', 'msas']:
                    if key in output_structure and isinstance(output_structure[key], DataStream):
                        output_structure[key] = output_structure[key].to_dict()

                # Convert TableInfo objects to dictionaries for JSON serialization
                if 'tables' in output_structure:
                    serializable_tables = {}
                    for name, table_info in output_structure['tables'].items():
                        if isinstance(table_info, TableInfo):
                            serializable_tables[name] = table_info.to_dict()
                        else:
                            serializable_tables[name] = table_info
                    output_structure['tables'] = serializable_tables
                
                # Build complete tool metadata
                tool_metadata = {
                    "tool_name": tool.TOOL_NAME,
                    "tool_class": tool.__class__.__name__,
                    "job_name": tool.job_name,
                    "execution_order": i,
                    "environments": tool.environments,
                    "output_structure": output_structure,
                    "configuration": {
                        "tool_parameters": tool.to_dict(),
                        "resources": tool.resources.copy(),
                        "pipeline_context": {
                            "pipeline_name": self.project,
                            "pipeline_job_name": self.job,
                            "pipeline_description": self.description
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
        gpu: GPU memory ("T4", "L4", "V100", "A100", "H100", "H200", "16GB", "24GB", "32GB", "80GB", "96GB", "32GB|80GB|96GB", "!L4", "gpu", "high-memory", None))
        memory: System RAM (e.g., "16GB", "32GB")
        time: Wall time limit (e.g., "24:00:00", "2:00:00", "1-00:00:00")

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


def Dependencies(job_ids: Union[str, List[str]]):
    """
    Set external SLURM job dependencies for the active pipeline.

    The first batch of this pipeline will wait for the specified jobs to
    complete successfully before starting.

    Must be called within a Pipeline context manager.

    Args:
        job_ids: Single job ID string or list of job ID strings
                e.g., "12345678" or ["12345678", "12345679"]

    Example:
        with Pipeline("Test", "Job", "Description"):
            Dependencies("12345678")  # Wait for job 12345678
            Resources(gpu="A100", memory="16GB", time="24:00:00")
            tool1 = SomeTool(...)

        # Or with multiple dependencies:
        with Pipeline("Test", "Job", "Description"):
            Dependencies(["12345678", "12345679"])
            Resources(gpu="A100", memory="16GB", time="24:00:00")
            tool1 = SomeTool(...)

    Raises:
        RuntimeError: If called outside a Pipeline context
    """
    pipeline = Pipeline.get_active_pipeline()
    if pipeline is None:
        raise RuntimeError(
            "Dependencies() must be called within a Pipeline context. "
            "Use: with Pipeline(...): Dependencies(...)"
        )
    pipeline.set_dependencies(job_ids)