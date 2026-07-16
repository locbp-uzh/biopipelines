# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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
from ._layout import INTERNAL_FOLDER
from .schedulers import get_backend, BATCH_SCHEDULERS
try:
    from .base_config import BaseConfig, ToolOutput
    from .combinatorics import Bundle, Each
    from .entities import *
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput
    from combinatorics import Bundle, Each
    from entities import *

import re as _re

_IDENTIFIER_RE = _re.compile(r"^[A-Za-z0-9._\-]+$")


def _validate_identifier(field: str, value: str) -> None:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{field!r} must be a non-empty string, got {value!r}")
    if not _IDENTIFIER_RE.match(value):
        raise ValueError(
            f"{field!r}={value!r} contains characters outside [A-Za-z0-9._-]. "
            "This value is interpolated into generated shell scripts and file "
            "paths; restrict it to ASCII letters, digits, '.', '_', '-'."
        )
    if value in (".", "..") or value.startswith("-"):
        raise ValueError(
            f"{field!r}={value!r} is not allowed: must not be '.' or '..' "
            "(path traversal) and must not start with '-' (would be parsed "
            "as a flag by sbatch and other tools)."
        )


# Re-export the shared helper so existing pipeline.py callers keep working.
# The canonical definition lives in base_config alongside _validate_freeform_string.
from .base_config import _escape_for_double_quotes  # noqa: F401

# Directive-injection guard; canonical definition lives in schedulers.py.
from .schedulers import _validate_directive_value


# Module-level context variable to track active pipeline for auto-registration
_active_pipeline: contextvars.ContextVar[Optional['Pipeline']] = contextvars.ContextVar('_active_pipeline', default=None)

class Pipeline:
    """
    Main pipeline orchestration class.
    
    Manages tools, dependencies, environment switching, and script generation
    for automated protein modeling workflows.
    """
    
    def __init__(self, project: str, job: str, description: str="Description missing", on_the_fly: Optional[bool]=None, local_output: Optional[bool]=None, config: Optional[str]=None, debug: bool=False, otf: Optional[bool]=None):
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
            otf: Alias for on_the_fly. If set (not None), takes precedence over on_the_fly.
            local_output: If True, write output to ./outputs/ (current working
                        directory) instead of the config-defined path.
                        If None (default), follows on_the_fly.
            config: Config variant name (e.g. "cluster", "colab", or any user-defined
                        site name matching `config.<variant>.yaml`). If None (default),
                        auto-detects: "colab" when running inside Google Colab,
                        otherwise "cluster". Must be set before any folders or
                        environments are resolved.
            debug: If True, save() prints per-tool outputs and the pipeline emits
                        a runtime snapshot under ``<output>/_debug_capture/``
                        (mamba env exports, pip freeze, system info) when
                        ``BIOPIPELINES_DEBUG=1`` is set in the script. Defaults
                        to False.
        """
        if ' ' in job: job=job.replace(' ','_') #It will create issues at runtime otherwise

        _validate_identifier("project", project)
        _validate_identifier("job", job)

        # Apply config variant before any ConfigManager-dependent code runs
        # (FolderManager, environment loading, etc. all read through ConfigManager).
        ConfigManager(variant=config)

        self.project = project
        self.job = job
        self.description = description
        self.debug = bool(debug)

        if otf is not None:
            on_the_fly = otf

        if on_the_fly is None:
            if os.environ.get("BIOPIPELINES_OTF") == "1":
                on_the_fly = True
            else:
                on_the_fly = self._detect_notebook()
        self.on_the_fly = on_the_fly

        if local_output is None:
            # Explicit override wins: single-node container backends (Modal, RunPod,
            # a plain Docker GPU box, ...) repoint biopipelines_output at a persistent
            # mount and must set BIOPIPELINES_LOCAL_OUTPUT=0 so results are NOT
            # diverted to the ephemeral cwd/outputs and lost on teardown.
            env_lo = os.environ.get("BIOPIPELINES_LOCAL_OUTPUT")
            if env_lo is not None:
                local_output = env_lo == "1"
            else:
                # Auto-enable local output for interactive/notebook runs, which usually
                # lack shared storage. EXCEPT on Colab: there the config's
                # biopipelines_output is deliberately repointed at mounted Drive for
                # persistence, and forcing local_output would overwrite it with the
                # ephemeral cwd/outputs (lost on runtime recycle). Let the config win.
                local_output = on_the_fly and ConfigManager().get_scheduler() != "colab"
        self.local_output = local_output

        self.folder_manager = FolderManager(project, job, local_output=local_output)
        self.folders = self.folder_manager.get_folders()
        # Include resolved container paths as container:<ToolName> keys
        for tool_name, container_path in self.folder_manager.get_containers().items():
            self.folders[f"container:{tool_name}"] = container_path
        # Expose bind mounts for container prefix construction
        self.folders["__container_binds__"] = self.folder_manager.get_container_binds()
        
        # Tool management
        self.tools = []
        self.tool_outputs = []
        self.execution_order = 0      # all tools, real run order
        self.public_step_order = 0    # public tools only, drives public folder/script names
        self.internal_order = 0       # internal tools only, drives .internal layout + ".NNN" names

        # Suffix for tool folder naming
        self.current_suffix = ""

        # Stack of public subfolder names pushed by the Folder() context manager.
        self._folder_stack = []

        # Script generation
        self.pipeline_script = ""
        self.job_script = ""
        self.scripts_generated = False

        # Batch and resource management
        self.batch_resources = []  # List of resource dicts, one per batch
        self.batch_start_indices = []  # Tool indices where each batch starts
        self.current_batch = -1  # Current batch index (-1 means no Resources() called yet)
        # Per-batch parent indices for the dependency DAG. batch_parents[i]
        # is the list of batch indices that batch i depends on (afterok).
        # The chain default is [i-1] for i>0 and [] for i=0; the Parallel
        # context manager populates this with multi-parent or shared-parent
        # values to express fan-out / fan-in.
        self.batch_parents = []  # type: List[List[int]]
        # Parallel list of "after" (started, not afterok) parents per batch.
        # Populated only by a Service block: the batch after the block waits
        # for the service daemon to START, not finish (afterok would mean
        # "wait for the long-running server to exit"). Emitted as a separate
        # `after:` term in the same --dependency line.
        self.batch_after_parents = []  # type: List[List[int]]

        # Parallel-block state. Set by `with Parallel():`; consumed by the
        # next Resources() call inside (sibling-with-shared-parent) and the
        # next Resources() call after the block exits (fan-in).
        self._parallel_anchor = None         # type: Optional[int]
        self._parallel_siblings = []         # type: List[int]
        self._parallel_after_anchor = None   # type: Optional[int]
        self._pending_post_parents = None    # type: Optional[List[int]]

        # Service-block state. Set on `with Service():` exit to the daemon's
        # batch index; consumed by the NEXT Resources() call, which records it
        # as an `after:` parent so that batch starts once the daemon is running.
        self._service_anchor = None          # type: Optional[int]
        self._pending_after_parent = None    # type: Optional[int]

        # External job dependencies
        self.external_dependencies = []  # List of external scheduler job IDs to depend on

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
        Exit context manager - generate job scripts or keep context alive.

        In normal mode: clears the active pipeline and generates the batch
        job scripts. In on_the_fly mode: keeps the pipeline active so
        subsequent notebook cells can continue adding tools.

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
            # Only generate scripts if no exception and Save() wasn't explicitly called
            if exc_type is None and not self._explicit_save_called:
                self.generate_job_scripts()
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
    
    def _assign_step_context(self, tool_config: BaseConfig):
        """Bump counters, create the step folder, and set the tool's pipeline context."""
        self.execution_order += 1
        if tool_config.internal:
            self.internal_order += 1
            folder_name = f"{self.internal_order:03d}_{tool_config.TOOL_NAME}"
            tool_output_folder = os.path.join(self.folders["output"], INTERNAL_FOLDER, folder_name)
            os.makedirs(os.path.join(self.folders["runtime"], INTERNAL_FOLDER), exist_ok=True)
            os.makedirs(os.path.join(self.folders["logs"], INTERNAL_FOLDER), exist_ok=True)
            public_step, internal_order = None, self.internal_order
        else:
            self.public_step_order += 1
            if self.current_suffix:
                folder_name = f"{self.public_step_order:03d}_{tool_config.TOOL_NAME}_{self.current_suffix}"
            else:
                folder_name = f"{self.public_step_order:03d}_{tool_config.TOOL_NAME}"
            tool_output_folder = os.path.join(self.folders["output"], *self._folder_stack, folder_name)
            public_step, internal_order = self.public_step_order, None
        os.makedirs(tool_output_folder, exist_ok=True)
        tool_config.set_pipeline_context(
            self, self.execution_order, tool_output_folder, self.current_suffix,
            public_step=public_step, internal_order=internal_order
        )

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
        self._assign_step_context(tool_config)

        # Configure inputs immediately so tool outputs are properly set
        tool_config.configure_inputs(self.folders)
        
        # Add to pipeline
        self.tools.append(tool_config)
        
        # Create and configure output object
        tool_output = ToolOutput(tool_config)

        # Immediately populate expected outputs using pure path construction
        try:
            expected_outputs = tool_config.get_output_files()
            # Materialize the full sub-layout (configuration/, execution/,
            # tables/, _extras/, and one folder per declared stream) so
            # tool authors never need to mkdir anything themselves. Runs
            # once per tool at config time; idempotent.
            tool_config._materialize_output_layout(expected_outputs)
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
        # Ensure Resources() has been called before adding tools.
        # Resource specs (gpu/memory/time/cpus) are only meaningful for batch
        # schedulers that consume them (slurm/lsf/pbs). In on-the-fly mode, or
        # when the active config uses a non-batch scheduler (colab / none), we
        # auto-initialize an empty resource batch so users don't have to
        # type Resources() for local or Colab runs.
        if self.current_batch == -1:
            scheduler = ConfigManager().get_scheduler()
            if self.on_the_fly or scheduler not in BATCH_SCHEDULERS:
                self.resources()
            else:
                raise RuntimeError(
                    "Resources() must be called before adding tools. "
                    "Use: with Pipeline(...): Resources(...); tool = Tool(...)"
                )

        # If a Parallel() block has just exited, the next tool addition
        # implicitly opens a new fan-in batch. Inherit the most-recent batch
        # resources (per user spec: "default to the resources of the
        # previous one"). _pending_post_parents is consumed inside
        # self.resources() so the new batch's parents are the siblings
        # rather than the chain default.
        if self._pending_post_parents is not None:
            self.resources()  # inherits from self.batch_resources[-1]

        # Merge current batch resources with tool resources
        current_resources = self.batch_resources[self.current_batch]
        for key, value in current_resources.items():
            if value is not None:
                tool_config.resources[key] = value

        # Set execution order and create step-numbered folder
        self._assign_step_context(tool_config)

        # Configure inputs immediately
        tool_config.configure_inputs(self.folders)

        # Add to pipeline
        self.tools.append(tool_config)

        # Create and configure output object
        tool_output = ToolOutput(tool_config)

        # Populate expected outputs and materialize the on-disk sub-layout
        # (canonical sub-dirs + one folder per declared stream). Tools
        # never need to mkdir anything themselves.
        try:
            expected_outputs = tool_config.get_output_files()
            tool_config._materialize_output_layout(expected_outputs)
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
        tool_script_path = os.path.join(self.folders["runtime"], f"{tool_config.script_basename}.sh")
        log_path = os.path.join(self.folders["logs"], f"{tool_config.script_basename}.log")

        # Generate the tool script
        script_content = tool_config.generate_script(tool_script_path)
        with open(tool_script_path, 'w') as f:
            f.write(script_content)
        os.chmod(tool_script_path, 0o755)

        print(f"\n{'='*60}")
        print(f"Running {tool_config.TOOL_NAME} (step {tool_config.execution_order})")
        print(f"{'='*60}")

        # Execute the script, streaming output in real-time
        tool_folder_log_path = os.path.join(tool_config.output_folder, "_log")
        with open(log_path, 'w') as log_file, open(tool_folder_log_path, 'w') as tool_folder_log:
            process = subprocess.Popen(
                ['bash', tool_script_path],
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
                tool_folder_log.write(line)
                tool_folder_log.flush()
            process.wait()

        if process.returncode != 0:
            raise RuntimeError(
                f"{tool_config.TOOL_NAME} failed with exit code {process.returncode}. "
                f"See log: {log_path}"
            )

        print(f"{'='*60}")
        print(f"{tool_config.TOOL_NAME} completed")
        print(f"{'='*60}\n")

    def validate_pipeline(self) -> bool:
        """
        Validate entire pipeline for consistency and compatibility.

        Only runs the deeper dependency / fan-out checks when the pipeline was
        constructed with ``debug=True``; in normal (non-debug) usage we just
        confirm the pipeline is non-empty.

        Returns:
            True if pipeline is valid

        Raises:
            ValueError: If pipeline has issues
        """
        if not self.tools:
            raise ValueError("Pipeline is empty")
        if not self.debug: return
        
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
    
    def save(self) -> str:
        """
        Generate and save pipeline execution script.

        Set ``debug=True`` on the ``Pipeline(...)`` constructor to print per-tool
        outputs here and to embed the runtime debug-capture block in the
        generated bash.

        Returns:
            Path to generated pipeline script
        """
        if not self.tools:
            raise ValueError("Cannot save empty pipeline")

        # Suppress auto-submit on context exit: calling save() directly inside
        # a `with Pipeline(...)` block is an explicit save just like Save().
        self._explicit_save_called = True

        if self.debug:
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

        # Snapshot active config + per-tool inventory for forensic reproducibility
        if self.debug:
            self._write_debug_capture_python_artifacts()

        # Export tool outputs for potential reuse with Load
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
            f'echo "Description: {_escape_for_double_quotes(self.description)}"',
            'echo'
        ]
        
        # Add tool configurations
        for tool in self.tools:
            if tool.internal:
                label = f"Internal {tool.internal_order:03d}: {tool.TOOL_NAME}"
            else:
                label = f"Step {tool.public_step:03d}: {tool.TOOL_NAME}"
            config_lines.extend([
                f'echo "{label}"',
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
        cfg = ConfigManager()
        scheduler_init_lines = cfg.get_scheduler_init()
        module_load = cfg.get_module_load_line()
        script_lines = [
            "#!/bin/bash",
            f"# Generated pipeline: {self.project}",
            f"# Tools: {', '.join([t.TOOL_NAME for t in self.tools])}",
            "",
            "umask 002  # Make all files group-writable by default",
            "",
        ]
        if scheduler_init_lines:
            # Initialise the host's shell environment (sourcing lmod's
            # profile, etc.) so 'module' and any other host primitives
            # exist when this script is executed standalone (i.e. outside
            # the slurm.sh wrapper, which would have done it itself).
            script_lines.append("# Scheduler init (machine.scheduler.init)")
            script_lines.extend(scheduler_init_lines)
            script_lines.append("")
        if module_load:
            # Needed when pipeline.sh is executed directly (outside the
            # slurm.sh wrapper) on hosts where mamba/conda are behind an
            # environment-module system (e.g. S3IT's miniforge3 module).
            script_lines += [
                "# Load cluster modules that expose the env manager on PATH",
                module_load,
                "",
            ]
        shell_hook = cfg.get_shell_hook_command()
        if shell_hook:
            script_lines += [
                "# Initialize environment manager (machine.env_manager.init)",
                shell_hook,
                "",
            ]

        if self.debug:
            script_lines.append('export BIOPIPELINES_DEBUG=1')
            script_lines.append(self._generate_debug_capture_block())
            script_lines.append("")

        script_lines += [
            "echo Configuration",
            f"{config_script} | tee {os.path.join(self.folders['output'], f'{self.project}_config.txt')}",
            "echo"
        ]

        # Add each tool execution (each tool handles its own environment activation)
        for tool in self.tools:
            # Generate tool-specific script
            tool_script_path = os.path.join(self.folders["runtime"], f"{tool.script_basename}.sh")
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
            log_file = os.path.join(self.folders["logs"], f"{tool.script_basename}.log")
            tool_folder_log = os.path.join(tool.output_folder, "_log")
            script_lines.extend([
                f"echo {tool.TOOL_NAME}",
                f"{tool_script_path} 2>&1 | tee {log_file} {tool_folder_log}",
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

    def _iter_pipeline_envs(self):
        """Unique environment names referenced by tools in this pipeline.

        Order is first-seen across self.tools so the export block is
        deterministic and easy to diff between runs.
        """
        seen = set()
        ordered = []
        for tool in self.tools:
            for env in getattr(tool, "environments", []) or []:
                if env and env not in seen:
                    seen.add(env)
                    ordered.append(env)
        return ordered

    def _iter_pipeline_containers(self):
        """List of (tool_name, image_path) for tools that use a container."""
        result = []
        for tool in self.tools:
            try:
                if tool.uses_container():
                    image = self.folders.get(f"container:{tool.TOOL_NAME}", "")
                    result.append((tool.TOOL_NAME, image))
            except Exception:
                continue
        return result

    def _generate_debug_capture_block(self) -> str:
        """Bash that snapshots the runtime onto disk under <output>/_debug_capture/.

        Runs on the node that executes the pipeline script (login when
        on_the_fly, compute under SLURM), so the captured GPU/driver/scheduler
        info reflects the actual execution context. Each command is wrapped
        so missing binaries (no nvidia-smi on a CPU node, no sinfo off SLURM,
        etc.) write a stub instead of failing the run.
        """
        cm = ConfigManager()
        env_manager = cm.get_env_manager()
        scheduler = cm.get_scheduler()
        try:
            container_executor = cm.get_container_executor()
        except Exception:
            container_executor = ""

        capture_root = os.path.join(self.folders["output"], "_debug_capture")
        envs_dir = os.path.join(capture_root, "environments")
        sys_dir = os.path.join(capture_root, "system")

        envs = self._iter_pipeline_envs()
        containers = self._iter_pipeline_containers()

        lines = [
            '# === BioPipelines debug capture ===',
            'if [ "$BIOPIPELINES_DEBUG" = "1" ]; then',
            f'  mkdir -p "{envs_dir}" "{sys_dir}"',
            f'  echo "Capturing runtime environment to {capture_root}"',
            '',
            '  # System info (best-effort: missing binaries write a stub)',
            f'  ( uname -a ) > "{sys_dir}/uname.txt" 2>&1 || echo "uname not available" > "{sys_dir}/uname.txt"',
            f'  ( command -v nvidia-smi >/dev/null && nvidia-smi ) > "{sys_dir}/nvidia-smi.txt" 2>&1 || echo "nvidia-smi not available" > "{sys_dir}/nvidia-smi.txt"',
            f'  ( {env_manager} --version ) > "{sys_dir}/env_manager.txt" 2>&1 || echo "{env_manager} not available" > "{sys_dir}/env_manager.txt"',
        ]

        if scheduler == "slurm":
            lines.append(f'  ( command -v sinfo >/dev/null && sinfo -V ) > "{sys_dir}/scheduler.txt" 2>&1 || echo "sinfo not available" > "{sys_dir}/scheduler.txt"')
        else:
            lines.append(f'  echo "scheduler={scheduler}" > "{sys_dir}/scheduler.txt"')

        if containers and container_executor:
            lines.append(f'  ( {container_executor} --version ) > "{sys_dir}/containers.txt" 2>&1 || echo "{container_executor} not available" > "{sys_dir}/containers.txt"')
        else:
            lines.append(f'  echo "no containers configured for any tool in this pipeline" > "{sys_dir}/containers.txt"')

        lines.append('')
        lines.append('  # Per-environment exports (one section per unique env across tools)')

        if env_manager == "pip":
            lines.append(f'  ( pip freeze ) > "{envs_dir}/pip.txt" 2>&1 || echo "pip not available" > "{envs_dir}/pip.txt"')
        else:
            for env in envs:
                env_yaml = f'{envs_dir}/{env}.yaml'
                env_pip = f'{envs_dir}/{env}.pip.txt'
                lines.append(f'  ( {env_manager} env export --no-builds -n {env} ) > "{env_yaml}" 2>&1 || echo "env export failed for {env}" > "{env_yaml}"')
                lines.append(f'  ( {env_manager} run -n {env} pip freeze ) > "{env_pip}" 2>&1 || echo "pip freeze failed for {env}" > "{env_pip}"')

        lines.append('  echo "Debug capture complete."')
        lines.append('fi')
        lines.append('# === end debug capture ===')
        return "\n".join(lines)

    def _write_debug_capture_python_artifacts(self):
        """Write Python-side debug artifacts at save() time.

        The bash block captures runtime state on the execution node; this
        writes the static, configuration-time artifacts (active config copy,
        per-tool inventory) immediately so they exist even if the pipeline
        is never actually executed.
        """
        capture_root = os.path.join(self.folders["output"], "_debug_capture")
        os.makedirs(capture_root, exist_ok=True)

        cm = ConfigManager()
        variant = cm.get_variant()

        # Snapshot the active config file
        try:
            src = cm._get_config_path()
            dst = os.path.join(capture_root, f"config.{variant}.yaml")
            shutil.copy2(src, dst)
        except Exception as e:
            with open(os.path.join(capture_root, "config_snapshot_error.txt"), "w") as f:
                f.write(f"Could not copy active config: {e}\n")

        # Per-tool inventory
        try:
            scheduler = cm.get_scheduler()
        except Exception:
            scheduler = "unknown"
        try:
            container_executor = cm.get_container_executor()
        except Exception:
            container_executor = ""

        tools_inventory = []
        for tool in self.tools:
            entry = {
                "tool_name": getattr(tool, "TOOL_NAME", tool.__class__.__name__),
                "tool_version": getattr(tool, "TOOL_VERSION", None),
                "class": tool.__class__.__name__,
                "environments": list(getattr(tool, "environments", []) or []),
                "container_image": "",
            }
            try:
                if tool.uses_container():
                    entry["container_image"] = self.folders.get(
                        f"container:{tool.TOOL_NAME}", ""
                    )
            except Exception:
                pass
            tools_inventory.append(entry)

        try:
            from . import __version__ as _bp_version
        except Exception:
            _bp_version = "unknown"

        inventory = {
            "project": self.project,
            "job": self.job,
            "variant": variant,
            "biopipelines_version": _bp_version,
            "env_manager": cm.get_env_manager() if hasattr(cm, "get_env_manager") else "",
            "scheduler": scheduler,
            "container_executor": container_executor,
            "tools": tools_inventory,
        }
        with open(os.path.join(capture_root, "tools.json"), "w") as f:
            json.dump(inventory, f, indent=2)

        # README
        readme = (
            "BioPipelines debug capture\n"
            "==========================\n"
            "Static artifacts written at save() time when Pipeline(debug=True):\n"
            "  - config.<variant>.yaml: copy of the active config file\n"
            "  - tools.json: biopipelines version, per-tool wrapper version,\n"
            "    environments, container images, scheduler\n"
            "Runtime artifacts written by the pipeline script when executed:\n"
            "  - environments/<env>.yaml + <env>.pip.txt: env exports per unique env\n"
            "  - system/{uname,nvidia-smi,env_manager,scheduler,containers}.txt\n"
            "Captures runtime state on the executing node (login if on_the_fly,\n"
            "compute under SLURM); not guaranteed to reinstall elsewhere.\n"
        )
        with open(os.path.join(capture_root, "README.txt"), "w") as f:
            f.write(readme)

    def generate_job_scripts(self, email: str = "auto"):
        """
        Generate batch job submission script(s) for the active scheduler.

        If multiple batches exist (multiple Resources() calls), generates:
        - Separate batch script for each batch
        - Dependency placeholders so the submit wrapper chains them

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
            raise RuntimeError("Resources() must be called before generating job scripts")

        # Handle auto email detection
        if email == "auto":
            email = ConfigManager().get_email()
            if email:
                print(f"Auto-detected email: {email}")
            else:
                print("Warning: No email configured. Disabling email notifications.")
                print("Set the 'machine.email' field in your config file.")

        email_line = self._backend().email_directive(email)

        num_batches = len(self.batch_resources)

        # Determine if single batch or multiple batches
        if num_batches == 1:
            self._generate_single_batch(email_line)
        else:
            self._generate_multi_batch(email_line, num_batches)

    def _backend(self):
        """Scheduler backend used for script generation.

        For a batch scheduler (slurm/lsf/pbs) this is its own backend. For
        non-batch schedulers (colab/none) we still emit an inert script (the
        submission step is skipped, and printing is suppressed) — fall back to
        the SLURM backend so that script keeps its historical #SBATCH form.
        """
        name = ConfigManager().get_scheduler()
        if name in BATCH_SCHEDULERS:
            return get_backend(name)
        return get_backend("slurm")

    def _is_batch_scheduler(self) -> bool:
        return ConfigManager().get_scheduler() in BATCH_SCHEDULERS

    def _scheduler_options_key(self) -> str:
        """Resource-dict key for extra directive kwargs of the active scheduler.

        ``slurm_options`` / ``lsf_options`` / ``pbs_options``. Falls back to
        ``slurm_options`` for non-batch schedulers (where it is never read).
        """
        if self._is_batch_scheduler():
            return self._backend().options_key
        return "slurm_options"

    def _batch_script_name(self, batch_idx: Optional[int] = None) -> str:
        """Generated script filename, stemmed by the active scheduler.

        SLURM keeps the historical ``slurm.sh`` / ``slurm_batch{n}.sh`` names;
        LSF/PBS use ``lsf.sh`` / ``pbs_batch{n}.sh`` etc. The submit wrapper
        globs the matching stem per scheduler. Non-batch schedulers
        (colab/none) reuse the ``slurm`` stem for their inert scripts.
        """
        name = ConfigManager().get_scheduler()
        stem = name if name in BATCH_SCHEDULERS else "slurm"
        if batch_idx is None:
            return f"{stem}.sh"
        return f"{stem}_batch{batch_idx + 1}.sh"

    def _generate_single_batch(self, email_line):
        """Generate single batch script for the active scheduler."""
        backend = self._backend()
        resources = self.batch_resources[0]
        _validate_directive_value("memory", resources["memory"])
        _validate_directive_value("time", resources["time"])
        gpu_line, gpu_warnings = backend.gpu_directive(resources["gpu"], resources.get("gpus", 1))
        gpu_setup = backend.gpu_setup(resources["gpu"])
        extra_directives = backend.extra_options(resources.get(backend.options_key, {}))
        header = backend.header_directives(resources["memory"], resources["time"], "job.out")
        dependency_line = backend.dependency_directive(list(self.external_dependencies), [])
        cm = ConfigManager()
        scheduler_init = cm.get_scheduler_init_block()
        module_load = cm.get_module_load_line()

        job_content = f"""#!/usr/bin/bash
{gpu_line}
{header}{extra_directives}{dependency_line}
{email_line}

# Make all files group-writable by default
umask 002
# Scheduler init (machine.scheduler.init): runs first so 'module' and
# any other host-level shell primitives are available before module load.
{scheduler_init}
{gpu_setup}
{module_load}

# Execute pipeline
{self.pipeline_script}
"""
        job_path = os.path.join(self.folders["runtime"], self._batch_script_name())
        with open(job_path, 'w') as f:
            f.write(job_content)
        os.chmod(job_path, 0o755)
        self.job_script = job_path

        # Batch scripts are meaningless without a real scheduler (Colab) — skip printing.
        if not self._is_batch_scheduler():
            return

        for w in gpu_warnings:
            print(f"Warning: {w}")
        print(f"Job script saved to: {job_path}")
        print("="*30+"Job"+"="*30)
        print(f"{self.project}: {self.job} ({self.folder_manager.job_id})")
        print("="*30+"Job Script"+"="*30)
        for line in job_content.split('\n'):
            if line != "":
                print(line)

    def _generate_multi_batch(self, email_line, num_batches):
        """Generate chained batch scripts with <JOBID_BATCH_NNN> dependency placeholders (one per parent in the DAG)."""
        backend = self._backend()
        cm = ConfigManager()
        scheduler_init = cm.get_scheduler_init_block()
        module_load = cm.get_module_load_line()
        warnings: List[str] = []

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

        # Generate per-batch scheduler scripts
        for batch_idx in range(num_batches):
            resources = self.batch_resources[batch_idx]
            _validate_directive_value("memory", resources["memory"])
            _validate_directive_value("time", resources["time"])
            gpu_line, gpu_warnings = backend.gpu_directive(resources["gpu"], resources.get("gpus", 1))
            warnings.extend(gpu_warnings)
            gpu_setup = backend.gpu_setup(resources["gpu"])
            extra_directives = backend.extra_options(resources.get(backend.options_key, {}))
            header = backend.header_directives(resources["memory"], resources["time"], f"job_batch{batch_idx + 1}.out")

            # Dependency directive for batch 2+ (internal DAG) and batch 1 (external).
            # Internal: emit one <JOBID_BATCH_NNN> placeholder per parent batch.
            # The submit script captures every batch's job id keyed by its
            # 1-based number and substitutes each placeholder when the dependent
            # batch is submitted (parents always have lower numbers, so they are
            # already in the captured map by the time we substitute).
            parents = self.batch_parents[batch_idx] if batch_idx < len(self.batch_parents) else []
            after_parents = self.batch_after_parents[batch_idx] if batch_idx < len(self.batch_after_parents) else []
            if parents or after_parents:
                afterok = [f"<JOBID_BATCH_{p + 1:03d}>" for p in parents]
                # Service: start once the daemon is RUNNING (after:), not finished.
                after = [f"<JOBID_BATCH_{p + 1:03d}>" for p in after_parents]
                dependency_line = backend.dependency_directive(afterok, after)
            elif batch_idx == 0 and self.external_dependencies:
                dependency_line = backend.dependency_directive(list(self.external_dependencies), [])
            else:
                dependency_line = ""

            batch_job_content = f"""#!/usr/bin/bash
{gpu_line}
{header}{extra_directives}{dependency_line}
{email_line}

# Make all files group-writable by default
umask 002
# Scheduler init (machine.scheduler.init): runs first so 'module' and
# any other host-level shell primitives are available before module load.
{scheduler_init}
{gpu_setup}
{module_load}

# Execute pipeline batch {batch_idx + 1}
{os.path.join(self.folders["runtime"], f"pipeline_batch{batch_idx + 1}.sh")}
"""
            batch_job_path = os.path.join(self.folders["runtime"], self._batch_script_name(batch_idx))
            with open(batch_job_path, 'w') as f:
                f.write(batch_job_content)
            os.chmod(batch_job_path, 0o755)

        self.job_script = os.path.join(self.folders["runtime"], self._batch_script_name(0))

        # Batches are meaningless without a real scheduler (Colab) — skip printing.
        if not self._is_batch_scheduler():
            return

        for w in warnings:
            print(f"Warning: {w}")
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
        print("# Submit batches manually (replace each <JOBID_BATCH_NNN> placeholder with the captured job ID for batch NNN):")
        for batch_idx in range(num_batches):
            print(f"# {os.path.join(self.folders['runtime'], self._batch_script_name(batch_idx))}")
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
            f'echo "Description: {_escape_for_double_quotes(self.description)}"',
            'echo',
            f'echo "Tools in this batch: {len(batch_tools)}"',
            'echo'
        ]

        for tool in batch_tools:
            if tool.internal:
                label = f"Internal {tool.internal_order:03d}: {tool.TOOL_NAME}"
            else:
                label = f"Step {tool.public_step:03d}: {tool.TOOL_NAME}"
            config_lines.extend([
                f'echo "{label}"',
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
        cfg = ConfigManager()
        scheduler_init_lines = cfg.get_scheduler_init()
        module_load = cfg.get_module_load_line()
        script_lines = [
            "#!/bin/bash",
            f"# Generated pipeline batch {batch_idx + 1}: {self.project}",
            f"# Tools: {', '.join([t.TOOL_NAME for t in batch_tools])}",
            "",
            "umask 002  # Make all files group-writable by default",
            "",
        ]
        if scheduler_init_lines:
            script_lines.append("# Scheduler init (machine.scheduler.init)")
            script_lines.extend(scheduler_init_lines)
            script_lines.append("")
        if module_load:
            script_lines += [
                "# Load cluster modules that expose the env manager on PATH",
                module_load,
                "",
            ]
        shell_hook = cfg.get_shell_hook_command()
        if shell_hook:
            script_lines += [
                "# Initialize environment manager (machine.env_manager.init)",
                shell_hook,
                "",
            ]

        if self.debug and batch_idx == 0:
            script_lines.append('export BIOPIPELINES_DEBUG=1')
            script_lines.append(self._generate_debug_capture_block())
            script_lines.append("")

        script_lines += [
            "echo Configuration",
            f"{config_script} | tee -a {os.path.join(self.folders['output'], f'{self.project}_config.txt')}",
            "echo"
        ]

        # Add each tool execution (each tool handles its own environment activation)
        for tool in batch_tools:
            # Tool script path
            tool_script_path = os.path.join(self.folders["runtime"], f"{tool.script_basename}.sh")
            log_file = os.path.join(self.folders["logs"], f"{tool.script_basename}.log")
            tool_folder_log = os.path.join(tool.output_folder, "_log")

            script_lines.extend([f"echo {tool.TOOL_NAME}", f"{tool_script_path} 2>&1 | tee {log_file} {tool_folder_log}", "echo"])

        # Final steps
        script_lines.extend(["echo", f"echo Batch {batch_idx + 1} done"])
        if batch_idx == len(self.batch_resources) - 1:  # Last batch
            script_lines.extend(["echo", "echo All jobs complete", f"echo Results in:", f"echo {self.folders['output']}"])

        return "\n".join(script_lines)

    def set_dependencies(self, job_ids: Union[str, List[str]]):
        """
        Set external scheduler job dependencies.

        The first batch of this pipeline will wait for these jobs to complete
        successfully before starting.

        Args:
            job_ids: Single job ID string or list of job ID strings
                    e.g., "12345678" or ["12345678", "12345679"]
        """
        if isinstance(job_ids, str):
            job_ids = [job_ids]
        self.external_dependencies.extend(job_ids)

    def resources(self, gpu: str = None, memory: str = None, time: str = None, cpus: int = None, gpus: int = None, **scheduler_options):
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
            cpus: CPUs per task → #SBATCH --cpus-per-task=<cpus>.
                - None: Inherits from previous batch (omitted on first batch)
            gpus: Number of GPUs → #SBATCH --gpus=<gpus> (combined with any
                gpu= model/memory constraint). Default behaviour (None) is 1 GPU
                whenever a gpu spec is set. Use gpus=2 to request two GPUs.
            **scheduler_options: Additional native directive parameters for the
                active scheduler, stored under its options key (slurm_options /
                lsf_options / pbs_options). On SLURM:
                - nodes=1 → #SBATCH --nodes=1
                - ntasks_per_node=1 → #SBATCH --ntasks-per-node=1
                - partition="gpu" → #SBATCH --partition=gpu
                On LSF/PBS the key/value is emitted as the scheduler's native
                flag (e.g. LSF q="normal" → #BSUB -q normal).

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
        # Extra kwargs are native directive knobs for the active scheduler;
        # store them under that scheduler's options key (slurm/lsf/pbs_options)
        # so the matching backend emits them. Defaults to slurm_options for
        # non-batch schedulers (colab/none), where it is never read.
        options_key = self._scheduler_options_key()

        # Reject newline/CR before anything else — these would inject extra
        # directive lines at script generation time.
        _validate_directive_value("gpu", gpu)
        _validate_directive_value("memory", memory)
        _validate_directive_value("time", time)
        _validate_directive_value("cpus", cpus)
        for opt_key, opt_value in scheduler_options.items():
            _validate_directive_value(f"{options_key}[{opt_key!r}]", opt_value)

        if gpus is not None and (not isinstance(gpus, int) or gpus < 1):
            raise ValueError(f"gpus must be a positive integer, got {gpus!r}")

        # Start new batch
        self.current_batch += 1
        self.batch_start_indices.append(len(self.tools))

        # Determine the new batch's parents in the dependency DAG. Two
        # orthogonal edge kinds are computed independently and may coexist
        # (e.g. siblings of a Parallel block that follows a Service):
        #
        #   after_parents (SLURM `after:`): a batch that follows a Service
        #     waits for the daemon to START, not finish. For a plain batch
        #     this is the single _pending_after_parent. For a Parallel block
        #     that follows a Service, EVERY sibling gets the edge, carried on
        #     _parallel_after_anchor. The daemon must never also be an afterok
        #     parent, so it is excluded from new_parents below.
        #
        #   new_parents (SLURM `afterok:`), in priority order:
        #     1. fan-in pending from a Parallel block that just exited.
        #     2. inside a Parallel block: parent is the anchor (the batch
        #        current when the block was entered). Record this sibling for
        #        later fan-in. An anchor that is the service daemon is dropped
        #        — that edge is expressed via after_parents instead.
        #     3. chain default: parent is the previous batch ([] for batch 0).
        after_parents = []  # type: List[int]
        if self._parallel_anchor is not None:
            if self._parallel_after_anchor is not None:
                after_parents = [self._parallel_after_anchor]
        elif self._pending_after_parent is not None:
            after_parents = [self._pending_after_parent]
            self._pending_after_parent = None

        if self._pending_post_parents is not None:
            new_parents = list(self._pending_post_parents)
            self._pending_post_parents = None
        elif self._parallel_anchor is not None:
            self._parallel_siblings.append(self.current_batch)
            if self._parallel_anchor >= 0 and self._parallel_anchor not in after_parents:
                new_parents = [self._parallel_anchor]
            else:
                new_parents = []
        elif after_parents:
            # Plain (non-parallel) batch after a Service: daemon is its only
            # predecessor via after:, so no afterok parent of its own.
            new_parents = []
        else:
            new_parents = [self.current_batch - 1] if self.current_batch > 0 else []
        self.batch_parents.append(new_parents)
        self.batch_after_parents.append(after_parents)

        # Fold first-class cpus into the options dict (maps to --cpus-per-task
        # on SLURM; passed through verbatim on other schedulers).
        if cpus is not None:
            scheduler_options = {**scheduler_options, "cpus": cpus}

        # Determine resources for this batch
        if self.current_batch == 0:
            # First batch: use provided values or defaults
            batch_res = {
                "gpu": gpu if gpu is not None else None,
                "gpus": gpus if gpus is not None else 1,
                "memory": memory if memory is not None else "15GB",
                "time": time if time is not None else "24:00:00",
                options_key: scheduler_options if scheduler_options else {}
            }
        else:
            # Subsequent batches: inherit from previous if not specified
            prev_res = self.batch_resources[self.current_batch - 1]
            prev_opts = prev_res.get(options_key, {})
            if scheduler_options:
                merged_opts = {**prev_opts, **scheduler_options}
            else:
                merged_opts = prev_opts
            batch_res = {
                "gpu": gpu if gpu is not None else prev_res["gpu"],
                "gpus": gpus if gpus is not None else prev_res.get("gpus", 1),
                "memory": memory if memory is not None else prev_res["memory"],
                "time": time if time is not None else prev_res["time"],
                options_key: merged_opts
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

        for tool in self.tools:
            n = tool.internal_order if tool.internal else tool.public_step
            tag = f".{n:03d}" if tool.internal else f"{n:03d}"
            summary.append(
                f"{tag}. {tool.TOOL_NAME} ({tool.environments}) -> {tool.job_name}"
            )

        return "\n".join(summary)
    
    def _export_tool_outputs(self):
        """
        Export tool output metadata to JSON files for potential reuse with Load.
        
        Creates a tool_outputs/ directory with JSON files for each tool containing
        complete output structure, configuration, and execution metadata.
        """
        from datetime import datetime
        
        # Create tool_outputs directory
        tool_outputs_dir = os.path.join(self.folders["output"], "ToolOutputs")
        os.makedirs(tool_outputs_dir, exist_ok=True)
        
        exported_count = 0
        
        for tool, tool_output in zip(self.tools, self.tool_outputs):
            try:
                # Get current output structure
                output_structure = tool.get_output_files()

                # Convert DataStream objects to dictionaries for JSON serialization
                from .datastream import DataStream
                from .base_config import TableInfo

                for key, value in list(output_structure.items()):
                    if isinstance(value, DataStream):
                        output_structure[key] = value.to_dict()

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
                    "execution_order": tool.execution_order,
                    "internal": tool.internal,
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
                
                # Filename mirrors the tool's script basename (internal tools nest under .internal/)
                output_path = os.path.join(tool_outputs_dir, f"{tool.script_basename}.json")
                os.makedirs(os.path.dirname(output_path), exist_ok=True)

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
        cpus: CPUs per task (SLURM --cpus-per-task; translated per scheduler)

    Example:
        with Pipeline("Test", "Job", "Description"):
            Resources(gpu="32GB", memory="16GB", time="24:00:00", cpus=8)
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
    Save pipeline configuration without generating job scripts.

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
    # save() itself sets _explicit_save_called so auto-submit is suppressed.
    pipeline.save()


def Dependencies(job_ids: Union[str, List[str]]):
    """
    Set external scheduler job dependencies for the active pipeline.

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
    if pipeline._parallel_anchor is not None:
        raise RuntimeError(
            "Dependencies() cannot be called inside a `with Parallel():` "
            "block. The block's structural rules already determine "
            "intra-block (sibling-with-shared-parent) and post-block "
            "(fan-in) dependencies; an explicit Dependencies() call would "
            "be ambiguous. Move the Dependencies() call outside the block."
        )
    if pipeline._service_anchor is not None:
        raise RuntimeError(
            "Dependencies() cannot be called inside a `with Service():` block."
        )
    pipeline.set_dependencies(job_ids)


class Parallel:
    """Context manager that runs the batches inside it as parallel siblings.

    Use it to express a fan-out / fan-in pattern. Every Resources() call
    inside the block opens a batch whose only dependency is the
    batch that was current immediately before the block (the "anchor"),
    so the sibling batches run in parallel rather than chained. The first
    Resources() call after the block exits opens a fan-in batch that
    waits for all siblings to finish. If the block sits at the top of the
    pipeline (no batch before it), the siblings simply have no upstream
    dependency.

    Each iteration MUST call Resources() to open its own sibling batch.
    Tools (including entity tools used as kwargs, e.g. Ligand("LIG")
    inside BoltzGen(ligand=Ligand("LIG"), ...)) then land in that batch.

    Example
    -------
    Run RFdiffusion 10 times in parallel and gather the outputs with Pool::

        with Pipeline(...):
            Resources(...)
            PreLoopTool(...)                  # batch 1
            runs = []
            with Parallel():
                for i in range(10):
                    Resources(gpu="A100")     # batches 2..11, all depend on batch 1
                    runs.append(RFdiffusion(...))
            Resources(...)                    # batch 12, depends on batches 2..11
            combined = Pool(runs=runs)        # 100 dense ids + pool.path

    Notes
    -----
    * Dependencies() inside the block raises — the structural rules
      already determine all dependencies.
    * Nesting Parallel() blocks is not supported in this version.
    """

    def __enter__(self):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is None:
            raise RuntimeError(
                "Parallel() must be used within a Pipeline context. "
                "Use: with Pipeline(...): with Parallel(): ..."
            )
        if pipeline._parallel_anchor is not None:
            raise RuntimeError(
                "Nested Parallel() blocks are not supported."
            )
        pipeline._parallel_anchor = pipeline.current_batch
        pipeline._parallel_siblings = []
        # If this block immediately follows a Service, every sibling waits for
        # the daemon to START via after: (not afterok on the anchor). Absorb
        # the single-use service edge so all siblings share it.
        pipeline._parallel_after_anchor = pipeline._pending_after_parent
        pipeline._pending_after_parent = None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is None:
            return False
        siblings = list(pipeline._parallel_siblings)
        after_anchor = pipeline._parallel_after_anchor
        pipeline._parallel_anchor = None
        pipeline._parallel_siblings = []
        pipeline._parallel_after_anchor = None
        if siblings:
            pipeline._pending_post_parents = siblings
        else:
            # Empty block: a transparent no-op. Restore the pending service
            # after: edge it absorbed in __enter__ so the next batch still
            # waits for the daemon. Otherwise fall back to chain default.
            if after_anchor is not None:
                pipeline._pending_after_parent = after_anchor
        return False  # don't suppress exceptions


class Service:
    """Context manager that launches its single batch as a detached daemon.

    Use it for a long-running service (e.g. an MSA server) that a later step
    consumes while it is still running. The block must contain exactly one
    batch (one ``Resources()`` + the daemon tool). That batch keeps its normal
    upstream ``afterok`` dependency, but it is NOT waited on by the chain: the
    first batch AFTER the block gets an ``after:`` dependency on it (SLURM
    naming; translated per scheduler) — "start once the daemon is RUNNING",
    not "after it finishes". ``afterok``
    would be wrong here, since the daemon only exits on its own idle-timeout.

    This removes the scheduling race where a separately-submitted server sits
    queued behind low priority while the client's allocation wall-clock runs
    out: with ``after:`` the client is not scheduled until the server is up.

    Example
    -------
    ::

        with Pipeline(...):
            Resources(...)
            seqs = LigandMPNN(...)            # batch A
            with Service():
                Resources(memory="900GB", cpus=32, time="24:00:00")
                MMseqs2Server(mode="cpu")     # batch B (daemon), afterok:A
            Resources(memory="16GB")
            msas = MMseqs2(sequences=seqs)    # batch C, after:B
            Resources(gpu="A100")
            Boltz2(..., msas=msas)            # batch D, afterok:C

    Notes
    -----
    * Exactly one batch inside the block (one ``Resources()`` + daemon tool).
    * The daemon is fully detached: nothing ``afterok``-waits on it, so the
      pipeline's success never hinges on its exit code. It self-terminates via
      its own idle-timeout once the consuming step has drained its queue.
    * Cannot be nested, cannot contain ``Dependencies()``, and cannot be used
      inside a ``Parallel()`` block.
    """

    def __enter__(self):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is None:
            raise RuntimeError(
                "Service() must be used within a Pipeline context. "
                "Use: with Pipeline(...): with Service(): ..."
            )
        if pipeline._service_anchor is not None:
            raise RuntimeError("Nested Service() blocks are not supported.")
        if pipeline._parallel_anchor is not None:
            raise RuntimeError("Service() cannot be used inside a Parallel() block.")
        # Remember how many batches existed on entry, to enforce exactly one.
        pipeline._service_anchor = pipeline.current_batch
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is None:
            return False
        anchor = pipeline._service_anchor
        pipeline._service_anchor = None
        if exc_type is not None:
            return False  # don't mask an in-block exception
        opened = pipeline.current_batch - anchor
        if opened != 1:
            raise RuntimeError(
                "A Service() block must contain exactly one batch (one "
                f"Resources() + the daemon tool); it opened {opened}."
            )
        # The daemon is the batch opened inside the block; the next batch
        # after the block gets an after: edge to it.
        pipeline._pending_after_parent = pipeline.current_batch
        return False  # don't suppress exceptions


class Folder:
    """Context manager that nests public tool outputs under a named subfolder.

    Purely organizational: it prepends a path segment to the output folder of
    every public tool created inside the block, leaving execution order,
    resources, batching, and dependencies untouched. The global step counter
    keeps running, so numbers still reflect true execution order::

        with Pipeline(...):
            Resources()
            Tool1()                 # 001_Tool1/
            with Folder("group"):
                Tool2()             # group/002_Tool2/
            Tool3()                 # 003_Tool3/

    Blocks nest (``with Folder("a"): with Folder("b"):`` -> ``a/b/...``).
    Internal tools ignore the stack; they always land under .internal/.

    In on-the-fly / Colab runs the block can be bound and zipped for download::

        with Folder("Results") as results:
            ...tools...
        results.download()   # zips Results/ and triggers a Colab download
    """

    def __init__(self, name: str):
        self.name = name
        self.path = None

    def __enter__(self):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is None:
            raise RuntimeError(
                "Folder() must be used within a Pipeline context. "
                "Use: with Pipeline(...): with Folder('...'): ..."
            )
        _validate_identifier("Folder name", self.name)
        if self.name == INTERNAL_FOLDER:
            raise ValueError(f"{INTERNAL_FOLDER!r} is reserved for framework-internal tools.")
        # Absolute path of this folder, including any enclosing Folder() blocks.
        self.path = os.path.join(pipeline.folders["output"], *pipeline._folder_stack, self.name)
        self._on_the_fly = pipeline.on_the_fly
        pipeline._folder_stack.append(self.name)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pipeline = Pipeline.get_active_pipeline()
        if pipeline is not None and pipeline._folder_stack:
            pipeline._folder_stack.pop()
        return False  # don't suppress exceptions

    def download(self) -> str:
        """Zip this folder and trigger a browser download in Colab.

        Meaningful only in on-the-fly / Colab runs, where tools execute as they
        are added so the folder already holds real files when the block exits.
        Returns the path to the created .zip. Outside Colab the archive is still
        written and its path printed.
        """
        if not self._on_the_fly:
            raise RuntimeError(
                "Folder.download() is only available in on-the-fly / Colab runs; "
                "in submit mode the outputs do not exist yet at this point."
            )
        if not os.path.isdir(self.path):
            raise FileNotFoundError(f"Folder {self.path!r} has no outputs to download.")
        archive = shutil.make_archive(self.path, "zip", self.path)
        try:
            from google.colab import files
            files.download(archive)
        except ImportError:
            print(f"Not running in Colab; archive written to {archive}")
        return archive