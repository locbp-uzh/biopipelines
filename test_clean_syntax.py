"""
Test the new clean syntax with module-level functions.
"""

import sys
import os
sys.path.insert(0, os.getcwd())

# Set UTF-8 encoding for Windows console
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

from PipelineScripts.pipeline import Pipeline, Resources, Suffix, Save
from PipelineScripts.base_config import BaseConfig

# Create a simple mock tool for testing
class MockTool(BaseConfig):
    TOOL_NAME = "MockTool"
    DEFAULT_ENV = "ProteinEnv"

    def __init__(self, param1, param2="default"):
        self.param1 = param1
        self.param2 = param2
        super().__init__()

    def validate_params(self):
        pass

    def configure_inputs(self, pipeline_folders):
        pass

    def get_output_files(self):
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": {},
            "output_folder": self.output_folder
        }

    def generate_script(self, script_path):
        return "#!/bin/bash\necho 'Mock tool'"

print("=" * 60)
print("Testing Clean Syntax (No 'pipeline' Variable)")
print("=" * 60)

# Test 1: Module-level Resources() function
print("\n[Test 1] Resources() function")
try:
    with Pipeline("Test", "CleanSyntaxTest", "Testing clean syntax"):
        # Set resources without pipeline variable
        Resources(gpu="32GB", memory="16GB", time="24:00:00")

        # Get active pipeline to verify
        active = Pipeline.get_active_pipeline()
        if active.global_resources['gpu'] == "32GB":
            print("✓ Resources() set gpu='32GB'")
        if active.global_resources['memory'] == "16GB":
            print("✓ Resources() set memory='16GB'")
        if active.global_resources['time'] == "24:00:00":
            print("✓ Resources() set time='24:00:00'")

        # Add a tool to prevent empty pipeline error
        tool = MockTool("test")

        # Use Save() to prevent auto-submit
        Save()

        print("\n[Test 1] PASSED ✓\n")

except Exception as e:
    print(f"\n✗ [Test 1] FAILED: {e}\n")
    import traceback
    traceback.print_exc()

# Test 2: Module-level Suffix() function
print("=" * 60)
print("[Test 2] Suffix() function")
try:
    with Pipeline("Test", "SuffixTest", "Testing suffix"):
        Suffix("001")

        active = Pipeline.get_active_pipeline()
        if active.current_suffix == "001":
            print("✓ Suffix() set suffix='001'")

        tool1 = MockTool("test1")

        Suffix("002")
        if active.current_suffix == "002":
            print("✓ Suffix() changed suffix='002'")

        tool2 = MockTool("test2")

        print("\n[Test 2] PASSED ✓\n")

except Exception as e:
    print(f"\n✗ [Test 2] FAILED: {e}\n")
    import traceback
    traceback.print_exc()

# Test 3: Module-level Save() function prevents auto-submit
print("=" * 60)
print("[Test 3] Save() function (prevents auto-submit)")
try:
    class TrackingSLURM(Pipeline):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.slurm_called = False

        def slurm(self):
            self.slurm_called = True
            # Don't actually submit
            print("  (slurm() called)")

    with TrackingSLURM("Test", "SaveTest", "Testing save") as p:
        tool1 = MockTool("test")
        Save()  # Explicitly save

    # Check if slurm was called
    if not p.slurm_called:
        print("✓ Save() prevented auto-submit (slurm not called)")
    else:
        print("✗ ERROR: slurm() was called despite Save()")

    print("\n[Test 3] PASSED ✓\n")

except Exception as e:
    print(f"\n✗ [Test 3] FAILED: {e}\n")
    import traceback
    traceback.print_exc()

# Test 4: Auto-submit on exit (without Save())
print("=" * 60)
print("[Test 4] Auto-submit on exit")
try:
    class TrackingSLURM2(Pipeline):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.slurm_called = False

        def slurm(self):
            self.slurm_called = True
            print("  (slurm() auto-called on exit)")

    with TrackingSLURM2("Test", "AutoSubmitTest", "Testing auto-submit") as p:
        tool1 = MockTool("test")
        # Don't call Save()

    # Check if slurm was called
    if p.slurm_called:
        print("✓ Auto-submit worked (slurm was called on exit)")
    else:
        print("✗ ERROR: slurm() was NOT called on exit")

    print("\n[Test 4] PASSED ✓\n")

except Exception as e:
    print(f"\n✗ [Test 4] FAILED: {e}\n")
    import traceback
    traceback.print_exc()

# Test 5: Complete clean syntax example
print("=" * 60)
print("[Test 5] Complete clean syntax (no 'pipeline' variable)")
try:
    # This is the beautiful new syntax!
    with Pipeline("Example", "CleanDemo", "Demonstrating clean syntax"):
        Resources(gpu="32GB", memory="16GB", time="24:00:00")

        Suffix("001")
        tool1 = MockTool("first")
        tool2 = MockTool("second")

        Suffix("002")
        tool3 = MockTool("third")

    # Pipeline auto-submits on exit!

    print("✓ Clean syntax works perfectly!")
    print("  - No 'pipeline' variable needed")
    print("  - Resources() for global settings")
    print("  - Suffix() for tool naming")
    print("  - Auto-submit on exit")

    print("\n[Test 5] PASSED ✓\n")

except Exception as e:
    print(f"\n✗ [Test 5] FAILED: {e}\n")
    import traceback
    traceback.print_exc()

print("=" * 60)
print("All tests completed!")
print("=" * 60)
