# BioPipelines Refactoring Proposal

## Executive Summary

This document outlines a comprehensive refactoring strategy to improve code elegance, reduce duplication, and enhance maintainability across the BioPipelines tool system using decorators, mixins, and helper classes.

## Current State Analysis

### Identified Issues

1. **Excessive Code Duplication**
   - Every tool implements similar `_initialize_file_paths()` and `_setup_file_paths()` methods
   - Input handling logic repeated across 30+ tools
   - Datasheet access patterns duplicated everywhere

2. **Verbose Boilerplate**
   - 50-100 lines of standard initialization code per tool
   - Repetitive input type checking (ToolOutput vs StandardizedOutput vs str vs list)
   - Redundant dependency management code

3. **Complex Datasheet Navigation**
   - Multiple formats to handle (dict, DatasheetContainer, legacy list)
   - Nested if-else chains for accessing datasheets
   - Fragile path extraction logic

4. **Script Generation Verbosity**
   - Repeated completion check headers/footers
   - Similar bash script patterns across tools
   - Manual string concatenation

## Proposed Solutions

### 1. Decorator-Based File Path Management

**Problem**: Every tool manually defines file paths in `_initialize_file_paths()` and `_setup_file_paths()`

**Solution**: Use class-level decorators to auto-generate file path management

```python
class FilePathDescriptor:
    """Descriptor for auto-managed file paths."""

    def __init__(self, filename_template, folder_key="output_folder"):
        self.filename_template = filename_template
        self.folder_key = folder_key
        self._attr_name = None

    def __set_name__(self, owner, name):
        self._attr_name = f"_{name}"

    def __get__(self, instance, owner):
        if instance is None:
            return self

        # Lazy initialization - compute path when first accessed
        if not hasattr(instance, self._attr_name):
            folder = getattr(instance, self.folder_key, "")
            if callable(self.filename_template):
                filename = self.filename_template(instance)
            else:
                filename = self.filename_template.format(
                    tool_name=instance.TOOL_NAME,
                    job_name=getattr(instance, 'name', 'job'),
                    pipeline_name=instance._extract_pipeline_name() if hasattr(instance, '_extract_pipeline_name') else 'pipeline'
                )
            setattr(instance, self._attr_name, os.path.join(folder, filename))

        return getattr(instance, self._attr_name)


# Usage in tools becomes:
class AlphaFold(BaseConfig):
    # Auto-managed file paths - no manual setup needed!
    queries_csv = FilePathDescriptor("{pipeline_name}_queries.csv")
    queries_fasta = FilePathDescriptor("{pipeline_name}_queries.fasta")
    main_datasheet = FilePathDescriptor("alphafold_results.csv")
```

### 2. Input Handler Mixin

**Problem**: Every tool repeats the same input type checking and extraction logic

**Solution**: Create a mixin that provides standardized input handling

```python
class InputHandlerMixin:
    """Mixin providing standardized input handling for all data types."""

    @staticmethod
    def resolve_input(input_param, expected_type='structures'):
        """
        Universal input resolver - handles all input types consistently.

        Args:
            input_param: Any input type (ToolOutput, StandardizedOutput, str, list, dict)
            expected_type: What type of data to extract ('structures', 'sequences', 'compounds')

        Returns:
            ResolvedInput object with .files, .ids, .datasheets attributes
        """
        if isinstance(input_param, ToolOutput):
            return ResolvedInput(
                files=input_param.get_output_files(expected_type),
                ids=input_param.get_output_files(f"{expected_type[:-1]}_ids"),  # structures -> structure_ids
                datasheets=input_param.datasheets,
                source=input_param
            )

        elif isinstance(input_param, StandardizedOutput):
            return ResolvedInput(
                files=getattr(input_param, expected_type, []),
                ids=getattr(input_param, f"{expected_type[:-1]}_ids", []),
                datasheets=input_param.datasheets,
                source=input_param
            )

        elif isinstance(input_param, list):
            return ResolvedInput(
                files=input_param,
                ids=[os.path.splitext(os.path.basename(f))[0] for f in input_param],
                datasheets=None,
                source=input_param
            )

        elif isinstance(input_param, str):
            return ResolvedInput(
                files=[input_param],
                ids=[os.path.splitext(os.path.basename(input_param))[0]],
                datasheets=None,
                source=input_param
            )

        else:
            raise ValueError(f"Unsupported input type: {type(input_param)}")


class ResolvedInput:
    """Container for resolved input data."""
    def __init__(self, files, ids, datasheets, source):
        self.files = files
        self.ids = ids
        self.datasheets = datasheets
        self.source = source
        self.is_tool_output = isinstance(source, ToolOutput)


# Usage in tools:
class ProteinMPNN(InputHandlerMixin, BaseConfig):
    def configure_inputs(self, pipeline_folders):
        # Old way: 50+ lines of if-elif-else
        # New way: 1 line!
        resolved = self.resolve_input(self.input_structures, 'structures')
        self.input_pdb_files = resolved.files
        self.structure_ids = resolved.ids
```

### 3. Datasheet Navigator Mixin

**Problem**: Complex nested logic to access datasheets from various formats

**Solution**: Unified datasheet navigation interface

```python
class DatasheetNavigatorMixin:
    """Mixin providing elegant datasheet navigation."""

    def get_datasheet(self, source, name=None, fallback_names=None):
        """
        Universal datasheet getter - works with all formats.

        Args:
            source: ToolOutput, StandardizedOutput, DatasheetContainer, dict, or list
            name: Preferred datasheet name (e.g., 'structures', 'sequences')
            fallback_names: List of fallback names to try if primary name not found

        Returns:
            DatasheetInfo object or path string
        """
        fallback_names = fallback_names or ['main', 'structures', 'sequences']
        names_to_try = [name] + fallback_names if name else fallback_names

        # Handle DatasheetContainer (modern format)
        if hasattr(source, '_datasheets'):
            for n in names_to_try:
                if n in source._datasheets:
                    return source._datasheets[n]

        # Handle dict format
        if isinstance(source, dict):
            for n in names_to_try:
                if n in source:
                    return source[n]

        # Handle ToolOutput/StandardizedOutput
        if hasattr(source, 'datasheets'):
            return self.get_datasheet(source.datasheets, name, fallback_names)

        # Handle legacy list format
        if isinstance(source, list) and source:
            return source[0]

        raise ValueError(f"Could not find datasheet '{name}' in {type(source)}")


# Usage:
class Filter(DatasheetNavigatorMixin, BaseConfig):
    def configure_inputs(self, pipeline_folders):
        # Old way: 30+ lines of nested if-else checking formats
        # New way: 1 line!
        self.input_datasheet = self.get_datasheet(
            self.data_input,
            name='filtered',
            fallback_names=['structures', 'sequences', 'main']
        )
```

### 4. Script Generator Decorators

**Problem**: Repetitive bash script generation with manual header/footer management

**Solution**: Method decorators for script sections

```python
def script_section(order=0, header=None):
    """Decorator to mark methods as script sections with auto-ordering."""
    def decorator(func):
        func._is_script_section = True
        func._section_order = order
        func._section_header = header
        return func
    return decorator


class BaseConfig(ABC):
    def generate_script(self, script_path: str) -> str:
        """
        Auto-generate script from decorated methods.
        Collects all methods marked with @script_section and assembles them.
        """
        sections = []

        # Collect all script section methods
        for name in dir(self):
            method = getattr(self, name)
            if callable(method) and getattr(method, '_is_script_section', False):
                sections.append(method)

        # Sort by section order
        sections.sort(key=lambda m: m._section_order)

        # Build script
        script_parts = [
            "#!/bin/bash",
            f"# {self.TOOL_NAME} execution script",
            "# Generated by BioPipelines pipeline system\n",
            self.generate_completion_check_header()
        ]

        for section in sections:
            header = section._section_header
            if header:
                script_parts.append(f"\n# {header}")
            script_parts.append(section())

        script_parts.append(self.generate_completion_check_footer())

        return "\n".join(script_parts)


# Usage in tools:
class AlphaFold(BaseConfig):
    @script_section(order=1, header="Prepare input sequences")
    def prepare_sequences(self):
        if "direct_sequence" in self.input_sources:
            return f"""cat > {self.queries_csv} << EOF
id,sequence
{self.name},{self.input_sources["direct_sequence"]}
EOF
"""
        # ... rest of logic

    @script_section(order=2, header="Run AlphaFold")
    def run_alphafold(self):
        return f"""
{self.colabfold_batch} {self.queries_csv} "{self.output_folder}" {self.af_options}
"""
```

### 5. Parameter Validator Decorator

**Problem**: Repetitive parameter validation in every `validate_params()` method

**Solution**: Declarative parameter validation

```python
def validate_param(param_name, validator):
    """Decorator to add parameter validators."""
    def decorator(cls):
        if not hasattr(cls, '_param_validators'):
            cls._param_validators = {}
        cls._param_validators[param_name] = validator
        return cls
    return decorator


class Validators:
    """Common validation functions."""

    @staticmethod
    def positive(value, name):
        if value <= 0:
            raise ValueError(f"{name} must be positive")

    @staticmethod
    def range_check(min_val, max_val):
        def validator(value, name):
            if not min_val <= value <= max_val:
                raise ValueError(f"{name} must be between {min_val} and {max_val}")
        return validator

    @staticmethod
    def not_empty(value, name):
        if not value:
            raise ValueError(f"{name} cannot be empty")


# Usage:
@validate_param('num_sequences', Validators.positive)
@validate_param('plddt_threshold', Validators.range_check(0, 100))
@validate_param('expression', Validators.not_empty)
class ProteinMPNN(BaseConfig):
    def validate_params(self):
        # Auto-validation happens in base class
        # Only add tool-specific complex validation here
        if self.fixed and self.redesigned:
            # Custom validation logic
            pass
```

### 6. Configuration Display Builder

**Problem**: Manual string building for configuration display

**Solution**: Declarative configuration display

```python
class ConfigDisplay:
    """Builder for configuration display."""

    def __init__(self):
        self._lines = []

    def add(self, key, value, formatter=None):
        """Add a config line with optional formatting."""
        if formatter:
            value = formatter(value)
        self._lines.append(f"{key}: {value}")
        return self

    def add_if(self, condition, key, value, formatter=None):
        """Conditionally add a config line."""
        if condition:
            self.add(key, value, formatter)
        return self

    def build(self):
        """Get all lines."""
        return self._lines


# Usage:
class ProteinMPNN(BaseConfig):
    def get_config_display(self):
        return (ConfigDisplay()
                .add("INPUT", self.input_structures)
                .add("NUM SEQUENCES", self.num_sequences)
                .add("FIXED", self.fixed or "None")
                .add("REDESIGNED", self.redesigned or "None")
                .add_if(self.plddt_threshold < 100,
                       "pLDDT THR", self.plddt_threshold)
                .build())
```

## Implementation Plan

### Phase 1: Core Infrastructure (Week 1)
1. Create `PipelineScripts/mixins/` directory
2. Implement base mixins:
   - `InputHandlerMixin`
   - `DatasheetNavigatorMixin`
   - `FilePathManager`
3. Add to `base_config.py`
4. Write comprehensive tests

### Phase 2: Decorator System (Week 2)
1. Implement decorators:
   - `@script_section`
   - `@validate_param`
   - `FilePathDescriptor`
2. Update `BaseConfig` to use decorators
3. Create example refactored tool

### Phase 3: Tool Migration (Week 3-4)
1. Refactor 2-3 tools as pilots:
   - `AlphaFold` (complex)
   - `Filter` (medium)
   - `SelectBest` (simple)
2. Gather feedback
3. Document migration guide
4. Migrate remaining tools

### Phase 4: Documentation & Testing (Week 5)
1. Update user manual
2. Add migration guide for custom tools
3. Comprehensive testing
4. Performance benchmarking

## Expected Benefits

### Code Reduction
- **50-70% reduction** in tool initialization code
- **30-40% reduction** in overall tool file sizes
- **90% reduction** in datasheet navigation code

### Maintainability
- Centralized input handling logic
- Consistent patterns across all tools
- Easier to add new tools
- Better IDE support with descriptors

### Reliability
- Standardized validation
- Fewer edge cases (centralized handling)
- Better error messages
- Type hints throughout

## Example: Before & After

### Before (ProteinMPNN - 627 lines)

```python
class ProteinMPNN(BaseConfig):
    def __init__(self, structures, ...):
        # 50+ lines of input type checking
        if isinstance(structures, ToolOutput):
            self.input_structures = structures
            self.input_is_tool_output = True
        elif isinstance(structures, StandardizedOutput):
            # ... more checks
        elif isinstance(structures, list):
            # ... more checks
        # ... etc

        self._initialize_file_paths()

    def _initialize_file_paths(self):
        self.parsed_pdbs_jsonl = None
        self.fixed_jsonl = None
        # ... 10+ more paths

    def _setup_file_paths(self):
        self.parsed_pdbs_jsonl = os.path.join(self.output_folder, "parsed_pdbs.jsonl")
        # ... 10+ more path setups

    def configure_inputs(self, pipeline_folders):
        # 100+ lines of if-elif-else
        if self.input_is_tool_output:
            tool_output = self.input_structures
            source_pdbs = tool_output.get_output_files("pdbs")
            # ... more extraction
        elif isinstance(self.input_structures, list):
            # ... handle list
        # ... etc
```

### After (ProteinMPNN - ~400 lines)

```python
class ProteinMPNN(InputHandlerMixin, DatasheetNavigatorMixin, BaseConfig):
    # Declarative file paths
    parsed_pdbs_jsonl = FilePathDescriptor("parsed_pdbs.jsonl")
    fixed_jsonl = FilePathDescriptor("fixed_pos.jsonl")
    seqs_folder = FilePathDescriptor("seqs")
    main_datasheet = FilePathDescriptor("proteinmpnn_results.csv")
    queries_csv = FilePathDescriptor("{pipeline_name}_queries.csv")

    @validate_param('num_sequences', Validators.positive)
    @validate_param('plddt_threshold', Validators.range_check(0, 100))
    def __init__(self, structures, ...):
        # Simple parameter storage - no input handling needed
        self.structures_input = structures
        self.num_sequences = num_sequences
        # ... just parameter storage
        super().__init__(**kwargs)

    def configure_inputs(self, pipeline_folders):
        # One line to resolve all input types
        resolved = self.resolve_input(self.structures_input, 'structures')
        self.input_pdb_files = resolved.files
        self.structure_ids = resolved.ids

        # One line to get datasheet regardless of format
        self.input_datasheet = self.get_datasheet(
            resolved.datasheets,
            name='structures',
            fallback_names=['main']
        )

    @script_section(order=1, header="Parse structures")
    def parse_structures(self):
        return f"python {self.parse_py} --input {self.input_directory}"

    @script_section(order=2, header="Run ProteinMPNN")
    def run_proteinmpnn(self):
        return f"python {self.pmpnn_py} --num_seq {self.num_sequences}"
```

**Result**: 36% fewer lines, much more readable, centralized logic

## Backward Compatibility

All changes will maintain backward compatibility:
- Existing tools continue to work unchanged
- Migration is opt-in
- Mixins don't override existing methods
- Decorators are additive

## Risks & Mitigation

| Risk | Mitigation |
|------|-----------|
| Breaking existing tools | Comprehensive testing, gradual rollout |
| Performance overhead | Benchmark decorators, use lazy evaluation |
| Learning curve | Detailed documentation, examples |
| Over-abstraction | Keep it simple, only abstract truly common patterns |

## Conclusion

This refactoring will significantly improve code quality, reduce duplication, and make the BioPipelines framework more elegant and maintainable while preserving all existing functionality.
