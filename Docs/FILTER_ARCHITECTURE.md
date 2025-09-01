# Filter Architecture Documentation

## Overview

The filter architecture provides a robust solution for filtering structures, sequences, and compounds in protein design pipelines while maintaining pipeline predictability and proper metadata tracking.

## Key Components

### 1. Base StructureFilter Class (`PipelineScripts/structure_filter.py`)

The abstract base class that all filters inherit from. Provides:
- Standardized interface with `apply_filter()` method
- Pipeline integration via `BaseConfig` inheritance
- Metadata tracking with `FilteringResult` objects
- Support for max_items limitation
- Filter-aware completion checking integration

### 2. Composite Filter Class (`PipelineScripts/filter.py`)

Combines multiple individual filters with flexible logic:
- **AND logic**: Items must pass all filters
- **OR logic**: Items pass if they satisfy any filter  
- **WEIGHTED logic**: Combines scores from multiple filters
- Score aggregation and ranking
- Nested filter composition

### 3. Individual Filter Implementations

#### AtomResidueDistanceFilter (`PipelineScripts/atom_residue_distance_filter.py`)
- Filters based on distances between atoms and residues
- Supports ligand-protein binding analysis
- Flexible selection syntax (`ligand.Cl`, `protein.D in TRGDTGH`)
- Expression-based filtering (`value<=5.0`, `value>2.5 and value<8.0`)

#### pLDDTFilter (`PipelineScripts/plddt_filter.py`)  
- Filters based on AlphaFold confidence scores
- Supports region-specific analysis (designed residues, binding sites)
- Multiple aggregation methods (mean, min, max, median)
- Datasheet reference integration

### 4. Enhanced Completion Checking (`HelpScripts/pipe_check_completion.py`)

Filter-aware validation that distinguishes between:
- **Critical files**: Datasheets and manifests (must exist)
- **Content files**: Structures/sequences (can be partially missing with warnings)
- **Filter manifests**: Automatic detection and reporting

### 5. Runtime Execution (`HelpScripts/pipe_structure_filter.py`)

Helper script for filter execution in bash scripts:
- Dynamic filter class loading
- Parameter serialization/deserialization
- Result file generation (manifests, datasheets, reports)

## Architecture Benefits

### 1. Maintains Pipeline Predictability
- Filters declare expected output structure at pipeline construction time
- `get_output_files()` provides file paths even before execution
- Downstream tools can be configured based on expected outputs

### 2. Handles Filtering Uncertainty
- Filter manifests track what was kept vs filtered
- Completion checking accounts for expected vs actual file counts
- Clear distinction between critical (must exist) and content files

### 3. Rich Metadata Tracking
- `FilteringResult` objects capture filtering criteria, scores, and statistics
- `StandardizedOutput` enhanced with filter metadata support
- Pass rates, item counts, and filtering history preserved

### 4. Flexible Composition
- Individual filters can be combined with different logic
- Score weighting allows sophisticated ranking
- Nested composition for complex filtering strategies

## Usage Examples

### Basic Individual Filter
```python
# Filter structures based on pLDDT confidence
quality_filter = pipeline.add(
    pLDDTFilter(
        expression='value>90',
        input=boltz.output,
        score_metric='mean'
    ),
    env="ProteinEnv"
)
```

### Composite Filter with Multiple Criteria
```python
# Combined distance and quality filtering
structure_filter = pipeline.add(
    Filter(
        AtomResidueDistanceFilter(
            atom='ligand.Cl',
            residue='protein.D in TRGDTGH', 
            expression='value<=5'
        ),
        pLDDTFilter(
            expression='value>90',
            selection="input.datasheets.sequences.designed_residues"
        ),
        input=boltz.output,
        combination_method="AND",
        max_items=5
    ),
    env="ProteinEnv"
)
```

### Weighted Scoring Filter
```python
# Combine multiple quality metrics with weights
quality_filter = pipeline.add(
    Filter(
        pLDDTFilter(expression='value>70'),
        BindingAffinityFilter(expression='value<-8.0'),
        input=predictions.output,
        combination_method="WEIGHTED",
        score_weights={
            "pLDDTFilter": 0.7,
            "BindingAffinityFilter": 0.3
        },
        max_items=10
    ),
    env="ProteinEnv"
)
```

## File Structure

### Filter Outputs
Each filter creates:
- `{job_name}_filter_manifest.json`: Complete filtering metadata
- `{job_name}_filtered_{type}.csv`: Datasheet with pass/fail results
- `{job_name}_filter_report.txt`: Human-readable summary

### Manifest Format
```json
{
  "kept_items": ["struct1.pdb", "struct2.pdb"],
  "filtered_items": ["struct3.pdb", "struct4.pdb"],
  "filter_criteria": {
    "filter_type": "plddt",
    "expression": "value>90",
    "score_metric": "mean"
  },
  "filter_scores": {
    "struct1.pdb": 92.5,
    "struct2.pdb": 88.1,
    "struct3.pdb": 75.3,
    "struct4.pdb": 45.2
  },
  "total_input": 4,
  "kept_count": 2,
  "filtered_count": 2,
  "pass_rate": 0.5
}
```

## Integration with Existing Tools

### Pipeline Integration
- Filters work seamlessly with existing pipeline tools
- Use standardized input/output format
- Support all existing pipeline features (environments, resources, etc.)

### Downstream Tool Compatibility  
- Tools automatically handle missing files gracefully
- Filter metadata available to downstream tools
- Completion checking provides warnings but doesn't fail pipeline

### Enhanced StandardizedOutput
```python
# Access filter information
if tool_output.output.is_filtered:
    print(f"Pass rate: {tool_output.output.get_filter_pass_rate():.1%}")
    print(f"Original count: {tool_output.output.get_original_items_count()}")
    print(f"Kept count: {tool_output.output.get_kept_items_count()}")
```

## Adding New Filter Types

To add a new filter:

1. **Inherit from StructureFilter**:
```python
class MyFilter(StructureFilter):
    def __init__(self, my_param, **kwargs):
        self.my_param = my_param
        super().__init__(**kwargs)
```

2. **Implement apply_filter()**:
```python
def apply_filter(self, items):
    kept_items = []
    filtered_items = []
    scores = {}
    
    # Your filtering logic here
    
    return FilteringResult(
        kept_items=kept_items,
        filtered_items=filtered_items,
        filter_criteria={"my_param": self.my_param},
        filter_scores=scores
    )
```

3. **Add to runtime helper**:
Update `HelpScripts/pipe_structure_filter.py` to recognize your new filter class.

## Testing

Run the integration test suite:
```bash
python test_filter_integration.py
```

This validates:
- Basic filter functionality
- Composite filter logic
- Output integration
- Completion checking
- Metadata tracking

The filter architecture is now fully integrated and ready for production use in your protein design pipelines!