#!/usr/bin/env python3
"""
Runtime filter execution helper script for new filter architecture.

This script is called by filter bash scripts to perform the actual filtering
logic at pipeline execution time using the new Filter/Criterion architecture.
"""

import os
import sys
import argparse
import json
import pandas as pd
import importlib
from typing import Dict, List, Any, Optional
from pathlib import Path


def load_criterion_results(results_file: str) -> Dict[str, Any]:
    """
    Load criterion evaluation results from JSON file.
    
    Args:
        results_file: Path to criterion results JSON file
        
    Returns:
        Dictionary with criterion results
    """
    try:
        with open(results_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading criterion results from {results_file}: {e}")
        raise


def combine_criterion_results(criterion_results: List[Dict[str, Any]], 
                                combination: str,
                                score_weights: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """
    Combine multiple criterion results using the specified combination method.
    
    Args:
        criterion_results: List of individual criterion results
        combination: Combination method ("AND", "OR", "WEIGHTED")
        score_weights: Weights for WEIGHTED combination
        
    Returns:
        Combined filter results
    """
    if not criterion_results:
        return {
            "kept_items": [],
            "filtered_items": [],
            "item_scores": {},
            "total_input": 0,
            "kept_count": 0,
            "filtered_count": 0,
            "pass_rate": 0.0,
            "combination_method": combination,
            "individual_results": []
        }
    
    # Get all items that were evaluated
    all_items = set()
    for result in criterion_results:
        all_items.update(result["kept_items"])
        all_items.update(result["filtered_items"])
    
    all_items = list(all_items)
    
    if combination == "AND":
        # Items must pass ALL criteria
        kept_items = []
        for item in all_items:
            passes_all = True
            for result in criterion_results:
                if item not in result["kept_items"]:
                    passes_all = False
                    break
            if passes_all:
                kept_items.append(item)
        
        filtered_items = [item for item in all_items if item not in kept_items]
        final_scores = {item: 1.0 if item in kept_items else 0.0 for item in all_items}
    
    elif combination == "OR":
        # Items pass if they satisfy ANY criterion
        kept_items = []
        for item in all_items:
            passes_any = False
            for result in criterion_results:
                if item in result["kept_items"]:
                    passes_any = True
                    break
            if passes_any:
                kept_items.append(item)
        
        filtered_items = [item for item in all_items if item not in kept_items]
        final_scores = {item: 1.0 if item in kept_items else 0.0 for item in all_items}
    
    elif combination == "WEIGHTED":
        # Weighted combination of scores
        if not score_weights:
            # Default equal weights
            score_weights = {result["criterion_class"]: 1.0/len(criterion_results) 
                           for result in criterion_results}
        
        # Calculate weighted scores
        final_scores = {}
        for item in all_items:
            weighted_score = 0.0
            total_weight = 0.0
            
            for result in criterion_results:
                criterion_class = result["criterion_class"]
                weight = score_weights.get(criterion_class, 0.0)
                
                if item in result["item_scores"]:
                    # Normalize score to 0-1 based on whether it passed the expression
                    raw_score = result["item_scores"][item]
                    normalized_score = 1.0 if item in result["kept_items"] else 0.0
                    weighted_score += weight * normalized_score
                    total_weight += weight
            
            if total_weight > 0:
                final_scores[item] = weighted_score / total_weight
            else:
                final_scores[item] = 0.0
        
        # Sort by score and keep items with score > threshold
        score_threshold = 0.5  # Could be configurable
        kept_items = [item for item in all_items if final_scores[item] > score_threshold]
        filtered_items = [item for item in all_items if item not in kept_items]
    
    else:
        raise ValueError(f"Unknown combination method: {combination}")
    
    return {
        "kept_items": kept_items,
        "filtered_items": filtered_items,
        "item_scores": final_scores,
        "total_input": len(all_items),
        "kept_count": len(kept_items),
        "filtered_count": len(filtered_items),
        "pass_rate": len(kept_items) / len(all_items) if all_items else 0.0,
        "combination_method": combination,
        "individual_results": criterion_results
    }


def apply_max_items_filter(combined_results: Dict[str, Any], max_items: Optional[int]) -> Dict[str, Any]:
    """
    Apply max_items filter to combined results.
    
    Args:
        combined_results: Combined filter results
        max_items: Maximum number of items to keep
        
    Returns:
        Updated results with max_items applied
    """
    if max_items is None or len(combined_results["kept_items"]) <= max_items:
        return combined_results
    
    # Sort kept items by score (highest first)
    kept_items_with_scores = [(item, combined_results["item_scores"][item]) 
                            for item in combined_results["kept_items"]]
    kept_items_with_scores.sort(key=lambda x: x[1], reverse=True)
    
    # Keep only top max_items
    final_kept = [item for item, score in kept_items_with_scores[:max_items]]
    additional_filtered = [item for item, score in kept_items_with_scores[max_items:]]
    
    # Update results
    combined_results["kept_items"] = final_kept
    combined_results["filtered_items"].extend(additional_filtered)
    combined_results["kept_count"] = len(final_kept)
    combined_results["filtered_count"] = len(combined_results["filtered_items"])
    combined_results["pass_rate"] = len(final_kept) / combined_results["total_input"] if combined_results["total_input"] > 0 else 0.0
    
    return combined_results


def save_filtering_results(result: Dict[str, Any], output_folder: str, job_name: str, 
                          filter_type: str, input_data: Dict[str, Any]):
    """
    Save filtering results to files.
    
    Args:
        result: Combined filter results dictionary
        output_folder: Output directory
        job_name: Job name for file naming
        filter_type: Type of data that was filtered
        input_data: Original input data for reference
    """
    # Save filter manifest
    manifest_file = os.path.join(output_folder, f"{job_name}_filter_manifest.json")
    with open(manifest_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    # Create filtered datasheet
    datasheet_file = os.path.join(output_folder, f"{job_name}_filtered_{filter_type}.csv")
    
    # Build datasheet based on filter type
    rows = []
    
    # Add kept items
    for item in result["kept_items"]:
        item_id = os.path.splitext(os.path.basename(item))[0] if os.path.sep in item else item
        rows.append({
            "id": item_id,
            "file_path": item,
            "filter_passed": True,
            "final_score": result["item_scores"].get(item, 1.0),
            "combination_method": result["combination_method"]
        })
    
    # Add filtered items
    for item in result["filtered_items"]:
        item_id = os.path.splitext(os.path.basename(item))[0] if os.path.sep in item else item
        rows.append({
            "id": item_id,
            "file_path": item,
            "filter_passed": False,
            "final_score": result["item_scores"].get(item, 0.0),
            "combination_method": result["combination_method"]
        })
    
    # Save as CSV
    df = pd.DataFrame(rows)
    df.to_csv(datasheet_file, index=False)
    
    # Generate filter report
    report_file = os.path.join(output_folder, f"{job_name}_filter_report.txt")
    with open(report_file, 'w') as f:
        f.write(f"Filter Results Report\n")
        f.write(f"====================\n\n")
        f.write(f"Job: {job_name}\n")
        f.write(f"Filter Type: {filter_type}\n")
        f.write(f"Combination Method: {result['combination_method']}\n")
        f.write(f"Criteria Count: {len(result['individual_results'])}\n\n")
        
        f.write(f"Overall Results:\n")
        f.write(f"  Total Input: {result['total_input']}\n")
        f.write(f"  Items Kept: {result['kept_count']}\n")
        f.write(f"  Items Filtered: {result['filtered_count']}\n")
        f.write(f"  Pass Rate: {result['pass_rate']:.1%}\n\n")
        
        # Individual criterion results
        if result["individual_results"]:
            f.write("Individual Criterion Results:\n")
            for i, individual_result in enumerate(result["individual_results"], 1):
                criterion_type = individual_result.get('criterion_type', 'unknown')
                criterion_class = individual_result.get('criterion_class', 'unknown')
                f.write(f"  {i}. {criterion_class} ({criterion_type}):\n")
                f.write(f"     Expression: {individual_result.get('expression', 'N/A')}\n")
                f.write(f"     Kept: {individual_result['kept_count']}/{individual_result['total_input']} ({individual_result['pass_rate']:.1%})\n")
            f.write("\n")
        
        if result["kept_items"]:
            f.write("Kept Items:\n")
            for item in result["kept_items"]:
                score = result["item_scores"].get(item, "N/A")
                f.write(f"  - {os.path.basename(item)} (score: {score})\n")
            f.write("\n")
        
        if result["filtered_items"] and len(result["filtered_items"]) <= 20:  # Don't list too many
            f.write("Filtered Items:\n")
            for item in result["filtered_items"]:
                score = result["item_scores"].get(item, "N/A")
                f.write(f"  - {os.path.basename(item)} (score: {score})\n")
        elif result["filtered_items"]:
            f.write(f"Filtered Items: {len(result['filtered_items'])} items (too many to list)\n")
    
    print(f"Filter results saved:")
    print(f"  Manifest: {manifest_file}")
    print(f"  Datasheet: {datasheet_file}")
    print(f"  Report: {report_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Runtime filter execution for new filter architecture"
    )
    
    parser.add_argument(
        "--input", 
        required=True,
        help="JSON string with standardized input data"
    )
    
    parser.add_argument(
        "--filter-type",
        required=True,
        choices=["structures", "sequences", "compounds"],
        help="Type of data to filter"
    )
    
    parser.add_argument(
        "--output-folder",
        required=True,
        help="Output directory for filter results"
    )
    
    parser.add_argument(
        "--job-name",
        required=True,
        help="Job name for file naming"
    )
    
    parser.add_argument(
        "--combination",
        required=True,
        choices=["AND", "OR", "WEIGHTED"],
        help="Combination method for multiple criteria"
    )
    
    parser.add_argument(
        "--criteria",
        required=True,
        help="JSON string with paths to criterion result files"
    )
    
    parser.add_argument(
        "--score-weights",
        help="JSON string with score weights for WEIGHTED combination"
    )
    
    parser.add_argument(
        "--max-items",
        type=int,
        help="Maximum number of items to keep"
    )
    
    parser.add_argument(
        "--filter-class",
        default="Filter",
        help="Filter class to use (Filter)"
    )
    
    args = parser.parse_args()
    
    try:
        # Parse input data
        input_data = json.loads(args.input)
        criteria_results_paths = json.loads(args.criteria)  # Now expects paths to result files
        
        score_weights = None
        if args.score_weights:
            score_weights = json.loads(args.score_weights)
        
        print(f"Starting combined filtering...")
        print(f"Filter type: {args.filter_type}")
        print(f"Combination: {args.combination}")
        print(f"Criteria results: {len(criteria_results_paths)}")
        print(f"Output folder: {args.output_folder}")
        print(f"Max items: {args.max_items}")
        
        # Load criterion results
        criterion_results = []
        for results_path in criteria_results_paths:
            print(f"Loading results from: {results_path}")
            result = load_criterion_results(results_path)
            criterion_results.append(result)
            print(f"  {result['criterion_class']}: {result['kept_count']}/{result['total_input']} kept")
        
        if not criterion_results:
            print("No criterion results to combine")
            sys.exit(1)
        
        # Combine criterion results
        print(f"Combining results using {args.combination} method...")
        combined_result = combine_criterion_results(criterion_results, args.combination, score_weights)
        
        # Apply max_items filter if specified
        if args.max_items:
            print(f"Applying max_items filter: {args.max_items}")
            combined_result = apply_max_items_filter(combined_result, args.max_items)
        
        print(f"Combined filtering complete: {combined_result['kept_count']}/{combined_result['total_input']} kept ({combined_result['pass_rate']:.1%})")
        
        # Save results
        os.makedirs(args.output_folder, exist_ok=True)
        save_filtering_results(
            combined_result, args.output_folder, args.job_name, 
            args.filter_type, input_data
        )
        
        print("Combined filter execution completed successfully")
        
    except Exception as e:
        print(f"Error during filtering: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()