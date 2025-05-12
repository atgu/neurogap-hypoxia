```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.stats import rankdata

def calculate_ratio(InBed, Expect):
    """Calculate ratio, return 0 if InBed is 0"""
    return 0.0 if InBed == 0 else InBed / Expect

def read_data(file_path, observed_bed, random_bed_prefix='regions_'):
    """Read data and compute observed/permutation ratios"""
    data = pd.read_csv(file_path, sep="\t")
    
    # Extract observed values
    observed_row = data[data['Bed_File'] == observed_bed]
    if observed_row.empty:
        raise ValueError(f"Observed bed file {observed_bed} not found")
    
    # Compute observed ratio
    observed_ratio = calculate_ratio(
        observed_row['InBed_Index_SNP'].values[0],
        observed_row['ExpectNum_of_InBed_SNP'].values[0]
    )
    
    # Compute permutation ratios
    perm_rows = data[data['Bed_File'].str.startswith(random_bed_prefix)]
    perm_ratios = np.array([
        calculate_ratio(row['InBed_Index_SNP'], row['ExpectNum_of_InBed_SNP'])
        for _, row in perm_rows.iterrows()
    ])
    
    return observed_ratio, perm_ratios

def calculate_significance(observed_ratio, perm_ratios):
    """Calculate statistical significance"""
    results = {}
    
    # Empirical p-value (proportion of permutations ≥ observed)
    results['empirical_p'] = np.mean(perm_ratios >= observed_ratio)
    
    # Bonferroni correction
    results['bonferroni_threshold'] = 0.05 / len(perm_ratios)
    results['bonferroni_sig'] = results['empirical_p'] < results['bonferroni_threshold']
    
    # FDR correction (Benjamini-Hochberg)
    all_ratios = np.concatenate([[observed_ratio], perm_ratios])
    ranks = rankdata(-all_ratios)  # Descending rank
    pseudo_p = ranks / (len(all_ratios) + 1)
    
    # Compute FDR threshold
    sorted_p = np.sort(pseudo_p[1:])
    fdr_thresholds = 0.05 * np.arange(1, len(perm_ratios)+1) / len(perm_ratios)
    results['fdr_threshold'] = sorted_p[np.max(np.where(sorted_p <= fdr_thresholds)[0])] if any(sorted_p <= fdr_thresholds) else 0
    results['fdr_sig'] = pseudo_p[0] <= results['fdr_threshold']
    
    return results

def plot_results(observed_ratio, perm_ratios, results, observed_bed, trait):
    """Generate visualization"""
    plt.figure(figsize=(12, 8))
    
    # Plot permutation distribution
    plt.hist(perm_ratios, bins=50, color='skyblue', edgecolor='black', alpha=0.7,
             label=f'Permutation ratios (n={len(perm_ratios)})')
    
    # Mark observed value
    plt.axvline(observed_ratio, color='red', linestyle='--', linewidth=2,
                label=f'Observed: {observed_ratio:.3f}')
    
    # Mark statistics
    plt.axvline(np.mean(perm_ratios), color='green', linestyle='-', linewidth=1,
                label=f'Mean: {np.mean(perm_ratios):.3f}')
    ci_low, ci_high = np.percentile(perm_ratios, [2.5, 97.5])

    # Add stats box
    stats_text = (
        f"Trait: {trait}\n"
        f"Observed ratio: {observed_ratio:.3f}\n"
        f"Empirical p: {results['empirical_p']:.3e}\n"
        f"Bonferroni sig.: {'YES' if results['bonferroni_sig'] else 'NO'}\n"
        f"FDR sig.: {'YES' if results['fdr_sig'] else 'NO'}"
    )
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
             ha='right', va='top', bbox=dict(facecolor='white', alpha=0.8))
    
    plt.xlabel('Ratio (InBed_Index_SNP/ExpectNum_of_InBed_SNP)')
    plt.ylabel('Frequency')
    plt.title(f'Ratio Analysis - {observed_bed} ({trait})')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Save plot
    output_png = f"Ethio_Kemri_{trait}_ratio.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    return output_png

def main():
    parser = argparse.ArgumentParser(description='Ratio analysis tool')
    parser.add_argument('--file_path', required=True, help='Input TSV file path')
    parser.add_argument('--observed_bed', required=True, help='Observed BED filename')
    parser.add_argument('--trait', required=True, help='Trait name (for output naming)')
    
    args = parser.parse_args()
    
    # Run analysis
    observed_ratio, perm_ratios = read_data(args.file_path, args.observed_bed)
    results = calculate_significance(observed_ratio, perm_ratios)
    plot_file = plot_results(observed_ratio, perm_ratios, results, args.observed_bed, args.trait)
    
    # Print results
    print(f"\n=== Ratio Analysis Results ===")
    print(f"Trait: {args.trait}")
    print(f"Observed ratio: {observed_ratio:.4f}")
    print(f"Permutation mean: {np.mean(perm_ratios):.4f}")
    print(f"95% CI: [{np.percentile(perm_ratios, 2.5):.4f}, {np.percentile(perm_ratios, 97.5):.4f}]")
    print(f"\nSignificance testing:")
    print(f"Empirical p-value: {results['empirical_p']:.3e}")
    print(f"Bonferroni threshold: {results['bonferroni_threshold']:.3e} (significant: {'YES' if results['bonferroni_sig'] else 'NO'})")
    print(f"FDR threshold: {results['fdr_threshold']:.3e} (significant: {'YES' if results['fdr_sig'] else 'NO'})")
    print(f"\nPlot saved to: {plot_file}")

if __name__ == "__main__":
    main()
```
