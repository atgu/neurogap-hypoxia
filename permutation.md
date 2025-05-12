```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.stats import rankdata


def calculate_ratio(InBed, Expect):
    try:
        InBed = float(InBed)
        Expect = float(Expect)
        if np.isnan(InBed) or np.isnan(Expect) or Expect == 0:
            return 0.0
        return InBed / Expect
    except:
        return 0.0

def read_data(file_path, observed_bed, random_bed_prefix='regions_'):
    """Read data and compute ratios + extract P-values (NA → 1)"""
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
    
    # Extract observed P-value (NA → 1)
    observed_p = observed_row['PValue'].values[0]
    observed_p = 1.0 if pd.isna(observed_p) else observed_p
    
    # Compute permutation ratios and P-values
    perm_rows = data[data['Bed_File'].str.startswith(random_bed_prefix)]
    perm_ratios = np.array([
        calculate_ratio(row['InBed_Index_SNP'], row['ExpectNum_of_InBed_SNP'])
        for _, row in perm_rows.iterrows()
    ])
    perm_ps = perm_rows['PValue'].fillna(1.0).values  # NA → 1
    
    return observed_ratio, perm_ratios, observed_p, perm_ps


def calculate_fdr_adjusted_pvalues(p_values):
    m = len(p_values)
    if m == 0:
        return np.array([])
    
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    
    adjusted_p = np.zeros(m)
    for i in range(m):
        adjusted_p[i] = sorted_p[i] * m / (i + 1)  # rank = i + 1
    
    for i in range(m-2, -1, -1):
        adjusted_p[i] = min(adjusted_p[i], adjusted_p[i+1])

    adjusted_p = np.minimum(adjusted_p, 1.0)
    
    final_adjusted_p = np.zeros(m)
    final_adjusted_p[sorted_indices] = adjusted_p
    
    return final_adjusted_p

def calculate_significance(observed_p, perm_ps, observed_ratio, perm_ratios):
    """Calculate statistical significance using P-values"""
    results = {}
    results['observed_p'] = observed_p
    
    # Empirical p-value
    results['empirical_p'] = np.mean(perm_ps <= observed_p)
    
    
    # FDR correction (Benjamini-Hochberg)
    all_p = np.concatenate([[observed_p], perm_ps])  
    adjusted_p = calculate_fdr_adjusted_pvalues(all_p)
    results['fdr_adjusted_p'] = adjusted_p[0]  
    
    return results

def plot_results(observed_ratio, perm_ratios, results, observed_bed, trait):
    """Plot ratio distribution with P-value annotations"""
    plt.figure(figsize=(12, 8))

    plt.hist(perm_ratios, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
    
    # Mark observed ratio (red dashed)
    plt.axvline(observed_ratio, color='red', linestyle='--', linewidth=2,
                label=f'Observed ratio: {observed_ratio:.3f}')
    
    # Mark mean (green solid)
    plt.axvline(np.mean(perm_ratios), color='green', linestyle='-', linewidth=1,
                label=f'Mean: {np.mean(perm_ratios):.3f}')
    
    # Add stats box (upper right)
    stats_text = (
        f"Observed P: {results['observed_p']:.3e}\n"
        f"FDR-adjusted P: {results['fdr_adjusted_p']:.3e}\n"
        f"Empirical P: {results['empirical_p']:.3e}" 
    )
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
             ha='right', va='top', bbox=dict(facecolor='white', alpha=0.8))

    
    plt.xlabel('Ratio (InBed_Index_SNP / ExpectNum_of_InBed_SNP)')
    plt.ylabel('Frequency')
    plt.title(f'Ratio Distribution - {observed_bed} ({trait})')
    plt.legend(loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Save plot
    output_png = f"Ethio_Kemri_{trait}_ratio.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    return output_png

def main():
    parser = argparse.ArgumentParser(description='Ratio analysis with P-value stats')
    parser.add_argument('--file_path', required=True, help='Input TSV file path')
    parser.add_argument('--observed_bed', required=True, help='Observed BED filename')
    parser.add_argument('--trait', required=True, help='Trait name (for output naming)')
    
    args = parser.parse_args()
    
    # Read data and compute ratios + P-values
    observed_ratio, perm_ratios, observed_p, perm_ps = read_data(args.file_path, args.observed_bed)
    
    # Calculate significance using P-values
    results = calculate_significance(observed_p, perm_ps, observed_ratio, perm_ratios)
    
    # Plot ratio distribution with P-value annotations
    plot_file = plot_results(observed_ratio, perm_ratios, results, args.observed_bed, args.trait)
    
    # Print results
    print(f"\n=== Ratio Analysis Results ===")
    print(f"Trait: {args.trait}")
    print(f"Observed ratio: {observed_ratio:.4f}")
    print(f"Observed P-value: {observed_p:.3e}")
    print(f"Permutation mean ratio: {np.mean(perm_ratios):.4f}")
    print(f"\nPlot saved to: {plot_file}")

if __name__ == "__main__":
    main()
```
