```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

# -----------------------------
# Compute enrichment ratio
# -----------------------------
def calculate_ratio(InBed, Expect):
    try:
        InBed = float(InBed)
        Expect = float(Expect)
        # Return 0 if InBed or Expect is NaN or Expect is 0
        if np.isnan(InBed) or np.isnan(Expect) or Expect == 0:
            return 0.0
        return InBed / Expect
    except:
        return 0.0

# -----------------------------
# Read TSV file and extract observed & permutation values
# -----------------------------
def read_data(file_path, observed_bed, random_bed_prefix='regions_'):
    data = pd.read_csv(file_path, sep="\t")
    
    # Extract observed row
    observed_row = data[data['Bed_File'] == observed_bed]
    if observed_row.empty:
        raise ValueError(f"Observed bed file {observed_bed} not found")
    
    # Calculate observed ratio
    observed_ratio = calculate_ratio(
        observed_row['InBed_Index_SNP'].values[0],
        observed_row['ExpectNum_of_InBed_SNP'].values[0]
    )
    
    # Handle NA p-values: convert to 1.0
    observed_p = observed_row['PValue'].values[0]
    observed_p = 1.0 if pd.isna(observed_p) else observed_p
    
    # Extract permutations and compute ratios
    perm_rows = data[data['Bed_File'].str.startswith(random_bed_prefix)]
    perm_ratios = np.array([
        calculate_ratio(row['InBed_Index_SNP'], row['ExpectNum_of_InBed_SNP'])
        for _, row in perm_rows.iterrows()
    ])
    perm_ps = perm_rows['PValue'].fillna(1.0).values  # Replace NA with 1
    
    return observed_ratio, perm_ratios, observed_p, perm_ps

# -----------------------------
# Significance testing
# -----------------------------
def calculate_significance(observed_p, perm_ps, observed_ratio, perm_ratios):
    results = {}
    results['observed_p'] = observed_p
    
    # Empirical p: fraction of permutations with p <= observed p
    results['empirical_p'] = np.mean(perm_ps <= observed_p)
    
    # FDR correction
    all_p = np.concatenate([[observed_p], perm_ps])
    rejected, qvals, _, _ = multipletests(all_p, alpha=0.05, method='fdr_bh')
    results['fdr_adjusted_p'] = qvals[0]  # q-value for observed_p
    
    return results

# -----------------------------
# Plot ratio histogram and annotate stats
# -----------------------------
def plot_results(observed_ratio, perm_ratios, results, observed_bed, trait):
    plt.figure(figsize=(12, 8))

    # Plot permutation ratio distribution
    plt.hist(perm_ratios, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
    
    # Mark observed ratio
    plt.axvline(observed_ratio, color='red', linestyle='--', linewidth=2,
                label=f'Observed ratio: {observed_ratio:.3f}')
    
    # Mark permutation mean
    plt.axvline(np.mean(perm_ratios), color='green', linestyle='-', linewidth=1,
                label=f'Mean: {np.mean(perm_ratios):.3f}')
    
    # Annotate p-values
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
    
    # Save figure
    output_png = f"Ethio_Kemri_{trait}_ratio.png"
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    return output_png

# -----------------------------
# Main function
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description='Ratio analysis with P-value stats')
    parser.add_argument('--file_path', required=True, help='Input TSV file path')
    parser.add_argument('--observed_bed', required=True, help='Observed BED filename')
    parser.add_argument('--trait', required=True, help='Trait name (for output naming)')
    
    args = parser.parse_args()
    
    # Load data and compute values
    observed_ratio, perm_ratios, observed_p, perm_ps = read_data(args.file_path, args.observed_bed)
    
    # Perform significance testing
    results = calculate_significance(observed_p, perm_ps, observed_ratio, perm_ratios)
    
    # Plot result and save figure
    plot_file = plot_results(observed_ratio, perm_ratios, results, args.observed_bed, args.trait)
    
    # Print result summary
    print(f"\n=== Ratio Analysis Results ===")
    print(f"Trait: {args.trait}")
    print(f"Observed ratio: {observed_ratio:.4f}")
    print(f"Observed P-value: {observed_p:.3e}")
    print(f"Permutation mean ratio: {np.mean(perm_ratios):.4f}")
    print(f"\nPlot saved to: {plot_file}")

# -----------------------------
# Run script
# -----------------------------
if __name__ == "__main__":
    main()
```
