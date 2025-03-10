### before annotating ancestral allele
```python
import hail as hl
hl.init()

ht_path = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"
ht = hl.read_table(ht_path)
eur_samples = ht.filter(ht.hgdp_tgp_meta.genetic_region == 'EUR')
eur_sample_list = eur_samples.s.collect()
eur_samples_set = hl.literal(set(eur_sample_list))

vcf_path = f"gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chrom}.vcf.bgz"
mt = hl.import_vcf(vcf_path, reference_genome='GRCh38')
mt_eur = mt.filter_cols(eur_samples_set.contains(mt.s))

mt_eur = mt_eur.drop(mt_eur.info)
mt_eur = mt_eur.annotate_rows(
    AC=hl.agg.sum(mt_eur.GT.n_alt_alleles()),
    AN=hl.agg.count_where(hl.is_defined(mt_eur.GT)) * 2
)
mt_eur = mt_eur.annotate_rows(AF=mt_eur.AC / mt_eur.AN)
mt_eur = mt_eur.annotate_rows(hwe=hl.agg.hardy_weinberg_test(mt_eur.GT))

# annotate testsnp - biallelic variant, hwe P > 1e-6, AF(0.05-0.95)
mt_eur = mt_eur.annotate_rows(
    testsnp=(
        (mt_eur.hwe.p_value > 1e-6) &
        (mt_eur.AF >= 0.05) & (mt_eur.AF <= 0.95) &
        (hl.len(mt_eur.alleles) == 2)
    )
)

# annotate singletons
mt_eur = mt_eur.annotate_rows(
    singleton=(
        (mt_eur.AC == 1) &
        (hl.len(mt_eur.alleles[0]) == 1) & (hl.len(mt_eur.alleles[1]) == 1)
    )
)

mt_eur_testsnp = mt_eur.filter_rows(mt_eur.testsnp).select_rows('rsid', 'filters').select_entries('GT')
mt_eur_singletons = mt_eur.filter_rows(mt_eur.singleton).select_rows().select_entries()

# extract singleton positions
singleton_positions = hl.literal(
    mt_eur_singletons.aggregate_rows(hl.agg.collect_as_set(mt_eur_singletons.locus.position))
)

# calculate upstream and downstream singleton counts
margin = 1000000
mt_eur_testsnp = mt_eur_testsnp.annotate_rows(
    upstream_start=mt_eur_testsnp.locus.position - margin,
    upstream_end=mt_eur_testsnp.locus.position,
    downstream_start=mt_eur_testsnp.locus.position,
    downstream_end=mt_eur_testsnp.locus.position + margin
)

mt_eur_testsnp = mt_eur_testsnp.annotate_rows(
    num_singleton_upstream=hl.len(hl.filter(lambda pos: (mt_eur_testsnp.upstream_start <= pos) & (pos < mt_eur_testsnp.upstream_end), singleton_positions)),
    num_singleton_downstream=hl.len(hl.filter(lambda pos: (mt_eur_testsnp.downstream_start <= pos) & (pos < mt_eur_testsnp.downstream_end), singleton_positions))
)

mt_eur_testsnp = mt_eur_testsnp.annotate_rows(
    diff_ratio=hl.if_else(
        hl.min(mt_eur_testsnp.num_singleton_upstream, mt_eur_testsnp.num_singleton_downstream) > 0,
        hl.abs(mt_eur_testsnp.num_singleton_upstream - mt_eur_testsnp.num_singleton_downstream) / hl.max(mt_eur_testsnp.num_singleton_upstream, mt_eur_testsnp.num_singleton_downstream),
        0
    )
)

# extract testsnps with the difference in upstream and downstream singleton counts below 15% 
mt_eur_testsnp = mt_eur_testsnp.filter_rows(
    (mt_eur_testsnp.diff_ratio < 0.15) &
    (mt_eur_testsnp.num_singleton_upstream > 0) &
    (mt_eur_testsnp.num_singleton_downstream > 0)
)

#output
output_path = f"gs://hgdp-1kg/hgdp_tgp/yshi/eur_testsnp_chr{chrom}.vcf.bgz"
hl.export_vcf(mt_eur_testsnp, output_path)
```

### anotating ancestral allele
```shell
/vep -i eur_testsnp_chr${chrom}.vcf.bgz \
--assembly GRCh38 \
--cache \
--offline \
--use_given_ref \
--dir_cache /humgen/atgu1/methods/yshi/vep \
--vcf \
--vcf_info_field AA \
--no_stats --minimal \
--plugin AncestralAllele,/humgen/atgu1/methods/yshi/vep/homo_sapiens_ancestor_GRCh38.fa.gz \
-o eur_testsnp_chr${chrom}_ancestor.vcf
```

### swap ancestral allele
```python
import hail as hl
import pandas as pd

hl.init()

vcf_path = f"gs://hgdp-1kg/hgdp_tgp/yshi/eur_testsnp_chr{chr_num}_ancestor.vcf"
long_txt_path = f"gs://hgdp-1kg/hgdp_tgp/yshi/chr{chr_num}_testsnp_long.txt"
wide_txt_path = f"gs://hgdp-1kg/hgdp_tgp/yshi/EUR_chr{chr_num}_testsnp.txt"

mt_full = hl.import_vcf(vcf_path, reference_genome='GRCh38', force_bgz=True)
mt_full = mt_full.annotate_rows(
    ancestral_allele=mt_full.info.AA[0].split('|')[0]
)

mt = mt_full.select_rows('rsid', 'ancestral_allele').select_entries('GT')

# filter out variant without ancestral allele 
mt = mt.filter_rows(mt.ancestral_allele != '-')

# swap 
mt = mt.annotate_rows(
    swap=mt.alleles[0] == mt.ancestral_allele
)

mt = mt.annotate_rows(
    ref=mt.alleles[0],  
    alt=mt.alleles[1]
)
    
mt = mt.annotate_rows(  
    ref_adj=hl.if_else(mt.swap, mt.alt, mt.ref),  
    alt_adj=hl.if_else(mt.swap, mt.ref, mt.alt)
)

def encode_gt(gt):
    return hl.case() \
    .when(gt.is_hom_ref(), 0) \
    .when(gt.is_het(), 1) \
    .when(gt.is_hom_var(), 2) \
    .default(hl.missing(hl.tint32))

mt = mt.annotate_entries(
    GT=encode_gt(mt.GT),
    gt_adj=hl.if_else(
    mt.swap,  
    2 - mt.GT.n_alt_alleles(),  
    mt.GT.n_alt_alleles()
    )
)

# extract position (without chromosome)
mt = mt.annotate_rows(position=mt.locus.position)

entries_table = mt.entries()
entries_table = entries_table.select('rsid', 'ref_adj', 'alt_adj', 'position', 'gt_adj')
entries_table.export(long_txt_path)

import pandas as pd

# Set chunk size
chunksize = 788000  

print("🚀 Processing all chromosomes...")

# Store chunk results
wide_df_list = []
prev_chunk = None  # Store previous chunk to prevent (`rsid`, `ref_adj`, `alt_adj`, `position`) from being split

# Read data in chunks
reader = pd.read_csv(long_txt_path, sep='\t', dtype={'gt_adj': 'Int64'}, chunksize=chunksize)

for chunk in reader:
    # Merge with previous chunk if necessary (prevent splitting of (`rsid`, `ref_adj`, `alt_adj`, `position`))
    if prev_chunk is not None:
        chunk = pd.concat([prev_chunk, chunk])

    # Identify the last SNP in the chunk
    last_snp = tuple(chunk.iloc[-1][['rsid', 'ref_adj', 'alt_adj', 'position']])

    # Find all rows with the same SNP
    mask = (chunk[['rsid', 'ref_adj', 'alt_adj', 'position']] == last_snp).all(axis=1)
    last_snp_rows = chunk[mask]

    # Remove duplicate rows within the current chunk
    chunk = chunk[~mask]

    # Fill missing rsid values
    chunk['rsid'] = chunk['rsid'].fillna(
        'chr:' + chunk['position'].astype(str) + ':' + chunk['ref_adj'] + ':' + chunk['alt_adj']
    )

    # Convert to wide format
    wide_chunk = chunk.pivot(
        index=['rsid', 'ref_adj', 'alt_adj', 'position'],  
        columns='s',  
        values='gt_adj'  
    )

    # Remove NaN values
    wide_chunk = wide_chunk.dropna()

    # Store processed chunk
    wide_df_list.append(wide_chunk)

    # Save last SNP rows for the next chunk
    prev_chunk = last_snp_rows

# Merge all chunks
if wide_df_list:
    final_wide_df = pd.concat(wide_df_list)

    # Sort by position
    final_wide_df = final_wide_df.sort_values(by='position', ascending=True)

    # Export final wide-format table
    final_wide_df.to_csv(wide_txt_path, sep='\t')
