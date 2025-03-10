hl.init()

# extract EUR IDs
ht_path = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"
ht = hl.read_table(ht_path)
eur_samples = ht.filter(ht.hgdp_tgp_meta.genetic_region == 'EUR')
eur_sample_list = eur_samples.s.collect()
eur_samples_set = hl.literal(set(eur_sample_list))

# extract EUR samples
vcf_path = f"gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr{chrom}.vcf.bgz"
mt = hl.import_vcf(vcf_path, reference_genome='GRCh38')    
mt_eur = mt.filter_cols(eur_samples_set.contains(mt.s))

# annotation
mt_eur = mt_eur.drop(mt_eur.info)
mt_eur = mt_eur.annotate_rows(
    AC=hl.agg.sum(mt_eur.GT.n_alt_alleles()),  
    AN=hl.agg.count_where(hl.is_defined(mt_eur.GT)) * 2,  
)
mt_eur = mt_eur.annotate_rows(AF=mt_eur.AC / mt_eur.AN)

# extract singleton & SNV
mt_eur_singletons = mt_eur.filter_rows(mt_eur.AC == 1)
mt_eur_singletons = mt_eur_singletons.filter_rows(
    (hl.len(mt_eur_singletons.alleles[0]) == 1) & (hl.len(mt_eur_singletons.alleles[1]) == 1)
)

# extract samples with singletons (every sample should have singletons)
mt_eur_singletons = mt_eur_singletons.filter_entries(
    hl.is_defined(mt_eur_singletons.GT) & mt_eur_singletons.GT.is_non_ref()
)

# output
singletons_table = mt_eur_singletons.entries()
sample_singletons = singletons_table.group_by(singletons_table.s).aggregate(
    positions=hl.agg.collect(singletons_table.locus.position)
)
output_path = f"sample_singletons_chr{chrom}.txt"
with open(output_path, "w") as f:
    for row in sample_singletons.collect():
        sample = row.s
        positions_sorted = sorted(row.positions) 
        f.write(f"{sample} {' '.join(map(str, positions_sorted))}\n")
