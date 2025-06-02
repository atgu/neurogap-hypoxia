### files preparation
prepare gene list
```shell
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
# all gene location
zcat gencode.v44.annotation.gtf.gz | \
awk '$3 == "gene" {
    match($0, /gene_id "([^"]+)"/, arr);
    print $1, $4-1, $5, arr[1]
}' OFS='\t' > all_genes.bed
# coding gene location
zcat gencode.v44.annotation.gtf.gz | \
awk '$3 == "gene" && $0 ~ /gene_type "protein_coding"/ { 
    match($0, /gene_id "([^"]+)"/, arr); 
    print $1, $4-1, $5, arr[1] 
}' OFS='\t' > protein_coding_genes.bed
```

