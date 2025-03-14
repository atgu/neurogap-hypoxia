import hailtop.batch as hb
import hail as hl

#####
import argparse
parse = argparse.ArgumentParser()
parse.add_argument('-c', type=str, help='chr', required=True)
args = parse.parse_args()
input_chr = int(args.c)
#####

# get file size to allocate resources for batch jobs accordingly
def get_file_size(file):
    file_info = hl.utils.hadoop_stat(file)
    size_bytes = file_info['size_bytes']
    size_gigs = size_bytes / (1024 * 1024 * 1024)
    return size_gigs

# function to convert bcf to vcf
def select_plink(b, plink_file, snp_file,pop_file):

	j = b.new_job(name="filter plink_file") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  # Docker image with bftools and htslib

	j.declare_resource_group(ofile={
		'bed': '{root}.bed',
		'bim': '{root}.bim',
		'fam': '{root}.fam'})

	j.command(
		f'''plink2 --bfile {plink_file} --extract {snp_file} --keep {pop_file} --make-bed --out {j.ofile}''') 

	return j

# function to convert plink to vcf
def convert_vcf(b, bfile_filtered):

	j = b.new_job(name="convert plink to vcf") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  # Docker image with bftools and plink2

	j.declare_resource_group(ofile={
		'log': '{root}.log',
		'vcf': '{root}.vcf'})

	j.command(
		f'''plink2 --bfile {bfile_filtered} --recode vcf --out {j.ofile}''')

	return j

def vep_ancestral(b, vcf, fa_gz):

	j = b.new_job(name="ancestral allele file") # define job
	j.cpu(4)
	j.storage('50Gi')
	j.image('docker.io/daisylyu/dlyu-vep:v1')  # Docker image with vep

	j.declare_resource_group(ofile={
		'txt': '{root}',
		'html':'{root}_summary.html'})

	j.command(
		f'''tar xzf /root/.vep/homo_sapiens_vep_112_GRCh38.tar.gz -C /root/.vep''')

	j.command(
		f'''vep -i {vcf} --plugin AncestralAllele,{fa_gz} --cache -o {j.ofile}''')

	return j

def allele_freq(b, plink_file):

	j = b.new_job(name="calculate allele frequency") 
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  

	j.declare_resource_group(ofile={
		'afreq': '{root}.afreq',
		'log': '{root}.log'})

	j.command(
		f'''plink2 -bfile {plink_file} --freq --out {j.ofile}''') 

	return j

def DDAF(b, AA_file, selected_freq, unselected1_freq, unselected2_freq):

	j = b.new_job(name="calculate DDAF") 
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-ddaf:v2')  

	j.command(
		f'''python3 calculateDDAF.py --AA {AA_file} --selected {selected_freq} \
		--unselected1 {unselected1_freq} --unselected2 {unselected2_freq} --out {j.ofile}''') 

	return j


if __name__ == '__main__':

	backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://neurogap-bge-imputed-regional/daisy/temp') # set up backend

	chr = input_chr

	b = hb.Batch(backend=backend, name=f'DDAF-chr{chr}') # define batch

	snp_file = b.read_input('gs://neurogap-bge-imputed-regional/daisy/hm3_fistcol.txt')

	pop_file_E=b.read_input('gs://neurogap-bge-imputed-regional/daisy/Ethiopia_IDs_pcafiltered_unrelated2.txt')
	pop_file_U=b.read_input('gs://neurogap-bge-imputed-regional/daisy/Uganda_IDs_pcafiltered_unrelated2.txt')
	pop_file_K=b.read_input('gs://neurogap-bge-imputed-regional/daisy/KEMRI_pca_filtered_unrelated.fam')

	plink_file=b.read_input_group(                                                                                              
    	bed=f'gs://neurogap-bge-imputed-regional/glimpse2/plink/per_chr/NeuroGAP_impted_INFO0.8_wrsIDs_chr{chr}.bed',                                 
        bim=f'gs://neurogap-bge-imputed-regional/glimpse2/plink/per_chr/NeuroGAP_impted_INFO0.8_wrsIDs_chr{chr}.bim',                                 
        fam=f'gs://neurogap-bge-imputed-regional/glimpse2/plink/per_chr/NeuroGAP_impted_INFO0.8_wrsIDs_chr{chr}.fam') 

	fa_gz_file=b.read_input(f'gs://neurogap-bge-imputed-regional/daisy/homo_sapiens_ancestor_GRCh38_zipped/homo_sapiens_ancestor_{chr}.fa.gz')

	plink_size = round(get_file_size(f'gs://neurogap-bge-imputed-regional/glimpse2/plink/per_chr/NeuroGAP_impted_INFO0.8_wrsIDs_chr{chr}.bed'))
	storage = round(plink_size * 1.5)

	plink_file_E= select_plink(b, plink_file,snp_file,pop_file_E)
	plink_file_U= select_plink(b, plink_file,snp_file,pop_file_U)
	plink_file_K= select_plink(b, plink_file,snp_file,pop_file_K)

	vcf_file_E=convert_vcf(b,plink_file_E.ofile)
	vcf_file_U=convert_vcf(b,plink_file_U.ofile)
	vcf_file_U=convert_vcf(b,plink_file_K.ofile)

	ancestral_allele=vep_ancestral(b,vcf_file_E.ofile.vcf, fa_gz_file)

	allele_freq_E=allele_freq(b,plink_file_E.ofile)
	allele_freq_U=allele_freq(b,plink_file_U.ofile)
	allele_freq_K=allele_freq(b,plink_file_K.ofile)

	DDAF=DDAF(b, ancestral_allele.ofile.txt, allele_freq_E.ofile.afreq, allele_freq_U.ofile.afreq, allele_freq_K.ofile.afreq)

	b.write_output(DDAF.ofile, f'gs://neurogap-bge-imputed-regional/daisy/DDAF_result/chr{chr}_DDAF_Ethiopia.csv')

	b.run(wait=False)

	backend.close()
