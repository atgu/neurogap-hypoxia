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
def filter_vcf(b, vcfgz_file, snp_file,pop_file):

	j = b.new_job(name="filter vcf_file") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  

	j.command(f'cp {snp_file} ./hm3_bcftool.txt')
	j.command('ls')
	j.command(
		f'''bcftools view --include 'ID=@hm3_bcftool.txt' -S {pop_file} -o {j.ofile} -O z --force-samples {vcfgz_file}''') 

	return j

# function to convert plink to vcf
def convert_plink(b, vcf_filtered):

	j = b.new_job(name="convert vcf to plink") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  # Docker image with bftools and plink2

	j.declare_resource_group(ofile={
		'bed': '{root}.bed',
		'bim': '{root}.bim',
		'fam': '{root}.fam'})

	j.command(
		f'''plink2 --vcf {vcf_filtered} --make-bed --double-id --out {j.ofile}''')

	return j

def generate_map(b, bim_file,chr,map): 

	j = b.new_job(name="generating mapfile") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-selscan')  # Docker image with selscan, makeMap.py, python2, python3

	j.command(
		f'''python makeMap.py --chr {chr} --genmap {map} --bim {bim_file} --map_bim map --out {j.ofile}''')

	return j

def edit_fafile(b,fagz_file,chr): 

	j = b.new_job(name="edit ancestral.fa file") # define job
	j.cpu(4)
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  # Docker image with selscan, makeMap.py, python2, python3
	
	j.command(
		f'''gunzip -c {fagz_file} > temp_file.fa''')
	j.command(
		f'''sed "1s/.*/>chr{chr}/" temp_file.fa > {j.ofile}''')

	return j

def swap_ancestral(b,vcf_file,fa_file): 

	j = b.new_job(name="make ref/alt ance/deri") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-plinkbcf')  # Docker image with selscan, makeMap.py, python2, python3

	j.command(
		f'''bcftools +fixref {vcf_file} -Ov -o {j.ofile} -- -f {fa_file} -m flip''')

	return j

# run selscan

def selscan_ihs(b, vcf_file_E, mapfile):

	j = b.new_job(name="running selscan_ihs") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-selscan') 

	j.declare_resource_group(ofile={
		'log': '{root}.ihs.log',
		'out': '{root}.ihs.out'})

	j.command(
		f'''selscan --ihs --vcf {vcf_file_E} --map {mapfile} --threads 4 --out {j.ofile}''')


	return j

def selscan_xpehh(b, vcf_file_E, vcf_file_U,mapfile):

	j = b.new_job(name="running selscan_xpehh") # define job
	j.cpu(4)
	j.storage(f'{storage}Gi')
	j.image('docker.io/daisylyu/dlyu-selscan') 

	j.declare_resource_group(ofile={
		'log': '{root}.xpehh.log',
		'out': '{root}.xpehh.out'})

	j.command(
		f'''selscan --xpehh --vcf {vcf_file_E}  --vcf-ref {vcf_file_U} --map {mapfile} --threads 4 --out {j.ofile}''')

	return j


if __name__ == '__main__':

	backend = hb.ServiceBackend(billing_project='diverse-pop-seq-ref', remote_tmpdir='gs://neurogap-bge-imputed-regional/daisy/temp') # set up backend

	chr = input_chr
	b = hb.Batch(backend=backend, name=f'selection metrics-chr{chr}') # define batch
	#read inputs
	snp_file = b.read_input('gs://neurogap-bge-imputed-regional/daisy/NeuroGap_data/hm3_bcftool.txt')

	map_file = b.read_input(f'gs://neurogap-bge-imputed-regional/daisy/recomb-hg38/plink.chr{chr}.GRCh38.map')

	pop_file_E=b.read_input('gs://neurogap-bge-imputed-regional/daisy/NeuroGap_data/Ethiopia_IDs_pcafiltered_unrelated.txt')
	pop_file_U=b.read_input('gs://neurogap-bge-imputed-regional/daisy/NeuroGap_data/Uganda_IDs_pcafiltered_unrelated.txt')
	pop_file_K=b.read_input('gs://neurogap-bge-imputed-regional/daisy/NeuroGap_data/Kemri_IDs_pcafiltered_unrelated.txt')

	vcf_file=b.read_input(f'gs://neurogap-bge-imputed-regional/glimpse2/INFO0.8_filtered/neurogap_final_merged_chr{chr}_INFO0.8_rsids.vcf.gz') 

	ancestral_fa=b.read_input(f'gs://neurogap-bge-imputed-regional/daisy/homo_sapiens_ancestor_GRCh38_zipped/homo_sapiens_ancestor_{chr}.fa.gz')
	#set size
	vcf_size = round(get_file_size(f'gs://neurogap-bge-imputed-regional/glimpse2/INFO0.8_filtered/neurogap_final_merged_chr{chr}_INFO0.8_rsids.vcf.gz'))
	storage = round(vcf_size * 1.5)
	
	#filter vcf
	vcf_file_E= filter_vcf(b, vcf_file,snp_file,pop_file_E)
	vcf_file_U= filter_vcf(b, vcf_file,snp_file,pop_file_U)
	vcf_file_K= filter_vcf(b, vcf_file,snp_file,pop_file_K)
	#generate map file
	plink_file_E=convert_plink(b,vcf_file_E.ofile)
	selscan_map=generate_map(b,plink_file_E.ofile.bim,chr,map_file)
	#make ref/alt ancestral/derived
	edited_fafile=edit_fafile(b,ancestral_fa,chr)
	vcf_ancestral_E=swap_ancestral(b,vcf_file_E.ofile,edited_fafile.ofile)
	vcf_ancestral_U=swap_ancestral(b,vcf_file_U.ofile,edited_fafile.ofile)
	vcf_ancestral_K=swap_ancestral(b,vcf_file_K.ofile,edited_fafile.ofile)	

	ihs_E= selscan_ihs(b, vcf_ancestral_E.ofile,selscan_map.ofile)
	xpehh_EU= selscan_xpehh(b, vcf_ancestral_E.ofile, vcf_ancestral_U.ofile,selscan_map.ofile)
	xpehh_EK= selscan_xpehh(b, vcf_ancestral_E.ofile, vcf_ancestral_K.ofile,selscan_map.ofile)

	b.write_output(ihs_E.ofile, f'gs://neurogap-bge-imputed-regional/daisy/NeuroGap_summary_stats/ihs_phased/chr{chr}_Ethiopian_phased.ihs')
	b.write_output(xpehh_EU.ofile, f'gs://neurogap-bge-imputed-regional/daisy/NeuroGap_summary_stats/xpehh_phased/chr{chr}_Ethi_Ugan_phased.xpehh')
	b.write_output(xpehh_EK.ofile, f'gs://neurogap-bge-imputed-regional/daisy/NeuroGap_summary_stats/xpehh_phased/chr{chr}_Ethi_Kemri_phased.xpehh')

	b.run(wait=False)

	backend.close()
