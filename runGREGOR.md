### GREGOR source code
link https://csg.sph.umich.edu/GREGOR/index.php/site/downloadFile?fileTag=source_codes

### Download reference files
download link https://csg.sph.umich.edu/GREGOR/index.php/site/download

remember to check the md5 after downloading

then merge the part files into one gz file
``` shell
 cat \
   GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz.part.00 \
   GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz.part.01 \
   GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz.part.02 \
   GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz.part.03 \
   > GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz
```
then extract this file 
```shell
 tar zxvf GREGOR.EUR.ref.r2.greater.than.0.7.tar.gz
```
You will get one directory which has the name "EUR", the path to this EUR folder will be your REF_DIR (the path do not include EUR)

### Prepare index file
clump the sum stats to get index SNPs
```shell
plink --bfile ref_panel --clump summary_statistitcs --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.1 --maf 0.1 --out clumped_sum_stats
```
the index file looks like 
```
rs10184514
rs7145828
```

### Prepare bed file 
the bed file looks like (better to be sorted based on chromosome and position)
```
1       2406440 2436439
1       6570061 6600060
1       7980061 8160060
1       11190058        11220057
1       35745602        35820601
```
the bed file index looks like 
```
/humgen/atgu1/methods/yshi/GregorRepro/Ethio_Ugan.bed
/humgen/atgu1/methods/yshi/GregorRepro/Ethio_Kemri.bed
/humgen/atgu1/methods/yshi/GregorRepro/TibetREGION_sorted.bed 
```

### Prepare the conf file
```shell
INDEX_SNP_FILE = /humgen/atgu1/methods/yshi/GregorRepro/index/GCST006046_indexSNP
BED_FILE_INDEX = /humgen/atgu1/methods/yshi/GregorRepro/random.bed.file.index
REF_DIR = /humgen/atgu1/methods/yshi/Ethiopian/
R2THRESHOLD = 0.7
LDWINDOWSIZE = 1000000
OUT_DIR = /humgen/atgu1/methods/yshi/GregorRepro/GCST006046/
MIN_NEIGHBOR_NUM = 500
BEDFILE_IS_SORTED = FALSE
POPULATION = EUR
TOPNBEDFILES = 1
JOBNUMBER = 1
BATCHTYPE = local
```

### Run on server
```shell
#! /bin/bash

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse

use UGER
use .perl-5.8.9

cd /humgen/atgu1/methods/yshi/Ethiopian/GREGOR/script
perl GREGOR.pl --conf /humgen/atgu1/methods/yshi/GregorRepro/conf_file/conf.file
```
