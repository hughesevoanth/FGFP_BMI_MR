################################
## FGFP MR
## by: David Hughes
## date: March 14 2022
################################

##################################
## Step 1
## -----------
## working in R
## extract rsids
##
##################################
## Source the paramater file with file paths
source("parameters/snp_lists.sh")
files = ls()
## yengo file paths
w = grep("yengo", ls() )
yengo_files = files[w]
## pruit file paths
w = grep("pruit", ls() )
pruit_files = files[w]

## Extract rsids
rsid = c()
for(file in yengo_files){
	d = read.table(eval(as.symbol(file)), header = TRUE, sep = "\t", as.is = TRUE)
	ids = d[,1]
	cat(paste0("rsids are of length ", length(ids), "\n"))
	rsid = c(rsid, ids)
}

for(file in pruit_files){
	d = read.table(eval(as.symbol(file)), header = TRUE, sep = " ", as.is = TRUE)
	ids = d[,1]
	ids = sapply(ids, function(x){ strsplit(x, split = ":")[[1]][1] })
	cat(paste0("rsids are of length ", length(ids), "\n"))
	rsid = c(rsid, ids)
}

length(rsid) ## 6483
length( unique(rsid) ) ## 5684
rsid = sort( unique(rsid) )
length(rsid) ## 5684

rsid = data.frame(rsid)
write.table(rsid, file = "rsids_2_extract_2_dosage.txt", row.names = FALSE, col.names = FALSE, sep = "", quote = FALSE)

##################################
## STEP 2
## ------------
## executed on BP
## via a SLURM script
##
##################################
#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G

module add apps/qctool/2.0.7

## reset working directory
cd /group/njt-group/DHtemp

## path to needed data
SNPS=/group/njt-group/DHtemp/rsids_2_extract_2_dosage.txt
SAMFILE=/user/work/dh16508/fgfp/Metabolon_GWAS/data/processed/genotypes/bgen/data.sample
BGENS=/user/work/dh16508/fgfp/Metabolon_GWAS/data/processed/genotypes/bgen/
OUTDIR=/group/njt-group/DHtemp/

### Produce a DOSAGE file
qctool -g ${BGENS}data_chr#.bgen -s ${SAMFILE} -og ${OUTDIR}yengo_snps.dosage -ofiletype dosage -os ${OUTDIR}data.sample -incl-rsids ${SNPS} 

### run snpstats
qctool -g ${BGENS}data_chr#.bgen -s ${SAMFILE} -incl-rsids ${SNPS} -snp-stats -osnp yengo.snpstats



