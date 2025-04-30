#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --output=log/03_make_plink_%A_%a.log
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --account=sscm015962
#SBATCH --array=1-22

## Recommended feedback
echo Running on host "$(hostname)"
echo Time is "$(date)"
echo Directory is "$(pwd)"
echo Slurm job ID is "${SLURM_JOBID}"
echo This jobs runs on the following machines:
echo "${SLURM_JOB_NODELIST}"

## Report which job is being run
echo "Running task ${SLURM_ARRAY_TASK_ID}"

## add qctools to my environment
module add qctool/2.2.0

## Move to the directory the job was submitted from
cd "${SLURM_SUBMIT_DIR}"

################
## source my paramater file
## - it defines the DATADIR and OUTDIR
################
source parameter_file.txt

################
## Define job variables
################
## bgen and outfile names
if [ "$SLURM_ARRAY_TASK_ID" -ge 1 ] && [ "$SLURM_ARRAY_TASK_ID" -le 9 ]; then
    BGEN=${DATADIR}geno/bgen_w_rsid/data_chr0${SLURM_ARRAY_TASK_ID}.bgen
    OUTFILE=${DATADIR}geno/plink/data_chr0${SLURM_ARRAY_TASK_ID}_yengo_pulit
else
    BGEN=${DATADIR}geno/bgen_w_rsid/data_chr${SLURM_ARRAY_TASK_ID}.bgen
    OUTFILE=${DATADIR}geno/plink/data_chr${SLURM_ARRAY_TASK_ID}_yengo_pulit
fi

## Sample file
SAMFILE=${DATADIR}geno/bgen_w_rsid/data.imputed.sample
rsSNPS=${DATADIR}rsids_2_extract.txt


### Produce a PLINK file
qctool_v2.2.0 -g $BGEN -s $SAMFILE -og $OUTFILE -ofiletype binary_ped -incl-rsids $rsSNPS


# ## name of files to merge
# # ls *.bed | sed "s/\.bed\$//" > files2merge.txt
# for file in *.bed; do echo ${file%.bed}.bed ${file%.bed}.bim ${file%.bed}.fam; done >> files2merge.txt

# ## merge all the chromosomes together
# plink --bfile data_chr01_yengo_pulit --merge-list files2merge.txt --make-bed --out yengo_pulit_allchrs --noweb


# ## merge all the chromosomes together
# cat data_chr01_yengo_pulit.fam > merged/yengo_pulit.fam
# for i in {1..9}; do cat data_chr0${i}_yengo_pulit.bim; done > merged/yengo_pulit.bim
# (echo -en "\x6C\x1B\x01"; for chr in {1..9}; do tail -c +4 data_chr0${i}_yengo_pulit.bed; done) > merged/yengo_pulit.bed

# for i in {10..22}; do cat data_chr${i}_yengo_pulit.bim; done >> merged/yengo_pulit.bim
# (echo -en "\x6C\x1B\x01"; for chr in {10..22}; do tail -c +4 data_chr${i}_yengo_pulit.bed; done) >> merged/yengo_pulit.bed

