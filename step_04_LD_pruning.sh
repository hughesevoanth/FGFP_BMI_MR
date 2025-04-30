#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --output=log/04_LD_pruning_%A_%a.log
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
module add plink/1.07
# module add qctool/2.2.0

## Move to the directory the job was submitted from
cd "${SLURM_SUBMIT_DIR}"

################
## source my paramater file
## - it defines the DATADIR and OUTDIR
################
source parameter_file.txt

## move to the out folder
#cd ${DATADIR}gwas_snp_lists/20241111/ldprune

################
## Define job variables
################

## File with the RSIDS
RSIDS=${DATADIR}gwas_snp_lists/20241111/yengo_bmi_cojo_rsids.txt


## bgen and outfile names
if [ "$SLURM_ARRAY_TASK_ID" -ge 1 ] && [ "$SLURM_ARRAY_TASK_ID" -le 9 ]; then
    BED=${DATADIR}geno/plink/data_chr0${SLURM_ARRAY_TASK_ID}_yengo_pulit
    OUTFILE=${DATADIR}gwas_snp_lists/20241111/ldprune/data_chr0${SLURM_ARRAY_TASK_ID}_yengo_bmi_cojo
else
    BED=${DATADIR}geno/plink/data_chr${SLURM_ARRAY_TASK_ID}_yengo_pulit
    OUTFILE=${DATADIR}gwas_snp_lists/20241111/ldprune/data_chr${SLURM_ARRAY_TASK_ID}_yengo_bmi_cojo
fi


### Produce a DOSAGE file
plink --noweb --bfile $BED --indep-pairwise 50 5 0.001 --allow-no-sex --extract $RSIDS --out $OUTFILE


