#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --output=log/03_make_dosage_%A_%a.log
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
    OUTFILE=${DATADIR}geno/dosage/data_chr0${SLURM_ARRAY_TASK_ID}_yengo_pulit.dosage
else
    BGEN=${DATADIR}geno/bgen_w_rsid/data_chr${SLURM_ARRAY_TASK_ID}.bgen
    OUTFILE=${DATADIR}geno/dosage/data_chr${SLURM_ARRAY_TASK_ID}_yengo_pulit.dosage
fi

## Sample file
SAMFILE=${DATADIR}geno/bgen_w_rsid/data.imputed.sample
rsSNPS=${DATADIR}rsids_2_extract.txt

OUTSAMFILE=${DATADIR}geno/dosage/data.sample


### Produce a DOSAGE file
qctool_v2.2.0 -g $BGEN -s $SAMFILE -og $OUTFILE -ofiletype dosage -os $OUTSAMFILE -incl-rsids $rsSNPS


## Concatenate files together
# cat data_chr01_yengo_pulit.dosage > yengo_pulit_snps.dosage
# for file in data_chr0{2..9}_yengo_pulit.dosage; do sed '1d' "$file" >> yengo_pulit_snps.dosage; done
# for file in data_chr{10..22}_yengo_pulit.dosage; do sed '1d' "$file" >> yengo_pulit_snps.dosage; done


