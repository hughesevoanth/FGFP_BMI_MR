#!/bin/bash

#SBATCH --job-name=test_job
#SBATCH --output=00_run_snpstats_%A_%a.log
#SBATCH --partition=_____
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --account=_____
#SBATCH --mpi=pmi2
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
    BGEN=${DATADIR}data_chr0${SLURM_ARRAY_TASK_ID}.bgen
    OUTFILE=${OUTDIR}data_chr0${SLURM_ARRAY_TASK_ID}.snpstats
else
    BGEN=${DATADIR}data_chr${SLURM_ARRAY_TASK_ID}.bgen
    OUTFILE=${OUTDIR}data_chr${SLURM_ARRAY_TASK_ID}.snpstats
fi

## Sample file
SAMFILE=${DATADIR}data.imputed.sample

### run snpstats
qctool_v2.2.0 -g $BGEN -s $SAMFILE -snp-stats -osnp $OUTFILE

