#!/bin/bash

#SBATCH --job-name=make_dosage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G

module add apps/qctool/2.0.7

## reset working directory
cd /group/njt-group/DHtemp

## path to needed data
SNPS=rsids_2_extract_2_dosage.txt
SAMFILE=data/processed/genotypes/bgen/data.sample
BGENS=data/processed/genotypes/bgen/
OUTDIR=results/

### Produce a DOSAGE file
qctool -g ${BGENS}data_chr#.bgen -s ${SAMFILE} -og ${OUTDIR}anthro_snps.dosage -ofiletype dosage -os ${OUTDIR}data.sample -incl-rsids ${SNPS} 

### run snpstats
qctool -g ${BGENS}data_chr#.bgen -s ${SAMFILE} -incl-rsids ${SNPS} -snp-stats -osnp yengo.snpstats

