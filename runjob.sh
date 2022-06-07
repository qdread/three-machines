#!/bin/bash
#SBATCH --job-name=job
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --partition=short
#SBATCH --time=1-00:00:00

cd /home/quentin.read/GitHub/ars-misc
module load r/4.1.2
Rscript2 ${scriptname}

