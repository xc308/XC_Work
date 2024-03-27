#!/bin/bash
#SBATCH --account=vjgo8416-xchen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --qos=turing


apptainer exec ../r-base_4.2.3.sif Rscript Bask_Test_1.R
