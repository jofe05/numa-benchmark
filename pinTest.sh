#!/bin/bash

#SBATCH --job-name="pin TEST"
#SBATCH -w penguin
#SBATCH --time=12:00:00
#SBATCH --partition=p_hpca4se 
#SBATCH -o /scratch/077-hpca4se-bioinf/jlenis/software/custom/%j_benchtrace.log 
#SBATCH --mem=120GB
#SBATCH -c 64 
#cd /scratch/077-hpca4se-bioinf/jlenis/software/pin/pin-2.14-71313-gcc.4.4.7-linux/source/tools/SimpleExamples/
cd /scratch/077-hpca4se-bioinf/jlenis/software/pin/pin-2.14-71313-gcc.4.4.7-linux/third-party/atomic-memory-trace/trace/
pin -t obj-intel64/trace.so -o /scratch/077-hpca4se-bioinf/jlenis/software/custom/$SLURM_JOB_ID.txt -- /scratch/077-hpca4se-bioinf/jlenis/software/custom/test_numa


