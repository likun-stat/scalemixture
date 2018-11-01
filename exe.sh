
#!/bin/bash

#PBS -A drh20_a_g_sc_default
#PBS -V
#PBS -l nodes=1:ppn=20
#PBS -l walltime=200:00:00
#PBS -N scalemixture

module purge
module load r/3.4

cd $PBS_O_WORKDIR/

Rscript imputation_animation_02.R
