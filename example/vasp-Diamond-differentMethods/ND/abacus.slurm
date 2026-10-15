#!/bin/bash

#SBATCH --job-name=e3rd-test-nd-vasp
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem=32GB
#SBATCH -p hzw
#SBATCH -o e3rd-%j.out

module load vasp/5.4.4

elastic3rd nd
#RunCASTEP.sh -np $SLURM_NTASKS Fe

#mpirun -np $SLURM_NTASKS abacus

#python test_geomopt_isifs.py

