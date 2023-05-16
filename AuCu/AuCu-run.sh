#!/bin/bash -l
#SBATCH -p nodes
#SBATCH -n 4
#SBATCH --mem-per-cpu=4200
#SBATCH --exclusive
module purge
module add impi sci/dft sci/qe_7.2
export MKL_NUM_THREADS=1
export ASE_ESPRESSO_COMMAND="/home/sci/opt/qe-7.2_impi/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"
python AuCu-mc.py
