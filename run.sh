#!/bin/sh
#SBATCH -J test
#SBATCH -p normal.q
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -o log
#SBATCH -e err
#SBATCH --nice=10000
python3 run_site.py -p target_receptor.pdb -tsm foldseek -t output
