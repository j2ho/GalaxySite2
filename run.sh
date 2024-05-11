#!/bin/sh
#SBATCH -J test
#SBATCH -p normal.q
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -o log
#SBATCH -e err
#SBATCH --nice=10000
cd test
python3 /home/j2ho/projects/GalaxySite2/run_site.py -p /home/j2ho/projects/casp/sars-cov2-mpro/monomer_af/af_model/ranked_0.pdb -tsm foldseek -t result
