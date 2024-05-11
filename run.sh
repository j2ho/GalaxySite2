#!/bin/sh
#SBATCH -J test_b
#SBATCH -p normal.q
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 24
#SBATCH -o log
#SBATCH -e err
#SBATCH --nice=10000
cd /home/j2ho/galaxysite_casp16/L4000_example
python3 /home/j2ho/galaxysite_casp16/run_site.py -p /home/j2ho/projects/casp/sars-cov2-mpro/monomer_af/af_model/ranked_0.pdb -tsm foldseek -t site_result
