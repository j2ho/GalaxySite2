import os
import sys 
import subprocess

BABEL = '/applic/OpenBabel/current/bin/obabel' # 2.4.1
BABEL_3 = '/applic/OpenBabel/3.1.1/bin/obabel' # 3.1.1
# default FP2
finger_print = ['-xfFP3', '-xfFP4', '-xfMACCS']

def calculate_similarity(temp_lig_name, query_lig):
    sdf_file = f'/home/j2ho/DB/site_sdfs/{temp_lig_name}.sdf'
    p1=sdf_file
    p2=query_lig
    cmd = '%s %s %s -ofpt 2> /dev/null'%(BABEL_3,p1,p2)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    g = proc.stdout.readlines()
    for ln in g:
        ln = ln.decode()
        if 'Tanimoto' in ln:
            tan=ln.split()[-1]
    return float(tan)
