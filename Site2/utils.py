import subprocess

BABEL = '/your/path/to/obabel' 

finger_print = ['-xfFP3', '-xfFP4', '-xfMACCS']

def calculate_similarity(temp_lig_name, query_lig, fptype='FP2'):
    # ligand name = CCD code 
    sdf_file = f'/path/to/ligand/sdf/db/{temp_lig_name}.sdf'
    p1=sdf_file
    p2=query_lig
    if fptype == 'FP2':
        cmd = '%s %s %s -ofpt 2> /dev/null'%(BABEL,p1,p2)
    else:
        cmd = '%s %s %s -ofpt -%s 2> /dev/null'%(BABEL,p1,p2,fptype)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    g = proc.stdout.readlines()
    for ln in g:
        ln = ln.decode()
        if 'Tanimoto' in ln:
            tan=ln.split()[-1]
    return float(tan)