#!/usr/bin/env python3

import os
import sys
import subprocess
import Galaxy
from Galaxy.core import FilePath
from .config import BABEL

metal_s = ["NA", "MG", "AL", "K", "CA", "MN", "MN3", "FE2", "FE", "CO", "3CO", "NI", "3NI",\
           "CU1", "CU", "ZN", "AG", "CD", "PT", "AU", "AU3", "HG"]

def prepare_lig_mol2(job, ligand, exclude_H=False):
    templ = ligand.templ_s[0]
    lig, stat = copy_ligand(templ.pdb_fn, templ.max_contact_lig)
    if not stat:
        sys.stdout.write('ERROR: Failed to extract ligand from template PDB.')
        raise RuntimeError
    fout = open('%s.pdb'%ligand.lig_name, 'wt')
    fout.write(lig)
    fout.close()
    if (exclude_H):
        Galaxy.tools.babel.run(job, '%s.pdb'%ligand.lig_name, lig_name=ligand.lig_name,\
                                    out_format='mol2', delete_H=exclude_H, re_run=True)

        cmd = '%s -imol2 %s -osmi'%(BABEL,FilePath('%s.mol2'%ligand.lig_name).relpath())
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        g = '%s'%proc.stdout.read().split()[0]
        if '.' in g:
            Galaxy.tools.chimera.DockPrep(job, '%s.pdb'%ligand.lig_name, lig_name=ligand.lig_name,\
                                     out_format='mol2',re_run=True, report_status=True)
    else:
        Galaxy.tools.chimera.DockPrep(job, '%s.mol2'%ligand.lig_name, lig_name=ligand.lig_name,\
                                     out_format='mol2',re_run=True, report_status=True)

    cmd = '%s -imol2 %s -osmi'%(BABEL,FilePath('%s.mol2'%ligand.lig_name).relpath())
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    g = '%s'%proc.stdout.read().split()[0]
    if '.' in g:
        sys.stdout.write('ERROR: Failed to extract ligand from template PDB.')
        raise RuntimeError

    line_s = [line.strip() for line in open('%s.mol2'%ligand.lig_name)]
    if len(line_s) < 3:
        os.remove('%s.mol2'%ligand.lig_name)
        Galaxy.tools.babel.run(job, '%s.pdb'%ligand.lig_name, lig_name=ligand.lig_name,\
                               out_format='mol2', delete_H=exclude_H, re_run=True)

    return Galaxy.core.FilePath('%s.mol2'%ligand.lig_name)


def extract_ligand(pdb_fn, lig_info):
    ligName = lig_info.split('_')[0]
    ligChain = lig_info.split('_')[1]
    ligNo = int(lig_info.split('_')[2])
    if len(ligChain) == 0: ligChain = ' '

    pdb = Galaxy.core.PDB(pdb_fn)[0]
    find = False
    include_metal = False
    for residue in pdb.get_residues(res_range=[ligNo]):
        if residue.isAtom(): continue
        if (residue.chainID() == ligChain) and (residue.resName() == ligName):

            a=residue.write_full().replace("'",' ')
            a=a.split('\n')
            for i, ln in enumerate(a):
                if len(ln[12:16].split()) > 1:
                    a[i] = '%s%-4s%s'%(ln[:12],''.join(ln[12:16].split()),ln[16:])
                if len(a[i]) > 78:
                    todel = len(a[i])-78
                    a[i]='%s%s'%(a[i][:66],a[i][66+todel:])
            a='\n'.join(a)
            wrt = '%s'%a

            find = True
            for ia in range(len(residue)):
                if residue._atmName[ia] in metal_s:
                    include_metal = True
                    break
            break
    return wrt, find, include_metal


def copy_ligand(pdb_fn, lig_info):
    ligName = lig_info.split('_')[0]
    ligChain = lig_info.split('_')[1]
    ligNo = int(lig_info.split('_')[2])
    if len(ligChain) == 0: ligChain = ' '

    wrt=''
    find=False
    for ln in open('%s'%pdb_fn).readlines():
        if not ln.startswith('HETATM'):
            continue
        if ln[21] != ligChain:
            continue
        if ln[17:20].strip() != ligName.strip():
            continue
        if int(ln[22:26]) != ligNo:
            continue
        wrt+=ln
        find = True
    return wrt, find