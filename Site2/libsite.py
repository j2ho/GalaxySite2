#!/usr/bin/env python3

#!/usr/bin/python

import os
import sys
import time
import Galaxy
import libGalaxy
from Galaxy.core import FilePath, define_n_proc, file_status
from . import site
from .prep import prepare_lig_mol2
from .bsite_pred import predict_bsite

def run(job, fa_fn=None, pdb_fn=None, n_proc=None, re_run=False, **kwargs):
    #
    define_n_proc(n_proc)
    #
    if pdb_fn == None and fa_fn is not None:
        sys.stdout.write('ERR: Cannot GalaxySite using sequence %s.\n'%fa_fn)
        if not isinstance(fa_fn, FilePath):
            fa_fn = FilePath(fa_fn)
        #
    elif pdb_fn is not None:
        if not isinstance(pdb_fn, FilePath):
            pdb_fn = FilePath(pdb_fn)
        sys.stdout.write('PROC: Running GalaxySite using structure %s.\n'%pdb_fn)
        pdb = Galaxy.core.PDB(pdb_fn)
        if fa_fn == None:
            # get sequence from pdb_fn
            fa_fn = FilePath('%s.fa'%job.title)
            pdb.write_fasta_file(fa_fn)
    #
    
    job.pdb_fn = pdb.pdb_fn
    job.fa_fn = fa_fn
    #
    out_f_s = run_site(job, pdb, re_run=re_run, **kwargs)
    #
    return out_f_s


def run_site(job, pdb, re_run=False, **kwargs):
    ligand_s = site.select_ligand(job, pdb, re_run=re_run, **kwargs)
    if 'run_dock' in kwargs:
        if kwargs['run_dock'] == False:
            return
    #
    out_f_s = {}
    for ligand in ligand_s:
        if not ligand.selected: continue
        job.mkdir('%s'%ligand.lig_name, tag='%s.HOME'%ligand.lig_name)
        out_f_s['%s.site.pdb'%ligand.lig_name] = FilePath('%s.site.pdb'%(ligand.lig_name))
        out_f_s['%s.bsite.dat'%ligand.lig_name] = FilePath('%s.bsite.dat'%(ligand.lig_name))
        if (not re_run) and file_status(out_f_s):
            sys.stdout.write('INFO: GalaxySite result for ligand %s is existing, and it will be re-used.'%(ligand.lig_name))
            job.chdir_prev()
            continue
        mol2_fn = prepare_lig_mol2(job, ligand, exclude_H=True)
        ligand.lig_fn = mol2_fn
        rsr_fn, rsr_pdb = site.extract_rsr(job, pdb, ligand)
        mol2_fn = prepare_lig_mol2(job, ligand, exclude_H=False)
        opt_s = {'weight_type': 'GalaxySite', \
                 'rsr_pdb': rsr_pdb.relpath(), \
                 'rsr_file': rsr_fn.relpath(),\
                 'soften_dock_E': 'yes',\
                 'include_initial_conf': 'yes', \
                 'csa_bank': '30 30', \
                 'csa_seed': '2  15  100  2', \
                 'csa_n_opr_s': '2  3  2  3  0  0', \
                 }
        optimized, opt_energy = Galaxy.galaxy.ligdock(job, pdb, mol2_fn, outfile_prefix=ligand.lig_name, \
                                          lig_name=ligand.lig_name, opt_s=opt_s, re_run=re_run)
        bsite_s = predict_bsite(pdb[0], optimized[0])
        #
        fout = open('%s.site.pdb'%ligand.lig_name, 'wt')
        fout.writelines(optimized[0].write_in_PDB())
        fout.close()
        #
        fout = open('%s.bsite.dat'%ligand.lig_name, 'wt')
        fout.writelines(bsite_s)
        fout.close()
        #
        job.chdir_prev()
    job.update(out_f_s)
    return out_f_s

