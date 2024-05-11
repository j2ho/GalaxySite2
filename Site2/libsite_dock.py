#!/usr/bin/env python3

import sys
import os
import Galaxy
import libGalaxy
from libGalaxy import SITE_DB_HH_PDB70_HOME, DB_PDB_HOME, SITE_DB_HOME
from .prep import prepare_lig_mol2
from .site import extract_rsr, Ligand
from .libfr_site import SiteTemplate 

# templ_pdb in upper()
# lig_name = '%s_%s_%d'%(residue.resName(), residue.chainID(), residue.resNo())


def run_dock(job, pdb_fn, templ_pdb, templ_chain, lig_name):
    templ = SiteTemplate('%s_%s'%(templ_pdb,templ_chain),[])
    templ.write('.')
    templ.tm = Galaxy.utils.TM_align(templ.pdb_fn,pdb_fn)
    print (templ.tm.tm)
    templ.max_contact_lig = lig_name
    templ.rewrite('.',lig_name)
    #
    ligand = Ligand(lig_name.split('_')[0])
    ligand.append_templ(templ)

    mol2_fn = prepare_lig_mol2(job, ligand, exclude_H=True)
    mol2_fn = prepare_lig_mol2(job, ligand, exclude_H=False)
    ligand.lig_fn = mol2_fn

    rsr_fn, rsr_pdb = extract_rsr(job, pdb_fn, ligand, re_run=False)
    opt_s = {'weight_type': 'GalaxySite', \
             'rsr_pdb': rsr_pdb.relpath(), \
             'rsr_file': rsr_fn.relpath(),\
             'soften_dock_E': 'yes',\
             'include_initial_conf': 'yes', \
             'csa_bank': '30 30', \
             'csa_seed': '2  15  100  2', \
             'csa_n_opr_s': '2  3  2  3  0  0', \
             }
    #
    Galaxy.galaxy.ligdock(job, pdb_fn, mol2_fn, outfile_prefix=ligand.lig_name, \
                                          lig_name=ligand.lig_name, opt_s=opt_s, re_run=False)
    return 
