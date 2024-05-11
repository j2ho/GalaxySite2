#!/usr/bin/env python3
# modified 2021.08 sohee

import sys
import os
import copy
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator 
from rdkit import DataStructs
import Galaxy
import libGalaxy
from libGalaxy import SITE_DB_HOME
from .libfr_site import SiteTemplate, search_site_template
from Galaxy.utils.libtm import TM_result
#cmp not supported in python3 / functools.cmp_to_key used 
import functools
import subprocess
import numpy as np
import time 

CUT_HA = 0.8
CUT_TBM = 0.6

TM_HA = 0.5
TM_TBM = 0.4
TM_FM = 0.3
TM_CUT = 0.2

AVG_PAIR_DIST_CUTOFF = 10.0
MAX_LIGAND = 3
# TMP
MAX_LIGAND_TMP = 1000
# TMP
MAX_TEMPLATE = 10
BSITE_CUTOFF = 4.0
type_dict = {1:'lipid',2:'non-biological',3:'ions/metal',4:'glycan',0:'general'}

BABEL = '/applic/OpenBabel/current/bin/obabel' # 2.4.1
BABEL_3 = '/applic/OpenBabel/3.1.1/bin/obabel' # 3.1.1
# default FP2
finger_print = ['-xfFP3', '-xfFP4', '-xfMACCS']

class Ligand:
    def __init__(self, ligName):
        self.is_valid = True
        #self.lig_name = ligName
        self.lig_name = ligName.split('_')[0]
        self.lig_type = type_dict[int(ligName.split('_')[1])]
        self.templ_s = []
        self.templ_score_s = [] 
        self.n_lig = 0
        self.score = 0.0
        self.selected = False
        #
        self.is_bsite = False 
        self.bsite_of = ''
        #
    def append_templ(self, templ):
        self.templ_s.append(templ)
        self.n_lig += 1
    def check_positional_variation(self):
        from Galaxy.core.vector import distance
        n_pair = 0
        dist = 0
        for i, x in enumerate(self.templ_s):
            for j, y in enumerate(self.templ_s):
                if i>=j: continue
                n_pair += 1
                dist += distance(x.lig_cntr, y.lig_cntr)
        if n_pair > 0:
            self.avg_pd = dist/float(n_pair)
        else:
            self.avg_pd = 0.0
        #
        if self.avg_pd > AVG_PAIR_DIST_CUTOFF:
            self.is_valid = False
    def get_score(self, query_lig_fn=None, fptype=None, simtype=None):
        self.templ_score_s = [] 
        for templ in self.templ_s:
            if not templ.is_valid: continue
            templ_score = templ.score
            if query_lig_fn is not None: 
                similarity = calculate_similarity(self.lig_name, query_lig_fn) 
                #templ_score_v2 = templ_score*similarity
                templ_score_v2 = similarity # just ligand similarity
            elif query_lig_fn is None: 
                templ_score_v2 = -1.0
            self.templ_score_s.append([templ,float(templ_score),float(templ_score_v2)])
    def report(self):
        wrt = ''
        wrt += 'LIGAND  %3s  %s\n'%(self.lig_name, self.lig_type)
        templ_s = self.templ_s
        for templ_score in templ_s:
            templ = templ_score[0]
            score = templ_score[1]
            score_lig = templ_score[2]
            if not templ.is_valid: continue
            wrt += 'TEMPL   %3s  %4.1f  %4.1f  %s_%s   '%(self.lig_name, score, score_lig, templ.pdb_id, templ.chain_id)
            wrt += '%8.3f %8.3f %8.3f\n'%(\
                 templ.lig_cntr[0], templ.lig_cntr[1], templ.lig_cntr[2])
        return wrt
    def read_info(self, line):
        x = line.strip().split()
        self.score = float(x[2])
    def read_templ(self, pdb, line, templ_home):
        x = line.strip().split()
        ligand_s = [self.lig_name]
        templ = SiteTemplate(x[2], ligand_s)
        templ.is_valid = True
        templ.tm = TM_result()
        templ.tm.tm = float(x[3])
        templ.max_contact_lig = x[4]
        templ.pdb_fn = Galaxy.core.FilePath('templ_s/%s_%s.pdb'%(x[2],x[4].split('_')[0]))
        # write_template
        if not os.path.exists('%s'%templ.pdb_fn):
            templ.write(templ_home)
            tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
            templ.tm = tm
            templ.rewrite(templ_home,self.lig_name) # write lig_name
        #
        self.templ_s.append(templ)

def calculate_similarity_rdkit(temp_lig_name, query_lig, fptype='ecfp', simtype='dice'):
    sdf_file = f'/home/j2ho/DB/site_sdfs/{temp_lig_name}.sdf'
    #sdf_file = f'/home/j2ho/junk/Site/foldtest/temp.sdf'
    lig_mol = Chem.SDMolSupplier(sdf_file)[0]
    if lig_mol is None:
        return 0.0
    q_mol = query_lig # mol object of query ligand 
    if fptype == 'ecfp': 
        fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2) 
        fp_t = fpgen.GetSparseCountFingerprint(lig_mol)
        fp_q = fpgen.GetSparseCountFingerprint(q_mol)
    if fptype == 'rdkit': 
        fp_t = Chem.RDKFingerprint(lig_mol)
        fp_q = Chem.RDKFingerprint(q_mol) 
    if simtype == 'dice': 
        similarity = DataStructs.DiceSimilarity(fp_t, fp_q)
    else: 
        similarity = DataStructs.FingerprintSimilarity(fp_t, fp_q)
    return similarity

def calculate_similarity(temp_lig_name, query_lig, fptype='FP2', simtype='tani'):
    sdf_file = f'/home/j2ho/DB/site_sdfs/{temp_lig_name}.sdf'
    p1=sdf_file
    p2=query_lig
    if fptype == 'FP2':
        cmd = '%s %s %s -ofpt 2> /dev/null'%(BABEL_3,p1,p2)
    else:
        cmd = '%s %s %s -ofpt -%s 2> /dev/null'%(BABEL_3,p1,p2,fptype)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    g = proc.stdout.readlines()
    for ln in g:
        ln = ln.decode()
        if 'Tanimoto' in ln:
            tan=ln.split()[-1]
    return float(tan)

def sort_by_score(data1, data2):
    #return cmp(data1.score, data2.score)
    return (data1.score > data2.score) - (data1.score < data2.score)

def select_ligand(job, pdb, re_run=False, **kwargs):
    out_f_s = {'Site.selected_lig.dat': Galaxy.core.FilePath('%s.selected_lig.dat'%job.title)}
    #
    if 'benchmark' in kwargs:
        benchmark = kwargs['benchmark']
    #
    if (not re_run) and Galaxy.core.file_status(out_f_s):
        job.mkdir('templ_s', tag='SITE_templ_HOME', cd=False)
        job.update(out_f_s)
        sys.stdout.write("INFO: Prediction of bound ligand is already done, and this result will be re-used.\n")
        ligand_s = read_selected_lig(out_f_s['Site.selected_lig.dat'], pdb, job['SITE_templ_HOME'])
        return ligand_s
    #
    if 'search_metal' in kwargs:
        search_metal = kwargs['search_metal']
    else:
        search_metal = False
    #
    if 'search_method' in kwargs:
        search_method = kwargs['search_method'] 
    else:
        search_method = 'foldseek'
    lig_fn = kwargs['lig_fn']
    fptype = kwargs['fptype']
    simtype = kwargs['simtype']
    #
    #job.mkdir('site_fr')
    templ_s = search_site_template(job, pdb.pdb_fn, job.fa_fn, search_metal=search_metal, search_method=search_method,
                                   re_run=re_run, benchmark=benchmark)
    #job.chdir()
    #
    if len(templ_s) == 0:
        sys.stdout.write('INFO: There is no template to detect binding site. Terminate.\n')
        fout = open('%s.selected_lig.dat'%job.title, 'wt')
        fout.write('')
        fout.close()
        return []
    #
    job.mkdir('templ_s', tag='SITE_templ_HOME', cd=False)
    for i in range(len(templ_s)):
        best_templ = templ_s[0]
        best_templ.write(job['SITE_templ_HOME'])
        if best_templ.is_valid:
            templ_s = templ_s[i:]
            break
    #
    if not (best_templ.is_valid):
        sys.stdout.write('INFO: There is no template to detect binding site. Terminate.\n')
        fout = open('%s.selected_lig.dat'%job.title, 'wt')
        fout.write('')
        fout.close()
        return []
    #
    if search_method == 'str':
        best_tm = Galaxy.utils.TM_align(best_templ.pdb_fn, pdb.pdb_fn)
        tm_cutoff = set_TMscore_cutoff(best_tm.tm)
    else:
        #best_tm = Galaxy.utils.TM_align(best_templ.pdb_fn, pdb.pdb_fn)
        #tm_cutoff = set_TMscore_cutoff(best_tm.tm)
        tm_cutoff = TM_CUT
    #
    # Check global structure similarity to query protein
    ligand_s = []
    ligName_list = []
    #
    n_valid = 0
    t1 = time.time()
    totaltm = 0
    for templ in templ_s:
        templ.write(job['SITE_templ_HOME'])
        if not templ.is_valid:
            continue
        #tm1 = time.time()
        tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
        templ.tm = tm
        if len(templ.tm.tr) == 0: 
            templ.is_valid = False
        #tm2 = time.time()
        #totaltm += tm2-tm1
        for lig_name in templ.ligand_s:
            if lig_name in ligName_list:
                done = False
                idx = ligName_list.index(lig_name)
                for tmp_templ in ligand_s[idx].templ_s:
                    if tmp_templ.clu == templ.clu:
                        if tmp_templ.pdb_id == templ.pdb_id: # only exclude if same uniprotid AND pdbid already exists 
                            done = True
                            break
                if done:
                    continue
                tmp = copy.deepcopy(templ)
                ligand_s[idx].append_templ(tmp) # (tmp)
            else:
                ligand = Ligand(lig_name)
                tmp = copy.deepcopy(templ) 
                ligand.append_templ(tmp) # (tmp)
                ligand_s.append(ligand)
                ligName_list.append(lig_name)
        #
        n_valid += 1
    t2 = time.time()
    print (f'copy templ to ligand templ_s: {t2-t1}sec',flush=True)
    #print (f'total tm {totaltm}sec',flush=True)
    if n_valid == 0:
        sys.stdout.write('INFO: There is no template to detect binding site. Terminate.\n')
        fout = open('%s.selected_lig.dat'%job.title, 'wt')
        fout.write('')
        fout.close()
        return []
    #
    query_lig = lig_fn
    if not query_lig == None: 
        if not len(query_lig.split('.')) > 1:
            sdf_file = f'/home/j2ho/DB/site_sdfs/{query_lig}.sdf'
            if not os.path.exists(sdf_file):
                sdf_link = f'https://files.rcsb.org/ligands/download/{query_lig}_ideal.sdf'
                os.system(f'wget {sdf_link} -O /home/j2ho/DB/site_sdfs/{query_lig}.sdf') 
            query_lig = sdf_file
    for ligand in ligand_s:
        sdf_file = f'/home/j2ho/DB/site_sdfs/{ligand.lig_name}.sdf'
        if not os.path.exists(sdf_file): 
            sdf_link = f'https://files.rcsb.org/ligands/download/{ligand.lig_name}_ideal.sdf'
            os.system(f'wget {sdf_link} -O /home/j2ho/DB/site_sdfs/{ligand.lig_name}.sdf') 
        ligand.get_score(query_lig_fn=query_lig, fptype=fptype, simtype=simtype)
    #ligand_s.sort(key=functools.cmp_to_key(sort_by_score), reverse=True)
    #
    new_ligand_s = [] 
    wrt_lig = []
    wrt_nonlig = []
    wrt_metal = []
    i_lig = 0
    #
    #
    from Galaxy.core.vector import distance
    t3 = time.time()
    final_templ_s = {}
    for ligand in ligand_s:
        lig_name = ligand.lig_name
        lig_type = ligand.lig_type
        for templ_score in ligand.templ_score_s:
            is_bsite = True
            templ = templ_score[0]
            score = templ_score[1]
            score_lig = templ_score[2]
            tmp=[]
            #tm1 = time.time()
            templ.get_max_contact_lig(lig_name)
            #if lig_name == 'CA': 
            #    if temmpl.pdb_id == '6YJM':
            #        sys.exit() 
            #tm2 = time.time()
            #totaltm += tm2-tm1
            if not templ.is_valid:
                ligand.is_valid = templ.is_valid
                continue
            templ_name = f'{templ.pdb_id}_{templ.chain_id}'
            #if not templ_name in final_templ_s:
            #    #tm1 = time.time()
            #    tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
                #tm2 = time.time()
                #totaltm += tm2-tm1
            #    templ.tm = tm
            #    final_templ_s[templ_name] = copy.deepcopy(templ)
            #else:
            #    templ = final_templ_s[templ_name]
            for elem in templ.contact_ligs: 
                lig_name = elem[1]
                templ.rewrite(job['SITE_templ_HOME'],lig_name) # write lig_name
                ligcntr = templ.get_lig_cntr(job['SITE_templ_HOME'],lig_name)
                tmp.append([templ, score, score_lig])
            #ligand.templ_s = tmp
                new_ligand_s.append([lig_name.strip(), lig_type, templ, score, score_lig, ligcntr])
    #
    t4 = time.time()
    print (f'check templ validity, rewrite aligned pdb: {t4-t3}sec',flush=True)
    
    new_ligand_s.sort(key=lambda x:x[3], reverse=True)
    for ligand in new_ligand_s:
        lig_name = ligand[0]
        lig_type = ligand[1]
        templ = ligand[2]
        score = ligand[3]
        score_lig = ligand[4]
        ligcntr = ligand[5]
        if lig_type == 'non-biological':
            wrt = ''
            wrt += 'LIGAND  %3s  %s  Temp  Temp*Lig\n'%(lig_name, lig_type)
            wrt += 'TEMPL   %3s  %4.1f  %4.1f  %s_%s   '%(lig_name, score, score_lig, templ.pdb_id, templ.chain_id)
            wrt += '%8.3f %8.3f %8.3f\n'%(\
                 ligcntr[0], ligcntr[1], ligcntr[2])
            wrt_nonlig.append(wrt)
        elif lig_type == 'ions/metal':
            wrt = ''
            wrt += 'LIGAND  %3s  %s  Temp  Temp*Lig\n'%(lig_name, lig_type)
            wrt += 'TEMPL   %3s  %4.1f  %4.1f  %s_%s   '%(lig_name, score, score_lig, templ.pdb_id, templ.chain_id)
            wrt += '%8.3f %8.3f %8.3f\n'%(\
                 ligcntr[0], ligcntr[1], ligcntr[2])
            wrt_metal.append(wrt)
        else: 
            wrt = ''
            wrt += 'LIGAND  %3s  %s  Temp  Temp*Lig\n'%(lig_name, lig_type)
            wrt += 'TEMPL   %3s  %4.1f  %4.1f  %s_%s   '%(lig_name, score, score_lig, templ.pdb_id, templ.chain_id)
            wrt += '%8.3f %8.3f %8.3f\n'%(\
                 ligcntr[0], ligcntr[1], ligcntr[2])
            wrt_lig.append(wrt)

    fout = open('%s.selected_lig.dat'%job.title, 'wt')
    fout.writelines(wrt_lig)
    fout.close()
    #
    fout = open('%s.selected_nonbio.dat'%job.title, 'wt')
    fout.writelines(wrt_nonlig)
    fout.close()
    #
    fout = open('%s.selected_metal.dat'%job.title, 'wt')
    fout.writelines(wrt_metal)
    fout.close()
    #
    os.system('rm *.pdb')
    job.update(out_f_s)
    #
    return ligand_s

def set_TMscore_cutoff(best_tm):
    if best_tm > CUT_HA:
        tm_cutoff = TM_HA
    elif best_tm > CUT_TBM:
        tm_cutoff = TM_TBM
    else:
        tm_cutoff = TM_FM
    return tm_cutoff

def read_selected_lig(selected_fn, pdb, templ_home):
    ligand_s = []
    with open('%s'%selected_fn) as fp:
        for line in fp:
            if line.startswith('LIGAND'):
                x = line.strip().split()
                ligand = Ligand(x[1])
                ligand_s.append(ligand)
                if 'SELECTED' in line:
                    ligand.selected= True
            elif line.startswith('INFO'):
                ligand.read_info(line)
            elif line.startswith('TEMPL'):
                ligand.read_templ(pdb, line, templ_home)
    return ligand_s

def extract_rsr(job, pdb, ligand, re_run=False):
    templ_fn_s = []
    tm_s = []
    for templ in ligand.templ_s:
#        templ_fn_s.append(templ.pdb_fn_lig[ligand.lig_name])
        templ_fn_s.append(templ.pdb_fn)
        tm_s.append(templ.tm.tm)
    #
    rsr_fn = Galaxy.galaxy.site_restraint(job, pdb, ligand.lig_fn, templ_fn_s, tm_s,\
                outfile_prefix='%s_rsr'%ligand.lig_name, lig_name=ligand.lig_name, re_run=re_run)
    #
    return rsr_fn
