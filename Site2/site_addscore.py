#!/usr/bin/env python3
# modified 2021.08 sohee

import sys
import os
import copy
from rdkit import Chem 
from rdkit import DataStructs
import Galaxy
import libGalaxy
from libGalaxy import SITE_DB_HOME
from .libfr_site import SiteTemplate, search_site_template
from Galaxy.utils.libtm import TM_result
#cmp not supported in python3 / functools.cmp_to_key used 
import functools

CUT_HA = 0.8
CUT_TBM = 0.6

TM_HA = 0.5
TM_TBM = 0.4
TM_FM = 0.3

AVG_PAIR_DIST_CUTOFF = 10.0
MAX_LIGAND = 3
# TMP
MAX_LIGAND_TMP = 1000
# TMP
MAX_TEMPLATE = 10
BSITE_CUTOFF = 4.0
type_dict = {1:'lipid',2:'non-biological',3:'ions/metal',4:'glycan',0:'general ligand'}

class Ligand:
    def __init__(self, ligName):
        self.is_valid = True
        #self.lig_name = ligName
        self.lig_name = ligName.split('_')[0]
        self.lig_type = type_dict[int(ligName.split('_')[1])]
        self.templ_s = []
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
    def get_score(self, lig_q=None, fptype=None, simtype=None):
        self.score = 0.0
        for templ in self.templ_s:
            if not templ.is_valid: continue
            self.score += templ.score
        if lig_q is not None: 
            similarity = calculate_similarity(self.lig_name, lig_q, fptype=fptype, simtype=simtype) 
            self.score = self.score*similarity
    def report(self, rank):
        self.rank = rank
        wrt = ''
        wrt += 'LIGAND  %3s  %s '%(self.lig_name, self.lig_type)
        if (not self.is_bsite) and  (self.bsite_of != ''):
            wrt += 'SAME_BSITE %s '%self.bsite_of
        if rank <= MAX_LIGAND:
            self.selected = True
            wrt += 'SELECTED  %2d %8.3f %8.3f %8.3f\n'%(rank,\
                 self.templ_s[0].lig_cntr[0], self.templ_s[0].lig_cntr[1], self.templ_s[0].lig_cntr[2])
        else:
            #wrt += '\n'
            wrt += '%8.3f %8.3f %8.3f\n'%(\
                 self.templ_s[0].lig_cntr[0], self.templ_s[0].lig_cntr[1], self.templ_s[0].lig_cntr[2])
        #
#        wrt += 'INFO    %3s  %8.2f   %8.3f  %3d\n'%(self.lig_name, self.score,\
#                                                    self.avg_pd, self.n_lig)
        wrt += 'INFO    %3s  %8.2f %3d\n'%(self.lig_name, self.score, self.n_lig)
        #
        for templ in self.templ_s:
            if not templ.is_valid: continue
#            wrt += 'TEMPL   %3s  %s_%s  %6.4f  %s\n'%(self.lig_name, templ.pdb_id, templ.chain_id,\
#                                                      templ.tm.tm, templ.max_contact_lig)
            try:
                wrt += 'TEMPL   %3s  %s_%s  %6.4f  %s\n'%(self.lig_name, templ.pdb_id, templ.chain_id,\
                                                          templ.tm.tm, templ.max_contact_lig)
            except:
                wrt += 'TEMPL   %3s  %s_%s\n'%(self.lig_name, templ.pdb_id, templ.chain_id)
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

def calculate_similarity(temp_lig_name, query_lig, fptype='ecfp', simtype='tan'):
    similarity = 0.8
    sdf_file = f'sdfs/{temp_lig_name}.sdf'
    lig_mol = Chem.SDMolSupplier(sdf_file)
    q_mol = Chem.MolFromMol2File(query_lig)
    fp_t = Chem.RDKFingerprint(lig_mol[0])
    fp_q = Chem.RDKFingerprint(q_mol) 
    similarity = DataStructs.FingerprintSimilarity(fp_t, fp_q) 
    return similarity

def sort_by_score(data1, data2):
    #return cmp(data1.score, data2.score)
    return (data1.score > data2.score) - (data1.score < data2.score)

def select_ligand(job, pdb, re_run=False, **kwargs):
    out_f_s = {'Site.selected_lig.dat': Galaxy.core.FilePath('%s.selected_lig.dat'%job.title)}
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
    job.mkdir('site_fr')
    templ_s = search_site_template(job, pdb.pdb_fn, job.fa_fn, search_metal=search_metal, search_method=search_method, re_run=re_run)
    job.chdir()
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
        #tm_cutoff = set_TMscore_cutoff(best_templ.score)
        best_tm = Galaxy.utils.TM_align(best_templ.pdb_fn, pdb.pdb_fn)
        tm_cutoff = set_TMscore_cutoff(best_tm.tm)
    else:
        best_tm = Galaxy.utils.TM_align(best_templ.pdb_fn, pdb.pdb_fn)
        tm_cutoff = set_TMscore_cutoff(best_tm.tm)
    #
    # Check global structure similarity to query protein
    ligand_s = []
    ligName_list = []
    #
    n_valid = 0
# DELETED
#        templ.check_position(pdb) # close to protein 
#        if not templ.is_valid:
#            continue
    for templ in templ_s:
        # template tm cutoff
#        templ.write(job['SITE_templ_HOME'])
#        try:
#            if templ.tm.tm < tm_cutoff:
#                cotninue
#        except:
#            tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
#            if tm.tm < tm_cutoff:
#                continue
#            templ.tm = tm
        templ.write(job['SITE_templ_HOME'])
        if not templ.is_valid:
            continue
        tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
        templ.tm = tm
        if tm.tm < tm_cutoff:
            continue
        #
        for lig_name in templ.ligand_s:
            if lig_name in ligName_list:
                done = False
                idx = ligName_list.index(lig_name)
                for tmp_templ in ligand_s[idx].templ_s:
                    if tmp_templ.clu == templ.clu:
                        done = True
                        break
                if done:
                    continue
                tmp = copy.deepcopy(templ)
                ligand_s[idx].append_templ(tmp)
            else:
                ligand = Ligand(lig_name)
                tmp = copy.deepcopy(templ)
                ligand.append_templ(tmp)
                ligand_s.append(ligand)
                ligName_list.append(lig_name)
        #
        n_valid += 1
    if n_valid == 0:
        sys.stdout.write('INFO: There is no template to detect binding site. Terminate.\n')
        fout = open('%s.selected_lig.dat'%job.title, 'wt')
        fout.write('')
        fout.close()
        return []
    #
    # Check the positional variation of the ligand
    job.mkdir('sdfs', tag='SITE_sdf_HOME', cd=False)
    for ligand in ligand_s:
#        ligand.check_positional_variation()
#        if not ligand.is_valid: continue
        sdf_file = f'sdfs/{ligand.lig_name}.sdf'
        if not os.path.exists(sdf_file): 
            sdf_link = f'https://files.rcsb.org/ligands/download/{ligand.lig_name}_ideal.sdf'
            os.system(f'wget {sdf_link} -O sdfs/{ligand.lig_name}.sdf') 
        ligand.get_score(lig_q=lig_fn, fptype=fptype, simtype=simtype)
    ligand_s.sort(key=functools.cmp_to_key(sort_by_score), reverse=True)
    #
    bsite = {}
    wrt = []
    i_lig = 0
    #
    nodock=[ln.strip() for ln in open('%s/nodock'%SITE_DB_HOME).readlines()]
    #
    from Galaxy.core.vector import distance
    for ligand in ligand_s:
        lig_name = ligand.lig_name
        num_heavy = []
#        if i_lig <= MAX_LIGAND:
        # TMP RETURN ALL BSITE
        if len(bsite) < MAX_LIGAND_TMP:
        # TMP
        #if len(bsite) < MAX_LIGAND:
            tmp=[]
            is_bsite = True
            for templ in ligand.templ_s:
#                templ.write(job['SITE_templ_HOME'])
#                if not templ.is_valid:
#                    ligand.is_valid = templ.is_valid
#                    continue
#                tm = Galaxy.utils.TM_align(templ.pdb_fn, pdb.pdb_fn)
#                templ.tm = tm
                templ.get_max_contact_lig(lig_name)
                if not templ.is_valid:
                    ligand.is_valid = templ.is_valid
                    continue
                templ.rewrite(job['SITE_templ_HOME'],lig_name) # write lig_name
                templ.check_position(pdb)
                if not templ.is_valid:
                    ligand.is_valid = templ.is_valid
                    continue
                #
                if templ.num_heavy not in num_heavy:
                    num_heavy.append(templ.num_heavy)
                #templ.rewrite() # write lig_name
                templ.get_lig_cntr()
                tmp.append(templ)
                if len(tmp) == MAX_TEMPLATE:
                    break
            ligand.templ_s = tmp
            #
            if not ligand.is_valid: continue
            if max(num_heavy) < 5:
                is_bsite = False
                ligand.is_valid = False
                continue
            if len(num_heavy) > 1:
                tmp = []
                for templ in ligand.templ_s:
                    if templ.num_heavy != max(num_heavy):
                        continue
                    tmp.append(templ)
                ligand.templ_s = tmp
            # check_positional_variation
            ligand.check_positional_variation()
            if not ligand.is_valid: continue
            #
            templ = ligand.templ_s[0]
            for tmp_lig in bsite:
                R1 = bsite[tmp_lig]
                R2 = templ.lig_cntr
                if distance(R1,R2) < BSITE_CUTOFF:
                    is_bsite= False
                    break
            if ligand.lig_name.strip() in nodock:
                is_bsite = False
            elif is_bsite:
                bsite[lig_name] = templ.lig_cntr
                ligand.is_bsite = True
            else:
                ligand.bsite_of = tmp_lig
            #
        if ligand.is_valid:
            # TMP RETURN ALL BSITE
            if len(bsite) > MAX_LIGAND:
                wrt.append(ligand.report(100))
                ligand.is_bsite = False
                continue
            # TMP
            if ligand.is_bsite:
                i_lig=i_lig+1
                wrt.append(ligand.report(i_lig))
            else:
                wrt.append(ligand.report(100))
    #
    fout = open('%s.selected_lig.dat'%job.title, 'wt')
    fout.writelines(wrt)
    fout.close()
    #
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
