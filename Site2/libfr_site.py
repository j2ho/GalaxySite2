#!/usr/bin/env python3

import sys
import os
import time
import numpy as np 
import gzip
import Galaxy
import libGalaxy
from libGalaxy import DB_PDB_HOME
from Galaxy.core import Template
from Galaxy.utils import get_mmcif

D_CONTACT = 5.0
D_CONTACT_QUERY = 5.0
RCSB_PATH='http://www.rcsb.org/pdb/files/%s.pdb'

SITE_DB_HOME = 'site_db'

clusters = {} 
with open('/home/j2ho/DB/uni_pdb_map/uni_to_pdb_chain.csv','r') as f: 
    lines = f.readlines()
for i, ln in enumerate(lines):
    pdbs = ln.strip().split(',')[1].split(';')
    for elem in pdbs: 
        clusters[elem] = i

label_to_auth = {} 
with open(f'{SITE_DB_HOME}/chain.list','r') as f: 
    lines = f.readlines()
for ln in lines:
    if ln.startswith('#'): continue
    x = ln.strip().split()
    label_to_auth[x[0]] = x[1]

###

class SiteTemplate:
    def __init__(self, pdb_id, ligand_s):
        self.is_valid = True # originally False
        self.pdb_id = pdb_id.split('_')[0]
        self.chain_id = pdb_id.split('_')[1]
        self.ligand_s = ligand_s
        self.ligand_nameonly_s = []
        for lig in ligand_s:
            self.ligand_nameonly_s.append(lig.split('_')[0])
        # 
        self.clu = ''
        self.resol = ''
        self.num_heavy = 0
        #
        self.max_contact_lig = ''
        idx = self.pdb_id[1:3]
        self.org_pdb_fn = Galaxy.core.FilePath('%s/%s/pdb%s.ent'%(DB_PDB_HOME, idx.lower(), self.pdb_id.lower()))
        self.pdb_fn_lig= {}
        self.lig_cntr = {}
    def __repr__(self):
        return f'{self.pdb_id}_{self.chain_id}'
    def get_max_contact_lig(self, lig_name=''):
        from Galaxy.core.vector import distance
        model = Galaxy.core.PDB(self.pdb_fn)[0]
        Rl_s = {}
        for residue in model.get_residues():
            if residue.isAtom(): continue
            #if not residue.resName() in self.ligand_s: continue
            if not residue.resName().strip() in self.ligand_nameonly_s: continue
            if lig_name != '':
                #if residue.resName() != lig_name:
                if residue.resName().strip() != lig_name:
                    continue
            key = '%s_%s_%d'%(residue.resName(), residue.chainID(), residue.resNo())
            Rl_s[key] = residue._R
        if len(Rl_s) == 0:
            self.is_valid = False
            return
        #
        n_contact = []
        for key in Rl_s:
            count = 0
            for R1 in Rl_s[key]: 
                for residue in model.get_residues():
                    if residue.isHetatm(): continue
                    for R2 in residue._R:
                        if distance(R1, R2) < D_CONTACT:
                            count += 1
            if count > 0: 
                n_contact.append((count, key))
        n_contact.sort(reverse=True)
        #
        #print (n_contact)
        self.contact_ligs = n_contact
        if len(n_contact) == 0:
            self.is_valid=False
        else: 
            self.max_contact_lig = n_contact[0][1]
            self.lig_cntr = center(Rl_s[self.max_contact_lig])
            self.is_valid = True
    def check_position(self, pdb):
        from Galaxy.core.vector import distance
        model = Galaxy.core.PDB(self.pdb_fn)[0]
        for residue in model.get_residues():
            if residue.isAtom(): continue
            if not residue.resName().strip() in self.ligand_nameonly_s: continue
            key = '%s_%s_%d'%(residue.resName(), residue.chainID(), residue.resNo())
            if key == self.max_contact_lig:
                Rl = residue._R
                break
        #
        model = pdb[0]
        count = 0
        for R1 in Rl: 
            for residue in model.get_residues():
                if residue.isHetatm(): continue
                for R2 in residue._R:
                    if distance(R1, R2) < D_CONTACT_QUERY:
                        count += 1
        if count == 0:
            self.is_valid = False
    def get_lig_cntr(self, templ_home, lig_name):
        new_pdb_fn = Galaxy.core.FilePath('%s/%s_%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id, lig_name.strip()))
#        model = Galaxy.core.PDB(self.pdb_fn)[0]
        model = Galaxy.core.PDB(new_pdb_fn)[0]
        Rs = []
        for residue in model.get_residues():
            if residue.isAtom(): continue
            key = '%s_%s_%d'%(residue.resName(), residue.chainID(), residue.resNo())
            #if key == self.max_contact_lig:
            if key == lig_name:
                Rs.append(residue._R)
                #Rl_s = residue._R
                break
        #self.lig_cntr[lig_name.strip()] = center(Rl_s)
        #print (lig_name, Rs) 
        coords = np.array(Rs)
        cntr = (np.mean(coords, axis=1))
        return cntr[0]
    def tr_lig_cntr(self):
        T = self.tm.tr[0]
        R = self.tm.tr[1]
        self.lig_cntr = R.dot(self.lig_cntr) + T
    def remove(self, templ_home): 
        self.pdb_fn = '%s/%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id)
        if os.path.exists(self.pdb_fn):
            os.system('rm %s'%self.pdb_fn)
    def write(self, templ_home):
        if not self.org_pdb_fn.status():
            lower_id = self.pdb_id.lower()
            os.system("wget -q -c %s"%(RCSB_PATH%lower_id))
            self.org_pdb_fn = Galaxy.core.FilePath('%s.pdb'%lower_id)
            if not self.org_pdb_fn.status():
                get_mmcif(lower_id, self.chain_id)
                self.chain_id = self.chain_id[-1]
                if not self.org_pdb_fn.status():
                    self.is_valid = False
                    return
        #
        wrt = []
        check_bb = False
        model = Galaxy.core.PDB(self.org_pdb_fn)[0]
        for residue in model.get_residues():
            if residue.isHetatm():
                wrt.append('%s'%residue.write_full())
            else:
                if not residue.check_bb():
                    continue
                if (residue.chainID() == self.chain_id) or (residue.chainID() == ' '):
                    wrt.append('%s'%residue)
                    check_bb = True
#        for ln in open('%s'%self.org_pdb_fn).readlines():
#            if ln.startswith('HETATM'):
#                wrt.append(ln)
#            if ln.startswith('ATOM'):
#                if ln[21] == self.chain_id:
#                    wrt.append(ln)
#                    check_bb=True
#            if ln.startswith('ENDMDL'):
#                break
        if not check_bb:
            self.is_valid = False
            return
        #
        self.pdb_fn = Galaxy.core.FilePath('%s/%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id))
        fout = open('%s'%self.pdb_fn, 'wt')
        fout.writelines(wrt)
        fout.close()
    def rewrite(self, templ_home, lig_name):
        model = Galaxy.core.PDB(self.pdb_fn)[0]
        model.tr(self.tm.tr[0], self.tm.tr[1])
        wrt = []
        for residue in model.get_residues():
            if residue.isAtom():
                if not residue.check_bb():
                    continue
                wrt.append('%s'%residue)
            else:
                key = '%s_%s_%d'%(residue.resName(), residue.chainID(), residue.resNo())
                #if key == self.max_contact_lig:
                if key == lig_name:
                   #wrt.append('%s'%residue)
                   self.num_heavy = len(residue.get_heavy())
                   #wrt.append('%s'%residue.write_full())
                   a=residue.write_full()#.replace("'",'1')
                   a=a.split('\n')
                   at_list = []
                   for i, ln in enumerate(a):
                       at=''.join(ln[12:16].replace("'",' ').split())
                       if "'" in ln[12:16]:
                           if at in at_list:
                               at='%s1'%at
                           a[i] = '%s%-4s%s'%(ln[:12],at,ln[16:])
                       at_list.append(at)
                       a[i]=a[i].replace("'",' ')
                       if len(a[i]) > 78:
                           todel = len(a[i])-78
                           a[i]='%s%s'%(a[i][:66],a[i][66+todel:])
                   a='\n'.join(a)
                   wrt.append('%s'%a)
#                   wrt.append('%s'%residue.write_full().replace("'",' '))
        #
        #self.pdb_fn_lig[lig_name] = Galaxy.core.FilePath('%s/%s_%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id, lig_name))
        #fout = open('%s'%self.pdb_fn_lig[lig_name], 'wt')
        #self.pdb_fn = Galaxy.core.FilePath('%s/%s_%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id, lig_name.strip()))
        new_pdb_fn = Galaxy.core.FilePath('%s/%s_%s_%s.pdb'%(templ_home, self.pdb_id, self.chain_id, lig_name.strip()))
        fout = open('%s'%new_pdb_fn, 'wt')
        fout.writelines(wrt)
        fout.close()
        #

def center(R_s):
    cntr = [0.0, 0.0, 0.0]
    for R in R_s:
        cntr[0] += R[0]
        cntr[1] += R[1]
        cntr[2] += R[2]
    for i in range(3):
        cntr[i] /= len(R_s)
    import numpy as np
    cntr = np.array(cntr)
    return cntr

def read_site_db(search_metal=False):
    DB_path = '%s/ligands.240510'%SITE_DB_HOME
    lines = open(DB_path).readlines()
    #
    db = {}
    for line in lines:
        x = line.strip().split()
        pdb_id = x[0]
        templ = SiteTemplate(pdb_id, x[2].split(','))
        #templ.clu = int(x[1])
        templ.clu = x[1]
        if pdb_id in clusters: 
            templ.clu = clusters[pdb_id] # same uniprot id = same cluster 
        else: 
            templ.clu = '' # if has no uniprot entry, no cluster assigned.. 
        db[pdb_id] = templ
        #try:
        #    templ.resol = float(x[2])
        #except:
        #    templ.resol = -1*100
    return db


def _run_tm(inp):
    tm = Galaxy.utils.TM_align(inp[0],inp[1])
    return tm

def checkyear(pdbid): 
    mmciffile = f'/store/AlphaFold/pdb_mmcif/mmcif_files/{pdbid}.cif'
    if not os.path.exists(mmciffile): 
        return False
    else: 
        for ln in open(mmciffile).readlines(): 
            if '_pdbx_database_status.recvd_initial_deposition_date' in ln: 
                year = ln.split()[1].split('-')[0]
        return int(year)

def foldseek_search(pdb_fn, exclude_pdb, n_proc, re_run, benchmark=None):
    templ_s = []
    cwd = os.getcwd()
    target_dir = os.path.dirname(cwd)
    log = f'{target_dir}/foldseek.log'
    if not os.path.exists(log):
#        os.system(f'/home/j2ho/bin/foldseek/bin/foldseek easy-search {pdb_fn} /home/j2ho/DB/foldseekDB/foldseekDB {log} {target_dir}/tmp -v 0 -e 0.01 --threads 8')
        os.system(f'/home/j2ho/bin/foldseek/bin/foldseek easy-search {pdb_fn} /home/j2ho/DB/foldseekDB/foldseekDB {log} {target_dir}/tmp -e 0.1 --threads 8 --max-seqs 5000')
    f = open(log, 'r') 
    for ln in f.readlines(): 
        x = ln.strip().split()
        if not benchmark == None: 
            pdbid = x[1][:4]
            year = checkyear(pdbid)
            if not year: 
                continue
            elif year >= int(benchmark):
                continue
        pdbid = x[1][:4].upper()
        chain = x[1].split('_')
        if len(chain)<2:
            pdb_id = pdbid 
        else: 
            pdb_id = f'{pdbid}_{chain[-1]}'
        if pdb_id in label_to_auth: 
            pdb_id = label_to_auth[pdb_id]
        fident = float(x[2]) 
        bit = float(x[-1])
        evalue = float(x[-2])
        if evalue == 0: 
            escore = 100000
        else:
            escore = -np.log(evalue) 
        #if fident > 0.85: 
        #    continue
        templ = Template(pdb_id, [1,100], [1,100])  
        templ.score = escore
        #if len(templ_s) < 30: 
        templ_s.append(templ) 
    return templ_s

def mmseq_search(job, exclude_pdb, n_proc, re_run, benchmark=None): 
    templ_s = []
    cwd = os.getcwd()
    target_dir = os.path.dirname(cwd)
    log = f'{target_dir}/mmseq.log'
    mmseqs = '/applic/mmseqs2/15-6f452/bin/mmseqs'
    if not os.path.exists(log):
        #os.system(f'mmseqs easy-search {target_dir}/{job.title}/{job.title}.fa /home/j2ho/DB/mmseqs_pdb_db/pdbdb {log} {target_dir}/tmp --threads 8')
        os.system(f'{mmseqs} easy-search {target_dir}/{job.title}/{job.title}.fa /home/j2ho/DB/mmseqs_pdb_db/pdbdb {log} {target_dir}/tmp -s 9.5 -e 0.1 --threads 8 --max-seqs 5000')
    f = open(log, 'r')
    for ln in f.readlines(): 
        x = ln.strip().split()
        pdbid = x[1][:4]
        mmciffile = f'/store/AlphaFold/pdb_mmcif/mmcif_files/{pdbid}.cif'
        if not os.path.exists(mmciffile): 
            continue
        if not benchmark == None: 
            year = checkyear(pdbid) 
            if year >= int(benchmark):
                continue
        pdbid = x[1][:4].upper()
        chain = x[1].split('_')
        if len(chain)<2:
            pdb_id = pdbid 
        else: 
            pdb_id = f'{pdbid}_{chain[-1]}'
        if pdb_id in label_to_auth: 
            pdb_id = label_to_auth[pdb_id]
        evalue = float(x[-2])
        escore = -np.log(evalue) 
        templ = Template(pdb_id, [1,100], [1,100])  
        templ.score = escore
        templ_s.append(templ)
    
    return templ_s

def search_site_template(job, pdb_fn, fa_fn, search_metal=False, search_method='seq',\
                         exclude_pdb=[], n_proc=None, re_run=False, benchmark=False):
    n_proc = Galaxy.core.define_n_proc(n_proc, multi_node = False)
    out_f_s = {'site.templ_s': Galaxy.core.FilePath('%s.templ_id_s'%job.title)}
    #
    if libGalaxy.BENCHMARK_MODE != 'None':
        exclude_pdb_blast = Galaxy.tools.fr.libfr_utils.get_exclude_pdb(job, fa_fn, n_proc=n_proc)
        for pdb in exclude_pdb_blast:
            if pdb not in exclude_pdb:
                exclude_pdb.append(pdb)
    # 
    if search_method == 'seq':
        templ_s = mmseq_search(job, exclude_pdb, n_proc, re_run, benchmark)

    elif search_method == 'foldseek':
        t1 = time.time()
        templ_s = foldseek_search(pdb_fn, exclude_pdb, n_proc, re_run, benchmark)
        t2 = time.time()
        print (f'foldseek search time: {t2-t1} sec')
    else:
        sys.stdout.write('ERROR: wrong search method : str or seq allowed %s given.\n'%search_method)
    #
    db = read_site_db(search_metal=search_metal)
    #
    filtered = []
    done = []
    for templ in templ_s:
        pdb_id = templ.pdb_id
        if pdb_id not in db: 
            continue
        #if templ.pdb_id not in clu:
        #    continue
        #for pdb_id in clu[templ.pdb_id]:
        #    if pdb_id not in db:
        #        continue
        if pdb_id in done:
            continue
        done.append(pdb_id)
        data = db[pdb_id]
        data.score = templ.score
        filtered.append(data)
    #
    for elem in filtered: 
        print (elem, elem.score) 
    filtered.sort(reverse=True, key=lambda templ: (templ.score,templ.resol))
    return filtered

