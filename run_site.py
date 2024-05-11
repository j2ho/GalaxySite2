#!/usr/bin/env python3

import os
import sys
os.environ['GALAXY_HOME']='/home/j2ho/Galaxy'
os.environ['GALAXY_PIPE_HOME']='/home/j2ho/GalaxyPipe'
sys.path.insert(0, '%s/lib'%'/home/j2ho/GalaxyPipe')
import Galaxy
import libGalaxy

def site():
    opt = Galaxy.core.ArgumentParser\
        (description='''GalaxySite: Ligand binding site prediction''')
    #
    opt.add_argument('-p', '--pdb', dest='pdb_fn', metavar='PDB', default=None, \
                         help='Protein structure file(PDB)')
    opt.add_argument('-l', '--ligand', dest='lig_fn', metavar='LIGAND', default=None, \
                         help='Residue name of selected ligand')
    opt.add_argument('-tsm', '--method', dest='search_method', metavar='METHOD', default='seq', \
                         help='Template search using sequence or structure')
    opt.add_argument('-d', '--run_dock', dest='run_dock', metavar='RUN_DOCK', default=False, \
                         help='Run docking on top3 binding sites')
    opt.add_argument('-fp', '--fptype', dest='fptype', metavar='FPTYPE', default='FP2', \
            help='Ligand similarity calculation, choose fingerprint type')
    opt.add_argument('-sim', '--simtype', dest='simtype', metavar='SIMTYPE', default='tani', \
                         help='Ligand similarity calculation, choose similarity metrics')
    opt.add_argument('-benchmark', '--benchmark', dest='benchmark', metavar='BENCHMARK', default=None, \
                         help='run benchmark mode, select year')
    if len(sys.argv) == 1:
        opt.print_help()
        return
    #
    fn = opt.parse_args()
    #
    if fn.pdb_fn == None:
        sys.stdout.write('ERROR')
    #
    if fn.title == None:
        title = Galaxy.core.fn.name(fn.pdb_fn)
    else:
        title = fn.title
    if opt['qsub']:
        opt.que_submit(title=title, log='%s.log.q'%title)
        return
    #
    pdb_fn = None
    if fn.pdb_fn is not None:
        pdb_fn = Galaxy.core.FilePath(fn.pdb_fn)
    if fn.run_dock == 'False':
        fn.run_dock=False
    #
    # Initialize a job
    job = Galaxy.initialize(title=title)
    #
    # GalaxySite
    out_f_s = Galaxy.Site2.run(job, fa_fn=None, pdb_fn=pdb_fn, lig_fn=fn.lig_fn,\
                              search_method=fn.search_method,\
                              run_dock=fn.run_dock, fptype=fn.fptype, simtype=fn.simtype, benchmark=fn.benchmark)
    #
    # Report
#    job.mkdir("result")
#    os.system('cp %s .'%job.pdb_fn)
#    for key in out_f_s:
#        os.system('cp %s .'%out_f_s[key])
#    job.chdir()

if __name__=='__main__':
    site()

