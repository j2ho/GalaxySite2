### Not a stand-alone tool YET 
1. needs Galaxy, GalaxyPipe
2. needs openbabel (see "# OpenBabel path" in Site2/site.py)
3. need to set site_db path in "libfr_site.py"
```bash
SITE_DB_HOME = '/home/j2ho/projects/GalaxySite2/Site2/site_db'
```
4. needs Alphafold mmcif pdbs, rcsb pdb DB ... (otherwise will try to download from rcsbpdb using "wget")
5. uses DB/site_sdfs from /home/j2ho (in our cluster)
6. needs foldseek or mmseq (depending on -tsm option / see search functions in Site2/libfr_site.py)
     "/home/j2ho/bin/foldseek/bin/foldseek"

### How to run ..
```bash
$ python3 run_site.py -p {protein PDB file} -tsm foldseek -t result
```
