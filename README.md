### This is distribution for lab members
### Not a stand-alone tool. should be ran in seoklab cluster3  
1. needs Galaxy modules
2. needs openbabel (see "# OpenBabel path" in Site2/site.py)
3. need to set site_db path in "libfr_site.py"
```bash
SITE_DB_HOME = '/path/to/site_db'
```
4. needs Alphafold mmcif pdbs, rcsb pdb DB ... (otherwise will try to download from rcsbpdb using "wget")
5. There are many custom paths and DBs needed: Full mmcif DB, ligand ideal sdf DB. (ask me if you are a lab member) 
6. needs foldseek or mmseq (depending on -tsm option / see search functions in Site2/libfr_site.py) installed

### How to run ..
```bash
$ python3 run_site.py -p {protein PDB file} -tsm foldseek -t result
```
