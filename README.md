### Not a stand-alone tool YET 
1. needs Galaxy, GalaxyPipe
2. needs openbabel (see "# OpenBabel path" in Site2/site.py) 
3. needs Alphafold mmcif pdbs, rcsb pdb DB ... (otherwise will try to download from rcsbpdb using "wget")
4. uses DB/site_sdfs from /home/j2ho (in our cluster)
5. needs foldseek or mmseq (depending on -tsm option / see search functions in Site2/libfr_site.py)
     "/home/j2ho/bin/foldseek/bin/foldseek" 

### How to run ..
```bash
$ python3 run_site.py -p {protein PDB file} -tsm foldseek -t result
```
