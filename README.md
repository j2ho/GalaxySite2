### This is distribution for lab members
### Not a stand-alone tool. Designed for seoklab cluster3  


**Setup Requirements:**

1. **Configure paths in config.py** - Update `Site2/config.py` with your specific paths:
   - `SITE_DB_HOME`: Path to site database
   - `SDF_DB`: Path to ligand SDF database
   - `RCSB_mmcif_DB`: Path to RCSB mmCIF files
   - `BABEL`: Path to OpenBabel executable
   - `mmseqs`/`foldseek`: Tool paths and database paths

2. **Dependencies:**
   - Galaxy and GalaxyPipe packages
   - OpenBabel
   - pdb files are downloaded on the fly from RCSB PDB database, mmcif is used for info parsing mainly
   - Custom databases: Full mmCIF DB, ligand ideal SDF DB (SDF can be downloaded on the fly as well, but directory should be set in `config.py`)
   - foldseek and/or mmseqs2 (depending on -tsm option)

### How to run ..
```bash
$ python3 run_site.py -p {protein PDB file} -tsm foldseek -t result
```
