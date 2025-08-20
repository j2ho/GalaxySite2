#!/usr/bin/env python3

# Database paths in the repository
SITE_DB_HOME = '/your/path/to/Site2/site_db' 

# Paths to external databases and tools
SDF_DB = '/your/path/to/ligand/sdf/db' # should have all sdf files you want to search
RCSB_mmcif_DB = '/your/path/to/rcsb/mmcif_files' # should have all mmcif files you want to search
mmseqs = '/your/path/to/mmseqs'
mmseqs_db = '/your/path/to/mmseqsDB'
foldseek = '/your/path/to/foldseek'
foldseek_db = '/your/path/to/foldseekDB'

# Tool paths
BABEL = '/your/path/to/obabel'

# Distance thresholds
D_CONTACT = 5.0
D_CONTACT_QUERY = 5.0

# TM cutoff for template selection
TM_CUT = 0.2

# RCSB URL pattern
RCSB_PATH = 'http://www.rcsb.org/pdb/files/%s.pdb'

# Fingerprint options
finger_print = ['-xfFP3', '-xfFP4', '-xfMACCS']

# Ligand type dictionary
type_dict = {1:'lipid',2:'non-biological',3:'ions/metal',4:'glycan',0:'general'}