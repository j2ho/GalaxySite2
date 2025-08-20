#!/usr/bin/env python3

import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Union

import Galaxy
import libGalaxy
from Galaxy.core import FilePath, define_n_proc, file_status
from . import site
from .prep import prepare_lig_mol2
from .bsite_pred import predict_bsite

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GalaxySiteError(Exception):
    pass


def run(job, fa_fn: Optional[str] = None, pdb_fn: Optional[str] = None, 
        n_proc: Optional[int] = None, re_run: bool = False, **kwargs) -> Dict[str, Any]:
    """
    Run GalaxySite 
    
    Args:
        job: Job object from run_site.py 
        fa_fn: Path to FASTA file (optional)
        pdb_fn: Path to PDB file (optional)
        n_proc: Number of processes to use
        re_run: Whether to re-run existing results
        **kwargs: Additional keyword arguments
        
    Returns:
        Dictionary of outputs
        
    Raises:
        GalaxySiteError: If neither PDB nor FASTA file is provided, or if PDB is required but not provided
    """
    if pdb_fn is None and fa_fn is None:
        raise GalaxySiteError("Either pdb_fn or fa_fn must be provided")
    
    define_n_proc(n_proc)
    
    if pdb_fn is None and fa_fn is not None:
        error_msg = f'ERROR: Cannot run GalaxySite using sequence only: {fa_fn}'
        logger.error(error_msg)
        raise GalaxySiteError(error_msg)
    
    if pdb_fn is not None:
        pdb_fn = _ensure_filepath(pdb_fn)
        logger.info(f'Processing GalaxySite using structure: {pdb_fn}')
        
        try:
            pdb = Galaxy.core.PDB(pdb_fn)
        except Exception as e:
            raise GalaxySiteError(f"Failed to load PDB file {pdb_fn}: {e}")
        
        if fa_fn is None:
            fa_fn = FilePath(f'{job.title}.fa')
            try:
                pdb.write_fasta_file(fa_fn)
                logger.info(f'Generated FASTA file: {fa_fn}')
            except Exception as e:
                raise GalaxySiteError(f"Failed to generate FASTA from PDB: {e}")
        else:
            fa_fn = _ensure_filepath(fa_fn)
    
    job.pdb_fn = pdb.pdb_fn
    job.fa_fn = fa_fn
    
    try:
        out_f_s = run_site(job, pdb, re_run=re_run, **kwargs)
        logger.info("GalaxySite analysis completed successfully")
        return out_f_s
    except Exception as e:
        logger.error(f"GalaxySite analysis failed: {e}")
        raise


def run_site(job, pdb, re_run: bool = False, **kwargs) -> Dict[str, FilePath]:
    """
    Running actual prediction process (site.py) 
    
    Args:
        job: Job object 
        pdb: PDB object containing protein structure
        re_run: Whether to re-run existing results
        **kwargs: Additional keyword arguments
        
    Returns:
        Dictionary mapping output file keys to FilePath objects
    """
    try:
        ligand_s = site.select_ligand(job, pdb, re_run=re_run, **kwargs)
        logger.info(f"Selected {len(ligand_s)} ligands for analysis")
    except Exception as e:
        raise GalaxySiteError(f"Failed to select ligands: {e}")
    
    if kwargs.get('run_dock', True) is False:
        logger.info("Docking disabled, returning early")
        return {}
    
    out_f_s = {}
    
    for ligand in ligand_s:
        if not ligand.selected:
            continue
            
        ligand_name = ligand.lig_name
        logger.info(f"Processing ligand: {ligand_name}")
        
        try:
            job.mkdir(ligand_name, tag=f'{ligand_name}.HOME')
            
            site_pdb_key = f'{ligand_name}.site.pdb'
            bsite_dat_key = f'{ligand_name}.bsite.dat'
            
            out_f_s[site_pdb_key] = FilePath(f'{ligand_name}.site.pdb')
            out_f_s[bsite_dat_key] = FilePath(f'{ligand_name}.bsite.dat')
            
            # Check if results already exist and skip if not re-running
            if not re_run and file_status(out_f_s):
                logger.info(f'Existing results found for ligand {ligand_name}, reusing...')
                job.chdir_prev()
                continue
            
            # Process ligand
            _process_single_ligand(job, pdb, ligand, re_run)
            
        except Exception as e:
            logger.error(f"Failed to process ligand {ligand_name}: {e}")
            job.chdir_prev()
            raise GalaxySiteError(f"Ligand processing failed for {ligand_name}: {e}")
        
        finally:
            job.chdir_prev()
    
    job.update(out_f_s)
    return out_f_s


def _process_single_ligand(job, pdb, ligand, re_run: bool) -> None:
    """
    Process a single ligand through the pipeline
    
    Args:
        job: Job object 
        pdb: PDB object containing protein structure
        ligand: Ligand object 
        re_run: Whether to re-run existing results
    """
    ligand_name = ligand.lig_name
    
    mol2_fn_no_h = prepare_lig_mol2(job, ligand, exclude_H=True)
    ligand.lig_fn = mol2_fn_no_h
    
    rsr_fn, rsr_pdb = site.extract_rsr(job, pdb, ligand)
    
    mol2_fn_with_h = prepare_lig_mol2(job, ligand, exclude_H=False)
    
    opt_s = _get_docking_options(rsr_pdb, rsr_fn)
    
    optimized, opt_energy = Galaxy.galaxy.ligdock(
        job, pdb, mol2_fn_with_h, 
        outfile_prefix=ligand_name,
        lig_name=ligand_name, 
        opt_s=opt_s, 
        re_run=re_run
    )
    
    bsite_s = predict_bsite(pdb[0], optimized[0])
    
    _write_output_files(ligand_name, optimized[0], bsite_s)
    
    logger.info(f"Successfully processed ligand: {ligand_name}")


def _get_docking_options(rsr_pdb, rsr_fn) -> Dict[str, str]:
    """
    Get docking configuration options.
    
    Args:
        rsr_pdb: Receptor PDB file
        rsr_fn: Restraints file
        
    Returns:
        Dictionary of docking options
    """
    return {
        'weight_type': 'GalaxySite',
        'rsr_pdb': rsr_pdb.relpath(),
        'rsr_file': rsr_fn.relpath(),
        'soften_dock_E': 'yes',
        'include_initial_conf': 'yes',
        'csa_bank': '30 30',
        'csa_seed': '2  15  100  2',
        'csa_n_opr_s': '2  3  2  3  0  0',
    }


def _write_output_files(ligand_name: str, optimized_structure, bsite_data) -> None:
    """
    Write output files
    
    Args:
        ligand_name: Name of the ligand
        optimized_structure: Optimized structure
        bsite_data: Binding site data
    """
    site_pdb_filename = f'{ligand_name}.site.pdb'
    try:
        with open(site_pdb_filename, 'w') as fout:
            fout.writelines(optimized_structure.write_in_PDB())
        logger.debug(f"Written structure file: {site_pdb_filename}")
    except Exception as e:
        raise GalaxySiteError(f"Failed to write structure file {site_pdb_filename}: {e}")
    
    bsite_dat_filename = f'{ligand_name}.bsite.dat'
    try:
        with open(bsite_dat_filename, 'w') as fout:
            fout.writelines(bsite_data)
        logger.debug(f"Written binding site file: {bsite_dat_filename}")
    except Exception as e:
        raise GalaxySiteError(f"Failed to write binding site file {bsite_dat_filename}: {e}")


def _ensure_filepath(filepath: Union[str, FilePath]) -> FilePath:
    """
    Ensure input is a FilePath object.
    
    Args:
        filepath: String path or FilePath object
        
    Returns:
        FilePath object
    """
    if not isinstance(filepath, FilePath):
        return FilePath(filepath)
    return filepath


def validate_inputs(fa_fn: Optional[str] = None, pdb_fn: Optional[str] = None) -> None:
    """
    Validate input file paths exist and are accessible.
    
    Args:
        fa_fn: Path to FASTA file
        pdb_fn: Path to PDB file
        
    Raises:
        GalaxySiteError: If required files don't exist or aren't accessible
    """
    if pdb_fn and not os.path.exists(pdb_fn):
        raise GalaxySiteError(f"PDB file not found: {pdb_fn}")
    
    if fa_fn and not os.path.exists(fa_fn):
        raise GalaxySiteError(f"FASTA file not found: {fa_fn}")


if __name__ == "__main__":
    pass
