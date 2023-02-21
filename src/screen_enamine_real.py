################################################################################
# Score all Enamine ligands. Points for:
# 1) including pyrimidine (common substructure of binders), 
# 2) being an anion (the known binder, ADPr, is an anion), 
# 3) scaffold_score, including a MCS similarity for global structure and a Tanimoto score for local similarity
################################################################################

import os
import time
import numpy as np
import signal
import logging
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
from pathos.multiprocessing import Pool
import pickle as pkl

from .utils import TimeoutException, timeout_handler

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def check_pyrimidine(ligand):
    """Check if a ligand contains a pyrimidine ring system
    Args:
        ligand (rdkit.Chem.rdchem.Mol): Ligand to check for pyramidine ring system.
    """
    pyrimidine = Chem.MolFromSmiles('c1cncnc1')
    if ligand.HasSubstructMatch(pyrimidine):
        return True
    else:
        return False

def check_anion(ligand):
    """Check if a ligand is an anion
    Args:
        ligand (rdkit.Chem.rdchem.Mol): Ligand to check for anion.
    """
    if Chem.rdmolops.GetFormalCharge(ligand) < 0:
        return True
    else:
        return False

def num_atoms_overlapping(ligand_scaffold, fragment_scaffold):
    """Compare the maximum common substructure between a ligand scaffold and a fragment scaffold and return the number of atoms in the MCS.
    Args:
        ligand (rdkit.Chem.rdchem.Mol): Ligand to check for fragment scaffold match.
        fragment_scaffold (rdkit.Chem.rdchem.Mol): Fragment scaffold to check for.
    Returns:
        float: Number of atoms in the maximum common substructure.
    """
    mcs = rdFMCS.FindMCS([ligand_scaffold, fragment_scaffold], completeRingsOnly=True, ringMatchesRingOnly=True)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    return len(mcs_mol.GetAtoms())

def calc_scaffold_score(ligand, fragments, fragment_fingerprints, fragment_scaffolds, tanimoto_weight=0.5, mcs_weight=0.5):
    """Calculate the scaffold score for a ligand, comprising an MCS score for global structure and a Tanimoto score for local similarity.
    Args:
        ligand (rdkit.Chem.rdchem.Mol): Ligand to calculate scaffold score for.
        fragments (list of rdkit.Chem.rdchem.Mol): List of fragments.
        fragment_fingerprints (list of rdkit.DataStructs.cDataStructs.ExplicitBitVect): List of fragment fingerprints.
        fragment_scaffolds (list of rdkit.Chem.rdchem.Mol): List of fragment scaffolds.
    Returns:
        float: Scaffold score.
    """
    # Calculate Tanimoto scores
    ligand_fingerprint = Chem.RDKFingerprint(ligand, fpSize=1024)
    tanimoto_scores = DataStructs.BulkTanimotoSimilarity(ligand_fingerprint, fragment_fingerprints)

    # Calculate MCS scores
    ligand_scaffold = GetScaffoldForMol(ligand)
    mcs_scores = []
    for fragment_scaffold in fragment_scaffolds:
        mcs_scores.append(num_atoms_overlapping(ligand_scaffold, fragment_scaffold))

    scaffold_scores = [mcs_score * mcs_weight + tanimoto_score * tanimoto_weight for mcs_score, tanimoto_score in zip(mcs_scores, tanimoto_scores)]
    # Get argmax of scaffold_scores
    max_score_idx = np.argmax(scaffold_scores)

    return fragments[max_score_idx], scaffold_scores[max_score_idx]


def compute_fragment_fingerprints_scaffolds(fragment_dir, merged_fragment_dir):
    """Compute fragment fingerprints and scaffolds.
    Args:
        ligand_hits_dir (str): Directory containing ligand hits.
    Returns:
        fragments (list of rdkit.Chem.rdchem.Mol): List of fragments.
        fragment_fingerprints (list of rdkit.DataStructs.cDataStructs.ExplicitBitVect): List of fragment fingerprints.
        fragment_scaffolds (list of rdkit.Chem.rdchem.Mol): List of fragment scaffolds.
    """
    # Precompute fragment fingerprints and scaffolds
    fragments = []
    fragment_fingerprints = []
    fragment_scaffolds = []

    for fragment_path in sorted(os.listdir(fragment_dir)):
        if os.path.isdir(os.path.join(fragment_dir, fragment_path)):
            continue
        # Load fragment from SDF
        fragment = Chem.SDMolSupplier(os.path.join(fragment_dir, fragment_path))[0]
        fragment.SetProp('_Name', fragment_path.split('.')[0])
        fragments.append(fragment)
        # Get morgan fingerprint
        fragment_fingerprints.append(Chem.RDKFingerprint(fragment, fpSize=1024))
        # Get Murcko scaffold
        fragment_scaffolds.append(GetScaffoldForMol(fragment))

    for fragment_path in sorted(os.listdir(merged_fragment_dir)):
        if os.path.isdir(os.path.join(merged_fragment_dir, fragment_path)):
            continue
        fragment = Chem.SDMolSupplier(os.path.join(merged_fragment_dir, fragment_path))[0]
        fragment.SetProp('_Name', fragment_path.split('.')[0])
        fragments.append(fragment)
        # Get morgan fingerprint
        fragment_fingerprints.append(Chem.RDKFingerprint(fragment, fpSize=1024))
        # Get Murcko scaffold
        fragment_scaffolds.append(GetScaffoldForMol(fragment))

    # De-dupe fragments
    fragments_deduped = []
    fragment_fingerprints_deduped = []
    fragment_scaffolds_deduped = []
    for fragment, fragment_fingerprint, fragment_scaffold in zip(fragments, fragment_fingerprints, fragment_scaffolds):
        if fragment_fingerprint not in fragment_fingerprints_deduped:
            fragments_deduped.append(fragment)
            fragment_fingerprints_deduped.append(fragment_fingerprint)
            fragment_scaffolds_deduped.append(fragment_scaffold)

    log_str = f'Number of fragments before deduplication: {len(fragments)}'
    logging.info(log_str)
    #log_cmd = f'echo "{log_str}" >> log.txt'
    #os.system(log_cmd)
    log_str = f"Number of fragments after deduplication: {len(fragments_deduped)}"
    logging.info(log_str)
    #log_cmd = f'echo "{log_str}" >> log.txt'
    #os.system(log_cmd)
    
    return fragments_deduped, fragment_fingerprints_deduped, fragment_scaffolds_deduped

def score_one_ligand(line, tanimoto_weight, mcs_weight, scaffold_score_weight, anion_score_weight, pyrimidine_score_weight, fragments, fragment_fingerprints, fragment_scaffolds):
    """Score one ligand.
    Args:
        line (str): Line from ligand hits file.
        tanimoto_weight (float): Weight of Tanimoto score in scaffold score.
        mcs_weight (float): Weight of MCS score in scaffold score.
        scaffold_score_weight (float): Weight of scaffold score in ligand score.
        anion_score_weight (float): Weight of anion score in ligand score.
        pyrimidine_score_weight (float): Weight of pyrimidine score in ligand score.
        fragments (list of rdkit.Chem.rdchem.Mol): List of fragments.
        fragment_fingerprints (list of rdkit.DataStructs.cDataStructs.ExplicitBitVect): List of fragment fingerprints.
        fragment_scaffolds (list of rdkit.Chem.rdchem.Mol): List of fragment scaffolds.
    Returns:
        ligand_name (str): Name of ligand.
        ligand_smiles (str): SMILES of ligand.
        base_fragment_smiles (str): SMILES of base fragment.
        ligand_score (float): Ligand score.
    """
    ligand_name = ""
    ligand_smiles = ""
    base_fragment_smiles = ""
    ligand_score = 0.0
    
    # Set a timeout for the merge
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(4)  # set the alarm for 4 seconds
    
    try:
        # Split line on tab
        cxsmiles, ligand_name, size = line.split('\t')
        ligand = Chem.MolFromSmiles(cxsmiles)
        ligand_smiles = Chem.MolToSmiles(ligand)
        base_fragment, scaffold_score = calc_scaffold_score(ligand=ligand, fragments=fragments, fragment_fingerprints=fragment_fingerprints, fragment_scaffolds=fragment_scaffolds, tanimoto_weight=tanimoto_weight, mcs_weight=mcs_weight)
        base_fragment_smiles = Chem.MolToSmiles(base_fragment)
        anion_score = check_anion(ligand)
        pyrimidine_score = check_pyrimidine(ligand)
        ligand_score = scaffold_score * scaffold_score_weight + anion_score * anion_score_weight + pyrimidine_score * pyrimidine_score_weight
        
    except TimeoutException as e:
        ligand_name = ligand_name + '_timeout'
        pass

    # reset the timeout
    signal.alarm(0)

    return ligand_name, ligand_smiles, base_fragment_smiles, ligand_score


def score_ligands_in_parallel(enamine_real_database_path, fragments, fragment_fingerprints, fragment_scaffolds, num_processes, scores_path):
    """Score all Enamine REAL ligands. Place results in a dataframe.
    Args:
        enamine_real_database_dir (str):  to Enamine REAL database.
        fragments (list of rdkit.Chem.rdchem.Mol): List of fragments.
        fragment_fingerprints (list of rdkit.DataStructs.cDataStructs.ExplicitBitVect): List of fragment fingerprints.
        fragment_scaffolds (list of rdkit.Chem.rdchem.Mol): List of fragment scaffolds.
        num_processes (int): Number of processes to use.
        scores_path (str): Path to save scores.
    """
    # Sync results every X ligands
    sync_every = 200 
    save_every = sync_every * 1
    # Score all ligands
    mcs_weight = 0.5
    tanimoto_weight = 0.5
    scaffold_score_weight = 1.
    # Choose these weights based on the range of scaffold scores observed
    anion_score_weight = 3.
    pyrimidine_score_weight = 3.
    ligand_names = []
    ligand_smiles_list = []
    base_fragment_smiles_list = []
    ligand_scores = []

    log_path = scores_path.replace('.pkl', '.log')
    # Check if checkpoint file exists
    if os.path.isfile(scores_path):
        with open(scores_path, 'rb') as f:
            checkpoint = pkl.load(f)
        ligand_names = checkpoint['ligand_names']
        ligand_smiles_list = checkpoint['ligand_smiles']
        base_fragment_smiles_list = checkpoint['base_fragment_smiles']
        ligand_scores = checkpoint['ligand_scores']
        

    with open(enamine_real_database_path, 'r') as f:
        # Skip first line
        next(f)
        
        # Skip lines that have already been processed
        for i in range(len(ligand_names)):
            next(f)
        # Define a function to score a single ligand
        def score_one_ligand_wrapper(line):
            return score_one_ligand(line=line, tanimoto_weight=tanimoto_weight, mcs_weight=mcs_weight, scaffold_score_weight=scaffold_score_weight, anion_score_weight=anion_score_weight, pyrimidine_score_weight=pyrimidine_score_weight, fragments=fragments, fragment_fingerprints=fragment_fingerprints, fragment_scaffolds=fragment_scaffolds)

        # Use a pool of processes to score ligands in parallel
        with Pool(num_processes) as pool:
            results = pool.imap_unordered(score_one_ligand_wrapper, f, chunksize=sync_every)

            batch_start_time = time.time()
            # iter_start_time = time.time()
            for i, (ligand_name, ligand_smiles, base_fragment_smiles, ligand_score) in enumerate(results):
                ligand_names.append(ligand_name)
                ligand_smiles_list.append(ligand_smiles)
                base_fragment_smiles_list.append(base_fragment_smiles)
                ligand_scores.append(ligand_score)
                # print("Iter time:", time.time() - iter_start_time)
                # Save checkpoint
                if i % (sync_every) == 0:
                    log_str = f"Finished batch {len(ligand_names) // sync_every}. Number of ligands scored: {len(ligand_names)}. Time for this batch: {time.time() - batch_start_time}"
                    logging.info(log_str)

                    if i % (num_processes * save_every) == 0:
                        checkpoint = {'ligand_names': ligand_names, 'ligand_smiles': ligand_smiles_list, 'base_fragment_smiles': base_fragment_smiles_list, 'ligand_scores': ligand_scores}
                        with open(scores_path, 'wb') as f:
                            pkl.dump(checkpoint, f)
                        log_str = "Checkpoint saved. Starting new batch"
                        logging.info(log_str)
                        #log_cmd = f'echo "{log_str}" >> log.txt'
                        #os.system(log_cmd)

                    batch_start_time = time.time()
                
                # iter_start_time = time.time()



def screen_enamine_real(enamine_real_database_path, fragment_dir, merged_fragment_dir, num_processes, scores_path):
    """Score all Enamine REAL ligands. Place results in a dataframe.
    Args:
        enamine_real_database_dir (str): Directory to Enamine REAL database.
        ligand_hits_dir (str): Directory containing fragment hits.
        results_path (str): Path to save results to.
    """
    # Precompute fragment fingerprints and scaffolds
    fragments, fragment_fingerprints, fragment_scaffolds = compute_fragment_fingerprints_scaffolds(fragment_dir, merged_fragment_dir)

    score_ligands_in_parallel(enamine_real_database_path, fragments, fragment_fingerprints, fragment_scaffolds, num_processes=num_processes, scores_path=scores_path)
