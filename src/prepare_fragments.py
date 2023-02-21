import os
from Bio.PDB import PDBParser, PDBIO, Superimposer, Selection
from rdkit import Chem

from fragmenstein import Victor, Monster
import pyrosetta
import signal
import glob

from .utils import TimeoutException, timeout_handler

AMINO_ACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def clean_pdb(in_dir, out_dir):
    """Some PDB files include ligands with residue names that are not 'LIG'. Clean PDB files by 
    renaming these to 'LIG'. Some ligands include 'HOH' and 'SO4' in the residue name. Leave these 
    as they are.
    Args:
        in_dir (str): Directory containing PDB files to clean.
        out_dir (str): Directory to save cleaned PDB files to.
    """
    for pdb_path in sorted(os.listdir(in_dir)):
        in_path = os.path.join(in_dir, pdb_path)
        out_path = out_dir + pdb_path.split('/')[-1]
        if not os.path.exists(out_path):
            with open(in_path, 'r') as f:
                lines = f.readlines()
            with open(out_path, 'w') as f:
                prev_line = ''
                for line in lines:
                    if (line.startswith('HETATM') or prev_line.startswith('HETATM')) and 'HOH' not in line and 'SO4' not in line:
                        line = line[:17] + 'LIG' + line[20:]
                    f.write(line)
                    prev_line = line

def get_atom(res):
    """Get an atom from a residue, preferring CA if it exists.
    Args:
        res (Bio.PDB.Residue): Residue to get atom from.
    Returns:
        Bio.PDB.Atom: Atom from residue.
    """
    if 'CA' in res:
        return res['CA']
    else:
        return res.get_list()[0]

def get_protein_residues(pdb_structure, chain):
    """Get the protein residues from a PDB structure.
    Args:
        pdb_structure (Bio.PDB.Structure): PDB structure to get protein residues from.
        chain (str): Chain to get protein residues from.
    Returns:
        list: List of Bio.PDB.Residue objects.
    """
    protein_residues = []
    for res in pdb_structure[0][chain]:
        if res.get_resname() in AMINO_ACIDS:
            protein_residues.append(res)
    return protein_residues

def align_pdb_to_ground_truth(in_dir, out_dir, ground_truth_path):
    """Align a directory of PDB files to a ground truth PDB file using BioPython and save the aligned PDB files to a new file.

    Args:
        pdb_dir (str): Path to the PDB files to align.
        ground_truth_path (str): Path to the ground truth PDB file.
        out_dir (str): Path to save the aligned PDB files to.
    """
    ground_truth_structure = PDBParser().get_structure('ground_truth', ground_truth_path)
    # Get residue objects from chain A of the ground truth structure
    ground_truth_protein_residues = get_protein_residues(ground_truth_structure, 'A')
    ground_truth_residue_nums = set([res.get_id()[1] for res in ground_truth_protein_residues])
    
    for pdb_path in sorted(os.listdir(in_dir)):
        in_path = os.path.join(in_dir, pdb_path)
        out_path = out_dir + pdb_path.split('/')[-1]

        if not os.path.exists(out_path):
            
            # Mac1 naming convention tells us whether the ligand is bound to chain A or B
            if '0B' in in_path:
                chain = 'B'
            elif '1B' in in_path:
                chain = 'B'
            elif '2B' in in_path:
                chain = 'B'
            else:
                chain = 'A'
            
            pdb_structure = PDBParser().get_structure('pdb', in_path)
            # Get protein atoms for which both proteins include that residue
            pdb_protein_residues = get_protein_residues(pdb_structure, chain)
            pdb_residue_nums = set([res.get_id()[1] for res in pdb_protein_residues])
            residues_in_common = ground_truth_residue_nums.intersection(pdb_residue_nums)
            pdb_protein_atoms = [get_atom(res) for res in pdb_protein_residues if res.get_id()[1] in residues_in_common]
            ground_truth_protein_atoms = [get_atom(res) for res in ground_truth_protein_residues if res.get_id()[1] in residues_in_common]
            # Align PDB to ground truth
            sup = Superimposer()
            # Align protein atoms
            sup.set_atoms(ground_truth_protein_atoms, pdb_protein_atoms)
            # Apply alignment to all atoms
            sup.apply(list(pdb_structure.get_atoms()))
            # Save aligned PDB
            io = PDBIO()
            io.set_structure(pdb_structure)
            io.save(out_path)

def extract_sdf(in_dir, out_dir, batch_1_pdb_dir, batch_1_sdf_dir, batch_2_sdf_path):
    """Transfer coordinates from aligned PDB file to mol object containing proper bonds.
    Args:
        in_dir (str): Directory containing aligned PDB files.
        out_dir (str): Directory to save SDF files to.
        sdf_path (str): Path to SDF file.
    """
    # Batch 1
    for pdb_cleaned_file in sorted(os.listdir(batch_1_pdb_dir)):
        # 1: Extract SDFs from 5 PDB_cleaned files that already have bonds
        if 'bound' in pdb_cleaned_file:
            continue
        # Get PDB file name
        pdb_cleaned_path = os.path.join(batch_1_pdb_dir, pdb_cleaned_file)
        # Get SDF file name
        basename = pdb_cleaned_path.split('/')[-1].split('.')[0]
        sdf_path = os.path.join(batch_1_sdf_dir, basename + '.sdf')
        if not os.path.exists(sdf_path):
            # Extract ligand from PDB file.
            mol_correct_bonds = Victor.extract_mol(name='tmp', filepath=pdb_cleaned_path)
            # Save SDF file
            Chem.SDWriter(sdf_path).write(mol_correct_bonds)

        # 2: Transfer coordinates from PDB_aligned file to SDF file
        pdb_aligned_file = pdb_cleaned_file
        pdb_aligned_path = os.path.join(in_dir, pdb_aligned_file)
        out_path = os.path.join(out_dir, basename + '.sdf')
        failed_path = os.path.join(out_dir, 'failed', basename + '.sdf')
        if not os.path.exists(out_path) and not os.path.exists(failed_path):
            # Extract ligand from PDB file. This has the correct coordinates
            pdb_mol = Victor.extract_mol(name='tmp', filepath=pdb_aligned_path)
            # Transfer coordinates from PDB to SDF
            # create a new conformer for mol2 and copy coordinates from mol1
            conf = Chem.Conformer(mol_correct_bonds.GetNumAtoms())
            correct_conf = pdb_mol.GetConformer()
            for i in range(mol_correct_bonds.GetNumAtoms()):
                try:
                    correct_pos = correct_conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, correct_pos)
                except:
                    with open(failed_path, 'w') as f:
                        f.write('')
                    continue

            mol_correct_bonds.RemoveConformer(0)
            mol_correct_bonds.AddConformer(conf)
            # Save as SDF
            Chem.SDWriter(out_path).write(mol_correct_bonds)


    # Batch 2
    # Load SDF file
    suppl = Chem.SDMolSupplier(batch_2_sdf_path)
    for mol in suppl:
        # Get PDB file name
        pdb_name = mol.GetProp('_Name')
        pdb_path = os.path.join(in_dir, pdb_name + '_bound.pdb')

        # Save the ligand to a PDB file using rdkit. Remove _bound from the end of the original file name if it is there.
        basename = pdb_path.split('/')[-1].split('.')[0].split('_bound')[0]
        out_path = os.path.join(out_dir, basename + '.sdf')
        failed_path = os.path.join(out_dir, 'failed', basename + '.sdf')

        # Save SDF file
        if not os.path.exists(out_path) and not os.path.exists(failed_path):
            # Extract ligand from PDB file. This has the correct coordinates
            pdb_mol = Victor.extract_mol(name='tmp', filepath=pdb_path)
            # Transfer coordinates from PDB to SDF
            # create a new conformer for mol2 and copy coordinates from mol1
            conf = Chem.Conformer(mol.GetNumAtoms())
            correct_conf = pdb_mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                try:
                    correct_pos = correct_conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, correct_pos)
                except:
                    with open(failed_path, 'w') as f:
                        f.write('')
                    continue

            mol.RemoveConformer(0)
            mol.AddConformer(conf)
            # Save as SDF
            Chem.SDWriter(out_path).write(mol)
    

# def cluster_ligands(ligand_paths):
#     """Cluster ligands spatially into four clusters using rdkit and sklearn

#     Args:
#         ligand_paths (list): List of paths to ligand PDB files.

#     Returns:
#         list: List of lists of ligand PDB files that are in the same cluster.
#     """
#     ligands = []
#     for ligand_path in ligand_paths:
#         ligand = Chem.MolFromPDBFile(ligand_path)
#         ligands.append(ligand)
#     ligand_coords = []
#     for ligand in ligands:
#         # Get center of mass of ligand
#         ligand_com = AllChem.ComputeCentroid(ligand)
#         ligand_coords.append(ligand_com)

#     kmeans = KMeans(n_clusters=4, random_state=0).fit(ligand_com)
#     clusters = [[] for i in range(4)]
#     for i, cluster in enumerate(kmeans.labels_):
#         clusters[cluster].append(ligand_paths[i])
#     return clusters

def fragmenstein_monster_merge_fragments(fragment_1, fragment_2):
    """Use Fragmenstein Monster to merge two fragments together.

    Args:
        fragment_1_path (str): Path to the first fragment.
        fragment_2_path (str): Path to the second fragment.

    Returns:
        rdkit.Chem.rdchem.Mol: Merged ligand.
    """
    monster = Monster(hits=[fragment_1, fragment_2])
    monster.combine()
    merged_mol = monster.positioned_mol
    monster.place(mol=merged_mol, merging_mode='none_permissive')
    merged_mol = monster.positioned_mol

    return merged_mol

def fragmenstein_victor_merge_fragments(fragment_1, fragment_2, reference_pdb_path):
    """Use Fragmenstein Victor to merge two fragments together.

    Args:
        fragment_1 (rdkit.Chem.rdchem.Mol): First fragment.
        fragment_2 (rdkit.Chem.rdchem.Mol): Second fragment.
        reference_pdb_path (str): Path to the reference PDB file.

    Returns:
        rdkit.Chem.rdchem.Mol: Merged ligand.
    """

    pyrosetta.init( extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')

    victor = Victor(hits=[fragment_1, fragment_2],
                    pdb_filename=reference_pdb_path,
                    covalent_resi=1) # if not covalent, just put the first residue or something.
    victor.combine()
    merged_mol = victor.positioned_mol
    victor.place(mol=merged_mol, merging_mode='none_permissive')
    merged_mol = victor.minimized_mol

    return merged_mol

def merge_fragments(fragment_dir, merged_fragment_dir, reference_pdb_path):
    """Use Fragmenstein to generate merged fragments by pairwise merging each pair of fragments

    Args:
        fragment_dir (str): Directory containing fragment PDB files.
        merged_light_dir (str): Directory to save the merged fragments to.
        reference_pdb_path (str): Path to the reference PDB file.
    """
    fragment_paths = sorted(glob.glob(fragment_dir + '*.sdf'))
    success_cnt = 0
    fail_cnt = 0
    for i, fragment_1_path in enumerate(fragment_paths):
        print("Fragment: ", i, "/", len(fragment_paths), "Success: ", success_cnt, "Fail: ", fail_cnt)

        # Now merge fragment_1 with all other fragments
        for fragment_2_path in fragment_paths[i+1:]:
            fragment1_name = fragment_1_path.split('/')[-1].split('.')[0]
            fragment2_name = fragment_2_path.split('/')[-1].split('.')[0]
            merged_fragment_name = fragment1_name + '_' + fragment2_name + '_merged.sdf'
            # Write merged fragment to PDB file with RDKit
            save_path = merged_fragment_dir + merged_fragment_name
            merge_timeout_save_path = os.path.join(merged_fragment_dir, 'error_timeout', merged_fragment_name)
            merge_connection_err_save_path = os.path.join(merged_fragment_dir, 'error_connection', merged_fragment_name)

            if not os.path.exists(save_path) and not os.path.exists(merge_timeout_save_path) and not os.path.exists(merge_connection_err_save_path):
                print('Merging: ', fragment1_name, fragment2_name)
                # Load fragment from sdf file to rdkit
                fragment_1 = Chem.SDMolSupplier(fragment_1_path)[0]
                fragment_2 = Chem.SDMolSupplier(fragment_2_path)[0]
                # Set a timeout for the merge
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(1)  # set the alarm for 1 second
                try:
                    # Try merging, but time out if the merge takes more than 1 second
                    merged_fragment = fragmenstein_monster_merge_fragments(fragment_1, fragment_2)
                    # merged_fragment = fragmenstein_victor_merge_fragments(fragment_1, fragment_2, reference_pdb_path)
                    success_cnt += 1
                    
                except TimeoutException as e:
                    fail_cnt += 1
                    # Write empty file to failed filepath
                    with open(merge_timeout_save_path, 'w') as f:
                        f.write('')
                    # reset the timeout
                    signal.alarm(0)
                    continue
                except (ConnectionError, AssertionError, Exception) as e:
                    fail_cnt += 1
                    # Write empty file to failed filepath
                    with open(merge_connection_err_save_path, 'w') as f:
                        f.write('')
                    # reset the timeout
                    signal.alarm(0)
                    continue

                # reset the timeout
                signal.alarm(0)
                # write to sdf
                Chem.SDWriter(save_path).write(merged_fragment)

    print('Success: ', success_cnt)
    print('Fail: ', fail_cnt)


def prepare_fragments(ground_truth_path, fragment_sdf_path, pdb_dir, cleaned_pdb_dir, aligned_pdb_dir, fragment_dir, merged_fragment_dir):
    """Create fragments dataset by merging known fragments together with Fragmenstein
    Args:
        ground_truth_path (str): Path to ground truth ligand-protein PDB file.
        fragment_sdf_path (str): Path to SDF file containing fragments.
        pdb_dir (str): Directory containing fragment-protein PDB files.
        aligned_pdb_dir (str): Directory to save aligned fragment-protein PDB files to.
        fragment_dir (str): Directory to save fragment PDB files to.
        merged_fragment_dir (str): Directory to save merged fragment PDB files to.
    """
    # Clean PDBs by relabeling ligands
    clean_pdb(pdb_dir, cleaned_pdb_dir)

    # Align PDB file to ground truth
    align_pdb_to_ground_truth(in_dir=cleaned_pdb_dir, out_dir=aligned_pdb_dir, ground_truth_path=ground_truth_path)
    
    batch_1_sdf_dir = os.path.dirname(os.path.dirname(fragment_sdf_path))
    # Extract coordinates from aligned PDB file, match atom ids to atoms from sdf file, and write to new sdf file
    extract_sdf(in_dir=aligned_pdb_dir, out_dir=fragment_dir, batch_1_pdb_dir=cleaned_pdb_dir, batch_1_sdf_dir=batch_1_sdf_dir, batch_2_sdf_path=fragment_sdf_path)

    # # Cluster fragments
    # clusters = cluster_ligands(args.ligand_dir)

    # Merge fragments
    merge_fragments(fragment_dir=fragment_dir, merged_fragment_dir=merged_fragment_dir, reference_pdb_path=ground_truth_path)