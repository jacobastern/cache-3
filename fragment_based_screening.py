
import argparse

from src.prepare_fragments import prepare_fragments
from src.screen_enamine_real import screen_enamine_real


parser = argparse.ArgumentParser(description='Screen Enamine REAL database against a docked fragment database.')

parser.add_argument('--enamine_real_database_path', help='Path to Enamine REAL database.', type=str, default='./data/0_raw_data/Enamine_REAL_350-3_lead-like_cxsmiles.cxsmiles')
parser.add_argument('--ground_truth_path', help='Path to ground truth ligand-protein PDB file.', type=str, default='./data/1_bound_fragments/Mac1-DLS-EU0034_0A_bound.pdb')
parser.add_argument('--fragment_sdf_path', help='Path to SDF file containing fragments.', type=str, default='./data/0_raw_data/Mac1/Mac1_combined.sdf')
parser.add_argument('--pdb_dir', help='Directory containing fragment-protein PDB files.', type=str, default='./data/1_bound_fragments/')
parser.add_argument('--cleaned_pdb_dir', help='Directory to save cleaned fragment PDB files to.', type=str, default='./data/2_bound_fragments_cleaned/')
parser.add_argument('--aligned_pdb_dir', help='Directory to save aligned fragment-protein PDB files to.', type=str, default='./data/3_bound_fragments_aligned/')
parser.add_argument('--fragment_dir', help='Directory to save fragment PDB files to.', type=str, default='./data/4_fragment_sdfs/')
parser.add_argument('--merged_fragment_dir', help='Directory to save merged fragment PDB files to.', type=str, default='./data/5_merged_fragment_sdfs/')
parser.add_argument('--fragment_scaffolds_dir', help='Directory to save fragment scaffolds to.', type=str, default='./data/6_ligand_scores/')
parser.add_argument('--scores_path', help='Path to save results to.', type=str, default='./data/6_ligand_scores/scores.pkl')
parser.add_argument('--num_processes', help='Number of processes to use.', type=int, default=100)



if __name__ == '__main__':
    args = parser.parse_args()
    
    # prepare_fragments(ground_truth_path=args.ground_truth_path, fragment_sdf_path=args.fragment_sdf_path, pdb_dir=args.pdb_dir, cleaned_pdb_dir=args.cleaned_pdb_dir, aligned_pdb_dir=args.aligned_pdb_dir, fragment_dir=args.fragment_dir, merged_fragment_dir=args.merged_fragment_dir)
    screen_enamine_real(enamine_real_database_path=args.enamine_real_database_path, fragment_dir=args.fragment_dir, merged_fragment_dir=args.merged_fragment_dir, scores_path=args.scores_path, num_processes=args.num_processes)
