# Fragment-based lead optimization

This code takes protein structures with bound ligand fragments (<250 daltons) and a ligand library (e.g. Enamine REAL) and filters the database to select ligands most similar to the fragments.

This is the workflow:
Step 0: Get known binding fragments - 205 binders from paper 1,  5 binders from paper 2 (5S4F, 5S4G, 5S4H, 5S4I, 5S4J)
Step 1: Increase size of binding fragments dataset by merging fragments with Fragmentstein
Step 2: Extract Bemis-Murcko (BM) scaffolds from fragments
Step 3: Score all Enamine ligands. Points for 1) including pyramidine (common substructure of binders), 2) being an anion (the known binder, ADPr, is an anion), 3) including a scaffold of a binding fragment, (measured using RDKit maximum common substructure on BM Scaffolds, with more points for including the scaffold of a larger binding fragment), 4) summed fragment tanimoto similarity score
Step 4: Screen the top L scoring ligands with xtb/US-HRE

```
docker build -t cache-3 -f dependencies/Dockerfile dependencies
docker run -dit --gpus all --rm --name cache -v $(pwd):/code cache-3
docker exec -it docking /bin/bash
cd code
python3 fragment_based_screening.py
```