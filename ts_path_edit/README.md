# Reaction generation development and comments.

Notes:
omol_mc_reactions_species.csv 
omol_mc_reactions.csv
reaction_complexes/
species_geometries/

Curated by Sam Blau (smblau@lbl.gov).

## RxN workflow generation:

1. Curate metals dataset to match oxidation state. Use mendeleev to detect spin. -> Done (1_curate_m_ligs_dataset_for_sampling/curate_m_ligds.ipynb)
2. Curate ligands dataset to match denticity/charge states. -> Done (1_curate_m_ligs_dataset_for_sampling/curate_m_ligds.ipynb)

3. Inspect all reactions and add meta-data column with index of metal atom (s) and oxidation states. -> Done (2_inspect_rxns/inspect_rxns.ipynb)
4. Write routine to detect/create possible ligand swaps given reaction complex. (init vs final) don't swap any ligands involved in bond breaking/forming. Make sure ligands swappable are only up to dent 2. (Check that charge detection is working correctly for all ligands in the complex.) -> Done (2_inspect_rxns/inspect_rxns.ipynb)

5. Create routine for sampling reproducibly from metal/ligand swaps. (Sample from 1-n ligands to swap, sample from metal swap) -> Done (3_swap_ligands/1_sampler_develop_testing.ipynb)
6. Create routine for implementing a selected metal/ligand swap. Swap metal(s). Update complex spin/charge. Delete ligand atoms from selected ligands to swap - ideally, save metal-coordination locations. Add ligands + create bonds - in 3D. See where UFF takes us. Maybe do constrained XTB relaxation. -> Done (3_swap_ligands/1_sampler_develop_testing.ipynb) -> Do UFF rather than XTB for speed.
7. Check indices of the newly generated reactions. They shouldn't be messed with. Check the 3D structures generated that not overlapping. -> Done (3_swap_ligands/1_sampler_develop_testing.ipynb) -> Do UFF rather than XTB for speed.

8. Down-select from possible reactions so that ~250K reactions generated as uniformly as possible across different templates. Flag the "test" set (e.g. no ligand swaps possible rxns). Include all rxns from templates with from those less than 250K/125. Take (250-all enumerated)/n_remaining samples from remaining rows if possible. -> Done (3_swap_ligands/2_sample_downselect.ipynb)
9. Run reaction generation script over all selected reactions so that as many as possible generated. -> Done (4_production/swap_production.py)


### To run initial/final structure generation and create combined sdf files.

# With initial/final structures and assigned bonds.

All you should need are the files in 4_production.

Adapt to your HPC.