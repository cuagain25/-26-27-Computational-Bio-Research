## Specific Aims
1. **Aim 1:** Analyze the population frequency and genotype distribution of ADRB1 rs1801253 using public genomic data (gnomAD).
2. **Aim 2:** Construct and validate Wild-Type ($\text{Arg}_389$) and Mutant ($\text{Gly}_389$) $\beta_1$-adrenergic receptor ($\beta_1\text{AR}$) structural models.
3. **Aim 3:** Perform molecular docking of multiple $\beta$-blockers (e.g., metoprolol, bisoprolol, carvedilol) against WT and mutant receptors using AutoDock Vina.
4. **Aim 4:** Quantitatively analyze binding affinities, interaction profiles, and generate genotype-specific drug rankings validated against published pharmacogenomic literature.

---

## Computational Workflow
1. **Population Data:** Extract allele/genotype frequencies via gnomAD.
2. **Structural Modeling:** Retrieve template structures from RCSB PDB and generate variant models.
3. **Ligand Preparation:** Preprocess $\beta$-blocker structures using RDKit.
4. **Molecular Docking:** Run batch docking simulations using AutoDock Vina.
5. **Analysis:** Evaluate binding scores, hydrogen/hydrophobic interactions, and structural variations.

## GitHub Repository Structure

Computational Bio Research/
│
├── README.md                 # Overview and documentation for the project
├── config.txt                # Configuration file for AutoDock Vina Grid Box parameters
│
├── data/
│   ├── raw/                  # Raw input datasets downloaded from databases
│   │   ├── 2VT4.pdb
│   │   ├── atenolol.sdf
│   │   └── metoprolol.sdf
│   │
│   └── processed/            # Preprocessed files ready for docking runs
│       ├── 2VT4_receptor.pdbqt
│       └── atenolol.pdbqt
│
├── figures/                  # PyMOL renders, structural captures, and analysis plots
│   └── 2VT4_active_site.png
│
└── scripts/                  # Python scripts and shell utilities for processing/analysis
    └── prep_molecules.py

