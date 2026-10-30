# Molecular Similarity Screening Tool

A lightweight cheminformatics tool built with **RDKit** to perform ligand-based virtual similarity screening using **Morgan Fingerprints** and **Tanimoto Coefficients**. 

## Purpose
Used as a preliminary filtering step in drug and agrochemical discovery pipelines to evaluate whether known active reference molecules already match or resemble entries in custom chemical datasets before deploying *de novo* generative models.

## Usage
1. Place your clean dataset CSV file in the directory.
2. Update the `sample_reference_smiles` variable in `screen.py` with your target active compound.
3. Run the script:
   ```bash
   python screen.py
