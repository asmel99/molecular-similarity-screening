import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def run_similarity_screening(dataset_path, reference_smiles, threshold=0.5):
    print(f"--- Loading Dataset from {dataset_path} ---")
    df = pd.read_csv(dataset_path)
    
    # 1. Parse the known reference molecule (e.g., a known active compound)
    ref_mol = Chem.MolFromSmiles(reference_smiles)
    if ref_mol is None:
        print("❌ Error: Invalid reference SMILES string provided.")
        return
    
    # Generate Morgan Fingerprint for the reference (Radius=2, Bits=2048)
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
    
    print(f"Comparing dataset against reference molecule (Threshold >= {threshold})...")
    
    scores = []
    for smiles in df['Canonical_SMILES']:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Generate Morgan Fingerprint for dataset molecule
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            # Calculate Tanimoto Coefficient
            score = DataStructs.TanimotoSimilarity(ref_fp, fp)
            scores.append(round(score, 3))
        else:
            scores.append(0.0)
            
    # 2. Append results to DataFrame
    df['Tanimoto_Score'] = scores
    
    # Filter by threshold and sort by highest similarity
    filtered_df = df[df['Tanimoto_Score'] >= threshold]
    sorted_df = filtered_df.sort_values(by='Tanimoto_Score', ascending=False)
    
    # 3. Output results
    print(f"\n✅ Screening Complete! Found {len(sorted_df)} matches above threshold.")
    if not sorted_df.empty:
        print("\nTop Matches:")
        print(sorted_df[['Canonical_SMILES', 'Tanimoto_Score']].head(10))
        
        # Save results to a new CSV file
        output_file = "screening_hits.csv"
        sorted_df.to_csv(output_file, index=False)
        print(f"\n📁 Saved filtered hits to '{output_file}'")
    else:
        print("No molecules met the similarity threshold. Try lowering it slightly.")

if __name__ == "__main__":
    # Example test: Replace this with your actual reference SMILES 
    # (e.g., a known antifungal/antibacterial scaffold)
    sample_reference_smiles = "Clc1ccc(C(c2ccccc2)(c3ccccc3)N4C=CN=C4)cc1" # Example structure
    
    # Path to your cleaned dataset file
    dataset_csv = "clean_test.csv"
    
    run_similarity_screening(dataset_csv, sample_reference_smiles, threshold=0.4)
