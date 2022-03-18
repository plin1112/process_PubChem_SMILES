## process_PubChem_SMILES

python scripts used to process the SMILES files obtained from PubChem

## General Information

 - sort_CID_SMILES_byElements.py split the CID-SMILES by element groups, and by complex types, ion or neutral, etc.
 - sort_CID_SMILES_byHeavyAtoms.py split the CID-SMILES by number of heavy atoms
 - sort_stat_CID_SMILES_ECFP_r3.py process each splitted CID-SMILES, calculate fingerprints and count and save the statistics (recommand to run through splitted element groups and by number of heavy atoms
 - csv2pqz_CID_SMILES_ECFP_r3.py split the statistics of fingerprints based on radius used and save in dataframe
 - sort_pick_split_smiles_RDKit.py further splitting SMILES presenting the unique fingerprint by charge and multiplicity
 - gen_xtbxyz_RDKit_from_csv.py Generate xyz files ready for xTB calculation from csv formatted file
