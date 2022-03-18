import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Descriptors
import csv

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Existing dataframe of fingerprint statistics in parquet format")
    parser.add_argument('-o', '--output', required=True, help="prefix of output files in  csv format")
    args = parser.parse_args()

    tdf = pd.read_parquet(args.input, columns=['count', 'CID', 'smiles'])
    aggregation_functions = {'count': 'sum', 'CID': 'first', 'smiles': 'first'}
    tdf_new = tdf.groupby('smiles', as_index=False).aggregate(aggregation_functions)
    tdf_new.sort_values(by='count', ascending=False, inplace=True)

    header = ['count', 'CID', 'smiles', 'charge', 'multi']
    cf_s = open(args.output+'_s.csv', 'a') # neutral singlet 
    writer_s = csv.DictWriter(cf_s, fieldnames=header)
    writer_s.writeheader()
    cf_cs = open(args.output+'_cs.csv', 'a') # charged singlet
    writer_cs = csv.DictWriter(cf_cs, fieldnames=header)
    writer_cs.writeheader()
    cf_d = open(args.output+'_d.csv', 'a') # neutral doublet
    writer_d = csv.DictWriter(cf_d, fieldnames=header)
    writer_d.writeheader()
    cf_cd = open(args.output+'_cd.csv', 'a') # charged doublet
    writer_cd = csv.DictWriter(cf_cd, fieldnames=header)
    writer_cd.writeheader()

    # the multiplicity is determined by even or odd number of valence electrons, higher spin state not considered.
    n_mols = {'s':0, 'd':0, 'cs':0, 'cd':0}
    for idx in tdf_new.index:
        row = tdf_new.loc[[idx]]
        cid = tdf_new['CID'][idx]
        smiles_str = tdf_new['smiles'][idx]
        mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
        mol = Chem.AddHs(mol)
        charge = Chem.GetFormalCharge(mol)
        # charge = 0 # can be set if known
        valence = Chem.Descriptors.NumValenceElectrons(mol)
        uhf = valence % 2
        new_dict = {}
        if charge != 0:
            if uhf == 1:
                n_mols['cd'] += 1
                new_dict = row.to_dict('records')[0]
                new_dict['charge'] = charge
                new_dict['multi'] = 2
                writer_cd.writerow(new_dict)
            else:
                n_mols['cs'] += 1
                new_dict = row.to_dict('records')[0]
                new_dict['charge'] = charge
                new_dict['multi'] = 1
                writer_cs.writerow(new_dict)
        else:
            if uhf == 1:
                n_mols['d'] += 1
                new_dict = row.to_dict('records')[0]
                new_dict['charge'] = 0
                new_dict['multi'] = 2
                writer_d.writerow(new_dict)
            else:
                n_mols['s'] += 1
                new_dict = row.to_dict('records')[0]
                new_dict['charge'] = 0
                new_dict['multi'] = 1
                writer_s.writerow(new_dict)

    print(n_mols)

    cf_s.close()
    cf_cs.close()
    cf_d.close()
    cf_cd.close()

