import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRotatableBonds 
import rdkit.Chem.Descriptors

MAX_mols_s = 6000
MAX_mols_d = 1000
MAX_mols_cs = 2000
MAX_mols_cd = 1000
MAX_confs = 30 

def generate_conformations(mol, max_confs=MAX_confs):
    rot_bond = CalcNumRotatableBonds(mol)
    confs = min(1 + 5*rot_bond, max_confs)
    AllChem.EmbedMultipleConfs(mol, numConfs=confs)
    return mol

def write_xtb_input_file(fragment, fragment_name, charge, uhf):
    number_of_atoms = fragment.GetNumAtoms()
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()]
    if charge != 0:
        if uhf == 1:
            fragment_name = 'xtbxyz_charged_doublet_' + fragment_name
        else:
            fragment_name = 'xtbxyz_charged_singlet_' + fragment_name
    else:
        if uhf == 1:
            fragment_name = 'xtbxyz_doublet_' + fragment_name
        else:
            fragment_name = 'xtbxyz_singlet_' + fragment_name

    for i,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+"_"+str(i)+".xyz"
        with open(file_name, "w") as file:
            file.write(str(number_of_atoms)+"\n")
            file.write("title\n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
                file.write(line)
            if charge !=0:
                if uhf == 1:
                    file.write("$set\n")
                    file.write("chrg "+str(charge)+"\n")
                    file.write("uhf "+str(1)+"\n")
                    file.write("$end")
                else:
                    file.write("$set\n")
                    file.write("chrg "+str(charge)+"\n")
                    file.write("$end")
            else:
                if uhf == 1:
                    file.write("$set\n")
                    file.write("uhf "+str(1)+"\n")
                    file.write("$end")
    return 

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Generate xyz files ready for xTB calculation from csv formatted file")
    parser.add_argument('-p', '--processed', nargs='*', default=None, help="previously processed csv files")
    args = parser.parse_args()

    if args.processed:
        exist_processed_cids = False
        for f in args.processed:
            processed_df = pd.read_csv(f, usecols =['count', 'CID', 'smiles', 'charge', 'multi'])
            new_processed_cids = set(processed_df['CID'].tolist())
            if exist_processed_cids:
                processed_cids = set.union(processed_cids, new_processed_cids)
            else:
                processed_cids = new_processed_cids

    tdf = pd.read_csv(args.input, usecols =['count', 'CID', 'smiles', 'charge', 'multi'])
    aggregation_functions = {'count': 'sum', 'CID': 'first', 'smiles': 'first', 'charge': 'first', 'multi': 'first'}
    tdf_new = tdf.groupby('smiles', as_index=False).aggregate(aggregation_functions)
    tdf_new.sort_values(by='count', ascending=False, inplace=True)

    n_mols = {'s':0, 'd':0, 'cs':0, 'cd':0}
    for idx in tdf_new.index:
        cid = tdf_new['CID'][idx]
        if cid in processed_cids:
            continue
        smiles_str = tdf_new['smiles'][idx]
        charge = tdf_new['charge'][idx]
        multi = tdf_new['multi'][idx]
        uhf = multi - 1
        if charge != 0:
            if uhf == 1:
                n_mols['cd'] += 1
                if n_mols['cd'] > MAX_mols_cd:
                    continue
            else:
                n_mols['cs'] += 1
                if n_mols['cs'] > MAX_mols_cs:
                    continue 
        else:
            if uhf == 1:
                n_mols['d'] += 1
                if n_mols['d'] > MAX_mols_d:
                    continue 
            else:
                n_mols['s'] += 1
                if n_mols['s'] > MAX_mols_s:
                    continue 

        mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
        mol = Chem.AddHs(mol)
        formula = CalcMolFormula(mol)
        mol = generate_conformations(mol)
        file_name = formula + '_CID' + str(cid)
        write_xtb_input_file(mol, file_name, charge, uhf)

    print(n_mols)
