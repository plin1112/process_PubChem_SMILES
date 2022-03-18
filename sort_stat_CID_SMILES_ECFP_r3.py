import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import csv

BUFFER_SIZE_LINES = 100000  # Maximum number of lines to buffer in memory
radius_size = 3

import argparse
import glob

parser = argparse.ArgumentParser(description="Processing SMILES for unique fingerprints")
parser.add_argument("-i", "--input", required=True, help="Input CID-SMILES filename")
parser.add_argument("-d", "--dicts", default=None, help="Existing dictionary with unique fingerprints in csv format")
args = parser.parse_args()

fname = args.input

if args.dicts:
    tdf = pd.read_csv(args.dicts)
else:
    tdf = pd.DataFrame(columns = ['fp', 'radius', 'count', 'CID', 'smiles'])

print(f'Start processing {fname} ...')

r = open(fname, "r")
f_uni_r0 = fname + '-unique-r0'
f_uni_r1 = fname + '-unique-r1'
f_uni_r2 = fname + '-unique-r2'
f_uni_r3 = fname + '-unique-r3'
f_ext = fname + '-extend'
f_fail = fname + '-failed'
f_dict = fname + '-dict'
w_uni_r0 = open(f_uni_r0, "w")
w_uni_r1 = open(f_uni_r1, "w")
w_uni_r2 = open(f_uni_r2, "w")
w_uni_r3 = open(f_uni_r3, "w")
w_ext = open(f_ext, "w")
w_fail = open(f_fail, "w")

buf_uni_r0 = ""
buf_uni_r1 = ""
buf_uni_r2 = ""
buf_uni_r3 = ""
buf_ext = ""
buf_fail = ""
fp_dict = {}

bufLines = 0
for lineIn in r:
    blocks = lineIn.split()
    cid = int(blocks[0])
    smiles_str = blocks[1]
    unique_fp_r0 = False
    unique_fp_r1 = False
    unique_fp_r2 = False
    unique_fp_r3 = False
    if '*' in smiles_str:
        buf_fail += lineIn
    else:
        mol = Chem.MolFromSmiles(smiles_str, sanitize=True)
        try:
            fp_info={}
            fp = AllChem.GetMorganFingerprint(mol, radius_size, bitInfo=fp_info)
            for k, v in fp_info.items():
                idx = tdf.index[tdf['fp'] == k]
                if idx.size > 0:
                    for at, ra in v:
                        tdf.loc[idx]['count'] += 1 
                elif k in fp_dict.keys():
                    for at, ra in v:
                        fp_dict[k]['count'] += 1
                else:
                    for at, ra in v:
                        fp_dict[k] = {'CID':cid, 'count':1, 'radius':ra, 'smiles':smiles_str}
                    if ra == 0:
                        unique_fp_r0 = True
                    elif ra == 1:
                        unique_fp_r1 = True
                    elif ra == 2:
                        unique_fp_r2 = True
                    else:
                        unique_fp_r3 = True
                    
            if unique_fp_r0:
                buf_uni_r0 += lineIn 
            elif unique_fp_r1:
                buf_uni_r1 += lineIn
            elif unique_fp_r2:
                buf_uni_r2 += lineIn
            elif unique_fp_r3:
                buf_uni_r3 += lineIn
            else:
                buf_ext += lineIn
        except:
            buf_fail += lineIn 

    bufLines += 1
    if bufLines >= BUFFER_SIZE_LINES:
        # Flush buffer to disk
        if buf_uni_r0:
            w_uni_r0.write(buf_uni_r0)
        if buf_uni_r1:
            w_uni_r1.write(buf_uni_r1)
        if buf_uni_r2:
            w_uni_r2.write(buf_uni_r2)
        if buf_uni_r3:
            w_uni_r3.write(buf_uni_r3)
        if buf_ext:
            w_ext.write(buf_ext)
        if buf_fail:
            w_fail.write(buf_fail)

        buf_uni_r0 = ""
        buf_uni_r1 = ""
        buf_uni_r2 = ""
        buf_uni_r3 = ""
        buf_ext = ""
        buf_fail = ""

        bufLines = 0

# Flush remaining buffer to disk
if buf_uni_r0:
    w_uni_r0.write(buf_uni_r0)
if buf_uni_r1:
    w_uni_r1.write(buf_uni_r1)
if buf_uni_r2:
    w_uni_r2.write(buf_uni_r2)
if buf_uni_r3:
    w_uni_r3.write(buf_uni_r3)
if buf_ext:
    w_ext.write(buf_ext)
if buf_fail:
    w_fail.write(buf_fail)

w_uni_r0.close()
w_uni_r1.close()
w_uni_r2.close()
w_uni_r3.close()
w_ext.close()
w_fail.close()

# dump fingerprint dictionary in pickle
with open(f_dict, 'wb') as f:
    # Pickle the 'fp_dict' dictionary using the highest protocol available.
    pickle.dump(fp_dict, f, pickle.HIGHEST_PROTOCOL)

# overwrite previous csv file
tdf.to_csv(args.dicts, index=False)

cf = open(args.dicts, 'a')
header = ['fp', 'radius', 'count', 'CID', 'smiles']
writer = csv.DictWriter(cf, fieldnames=header)

for k, v in fp_dict.items():
    cid = v['CID']
    count = v['count']
    radius = v['radius']
    smiles_str = v['smiles']
    new_row = {'fp':k, 'radius':radius, 'count':count, 'CID':cid, 'smiles':smiles_str}
    writer.writerow(new_row)

cf.close()

r.close()
