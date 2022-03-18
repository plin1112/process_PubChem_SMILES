import time
startTime = time.time()

import pandas as pd
import argparse
import glob
import os

parser = argparse.ArgumentParser(description="Processing analyzed fingerprints in csv file and save to dataframe")
parser.add_argument("-i", "--input", required=True, help="Existing statistics of fingerprints in csv format")
args = parser.parse_args()

fname = args.input
tdf = pd.read_csv(fname)

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))

tdf0 = tdf[tdf['radius'] == 0]
tdf1 = tdf[tdf['radius'] == 1]
tdf2 = tdf[tdf['radius'] == 2]
tdf3 = tdf[tdf['radius'] == 3]

tdf0.sort_values(by='count', ascending=False, inplace = True)
tdf1.sort_values(by='count', ascending=False, inplace = True)
tdf2.sort_values(by='count', ascending=False, inplace = True)
tdf3.sort_values(by='count', ascending=False, inplace = True)

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))

basename = os.path.basename(fname)

pqfile0 = basename + "_fp0.pqz"
pqfile1 = basename + "_fp1.pqz"
pqfile2 = basename + "_fp2.pqz"
pqfile3 = basename + "_fp3.pqz"

print(f"=> saving as parquet format at {pqfile0}")
tdf0.to_parquet(pqfile0, compression='gzip')
result = pd.read_parquet(pqfile0)
print(result.info())
print(result.head())

print(f"=> saving as parquet format at {pqfile1}")
tdf1.to_parquet(pqfile1, compression='gzip')
result = pd.read_parquet(pqfile1)
print(result.info())
print(result.head())

print(f"=> saving as parquet format at {pqfile2}")
tdf2.to_parquet(pqfile2, compression='gzip')
result = pd.read_parquet(pqfile2)
print(result.info())
print(result.head())

print(f"=> saving as parquet format at {pqfile3}")
tdf3.to_parquet(pqfile3, compression='gzip')
result = pd.read_parquet(pqfile3)
print(result.info())
print(result.head())

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))


