import pandas as pd
import numpy as np
import math

# This code takes SANDPUMA predictions output that is merged into a single TSV file
# and creates a proper table containing ID of the input genome, locus tag, module number and list of predictions
# of each method with full names of substrates

full_names_dict = {}

# Read a file containing short names of substrates (provided by SANDPUMA) and their full names (annotated manually)
with open('specificities_dict.txt', 'r') as file:
    for line in file:
        short_name, full_name = line.strip().split('\t')
        full_names_dict[short_name] = full_name

# parse SANDPUMA output
df = pd.read_csv('sandpuma_basic_output.tsv', sep='\t', header=None, names=['genomic_id', 'locus_tag', 'method', 'prediction'])
df[['locus_tag', 'module']] = df['locus_tag'].str.rsplit('_', n=1, expand=True)

df['prediction'] = df['prediction'].astype(str).str.replace('|', ',')
df['prediction'] = np.where(df['prediction'].isnull() | (df['prediction'] == 'None'), np.nan, df['prediction'])
df['prediction'] = df['prediction'].apply(lambda x: [i.strip() for i in x.split(',')] if pd.notnull(x) else [])

df['prediction'] = df['prediction'].apply(lambda x: [full_names_dict.get(name, name) for name in x])

df['module'] = df['module'].str.replace('A', '').astype(int)

valid_methods = ['ASM', 'pHMM', 'SANDPUMA', 'SVM']

df = df[df['method'].astype(str).isin(valid_methods)]

new_df = pd.pivot_table(df, index=['genomic_id', 'locus_tag', 'module'], columns='method', values='prediction', aggfunc=lambda x: x)
new_df = new_df.reset_index()

new_df.to_csv('sandpuma_predictions.tsv', sep='\t', index=False)
