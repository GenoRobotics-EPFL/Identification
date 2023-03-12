import pandas as pd


genes = pd.DataFrame( columns= ["Family", "species", "sequence"])
genes.loc[len(genes)] = ['Homo', 'Homo Sapiens', 'ATGCGTAGCTAGCT']
print(genes)