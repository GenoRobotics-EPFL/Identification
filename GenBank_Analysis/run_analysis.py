from analysis_functions import *
import os.path as ospath
import sys
sys.path.insert(0,os.getcwd())
print(sys.path)

from utils import *

path = download_database("matk", 750, 751)
print("done download")
print(path)
gene_df = parse_data(path)
print("percentage of unverified sequences:", calculate_unverified_percentage(path))
print("mean ratio of nucleotides")
ratios = nucleotide_ratio(gene_df)
print("A: ", ratios.A.mean())
print("T: ", ratios["T"].mean())
print("C: ", ratios.C.mean())
print("G: ", ratios.G.mean())
show_length_graph(path)