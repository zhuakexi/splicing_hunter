import sys
import pickle
from parser import bed_parser, bin_parser, con_parser
from block_search import build_index, block_search
import random
ref_name = sys.argv[1]
bin_name = sys.argv[2]
cell_name = sys.argv[3]
chromsomes, ref_df = bed_parser(ref_name,"on")
cell = con_parser(cell_name)
bins = bin_parser(bin_name,"on")
bin_index = build_index(bins, 10000, chromsomes)
with open("ref/bin_10k_index","wb") as f:
    pickle.dump(bin_index,f)

