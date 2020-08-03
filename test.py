import sys
import pickle
ref_name = sys.argv[1]
bin_name = sys.argv[2]
cell_name = sys.argv[3]
from parser import bed_parser, bin_parser, con_parser
from block_search import build_index, block_search
import random
chromsomes, ref_df = bed_parser(ref_name,"on")
cell = con_parser(cell_name)
bins = bin_parser(bin_name,"on")
''' build index using 1/x fold sampling for test
sample1 = {name:random.sample(chromsomes[name],len(chromsomes[name])//10) for name in chromsomes}
bin_index = build_index(bins, 10000, sample1)
with open("ref/bin_10k_d10_index","wb") as f:
    pickle.dump(bin_index,f)
'''
# load directly from pickled bin_index
with open("ref/bin_10k_d10_index","rb") as f:
    bin_index = pickle.load(f)

result = block_search(bin_index, 10000, cell)
