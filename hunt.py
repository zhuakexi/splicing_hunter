import sys
import pickle

from parser import con_parser
from block_search import block_search
import random

BINSIZE = 10000

index_name = sys.argv[1]
cell_name = sys.argv[2]
cell = con_parser(cell_name)

# load directly from pickled bin_index
with open(index_name,"rb") as f:
    bin_index = pickle.load(f)
result = block_search(bin_index, BINSIZE, cell)
print("\n".join([str(i) for i in result]))
sys.stderr.write("block_search total questionable contacts: " + str(len(result)) + "\n")
