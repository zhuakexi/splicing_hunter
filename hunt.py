import sys
import pickle

from parser import pairs_parser
from block_search import block_search
import random

def clean_splicing_main(index_name, binsize, cell_name):
    # load directly from pickled bin_index
    cell = pairs_parser(cell_name)
    with open(index_name,"rb") as f:
        bin_index = pickle.load(f)
    return block_search(bin_index, BINSIZE, cell)
if __name__ == "__main__":
    BINSIZE = 10000

    index_name = sys.argv[1]
    BINSIZE=int(sys.argv[2])
    cell_name = sys.argv[3]
    
    result = clean_splicing_main(index_name, BINSIZE, cell_name)
    print(result)
