import sys
import pickle
import random

from parser import bed_parser, bin_parser, con_parser
from block_search import build_index, block_search

BINSIZE = 10000
ref_name = sys.argv[1]
bin_name = sys.argv[2]
sample_rate = int(sys.argv[3])
try:
    index_file_name = sys.argv[4]
except IndexError:
    index_file_name = "ref/bin_10k_d{}_index".format(sample_rate)
chromsomes, ref_df = bed_parser(ref_name,"on")
bins = bin_parser(bin_name,"on")

#build 1/x index
## sample different num according to exon number
sample1 = {name:random.sample(chromsomes[name],len(chromsomes[name])//sample_rate) for name in chromsomes}
bin_index = build_index(bins, BINSIZE, sample1)
with open(index_file_name,"wb") as f:
    pickle.dump(bin_index,f)

