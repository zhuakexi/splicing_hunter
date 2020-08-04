import sys
import pickle
import random
from concurrent import futures

from parser import bed_parser, bin_parser, con_parser
from block_search import build_chromsome_index, block_search
ref_name = sys.argv[1]
bin_name = sys.argv[2]
cell_name = sys.argv[3]
chromsomes, ref_df = bed_parser(ref_name,"on")
cell = con_parser(cell_name)
bins, bin_df = bin_parser(bin_name,"on")

#build 1/x index
## sample different num according to exon number
sample1 = {name:random.sample(chromsomes[name],len(chromsomes[name])//1000) for name in chromsomes}
##get grouped bin by chr_name, get separate reference exon list by chromsome
divide_bins = {key:value for key, value in bin_df.groupby(0)}
binsize = 10000
input = [(chr_name,divide_bins[chr_name], binsize, sample1[chr_name]) for chr_name in divide_bins]
print(len(input))
print(len(input[0]))
print(input[0][0])
'''
with futures.ProcessPoolExecutor() as executor:
    executor.map(build_chromsome_index, input)
bin_index = build_index(bins, 10000, sample1)

with open("ref/bin_10k_d2_index","wb") as f:
    pickle.dump(bin_index,f)
'''
'''
#build full index
bin_index = build_index(bins, 10000, chromsomes)
with open("ref/bin_10k_index","wb") as f:
    pickle.dump(bin_index,f)
'''
