import sys
import pickle
import random
import time
from concurrent import futures

from parser import bed_parser, bin_parser, con_parser
from block_search import build_chromsome_index, block_search

BINSIZE = 10000
ref_name = sys.argv[1]
bin_name = sys.argv[2]
index_file_name = sys.argv[3]
'''
sample_rate = int(sys.argv[3])
try:
    index_file_name = sys.argv[4]
except IndexError:
    index_file_name = "ref/bin_10k_d{}_index".format(sample_rate)
'''
chromsomes, ref_df = bed_parser(ref_name,"on")
bins = bin_parser(bin_name,"on")

'''
#build 1/x index
## sample different num according to exon number
sample1 = {name:random.sample(chromsomes[name],len(chromsomes[name])//sample_rate) for name in chromsomes}
t0 = time.time()
## generate input list for map 
input = [(chr_name, bins[chr_name], BINSIZE, sample1[chr_name]) for chr_name in bins]
## multiplex mapping
with futures.ProcessPoolExecutor() as executor:
    res = executor.map(build_chromsome_index, input)
## regenerate the dictionary from tuple-list
bin_index = {chr_name:p_bin_index for chr_name, p_bin_index in res}
sys.stderr.write("multi_build_index: finish in %s \n" % str(time.time()-t0))
## write the index to file
with open(index_file_name,"wb") as f:
    pickle.dump(bin_index,f)
'''

#build full index
t0 = time.time()
input = [(chr_name, bins[chr_name], BINSIZE, chromsomes[chr_name]) for chr_name in bins]
with futures.ProcessPoolExecutor() as executor:
    res = executor.map(build_chromsome_index, input)
bin_index = {chr_name:p_bin_index for chr_name, p_bin_index in res}
sys.stderr.write("multi_build_index: finish in %s \n" % str(time.time()-t0))
with open(index_file_name,"wb") as f:
    pickle.dump(bin_index,f)
