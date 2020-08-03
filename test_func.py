# test functions in block search


''' 
from block_search import filt_in_exon
from block_search import filt_relate_exon
'''
''' test filt_in_exon
exons = [(1,2),(3,4),(5,6),(4,7)]
locus = 5.5
print(filt_in_exon(locus, exons))
'''
''' test filt_relate_exon using manual data
exons = [(1,2),(3,4),(5,6),(4,7)]
bin = (1.5,3.5)
print(filt_relate_exon(bin,exons))
bin = (1,2)
print(filt_relate_exon(bin,exons))
'''
''' test set_list(deleted)
from block_search import set_list
a = []
set_list(a,4,["hello world"])
print(a)
set_list(a,3,200)
print(a)
set_list(a,5,["hi"])
print(a)
set_list(a,1000000,["hello world"])
print(len(a))
'''
''' test build_index using random generated data
import random
def random_exons(num, max_length, lower, upper)->"list of (start,end,'nonsense')":
    starts = random.sample(range(lower, upper-max_length), num)
    return [(start, random.randint(start, start+max_length), str(random.randint(0,num))) for start in starts]

from block_search import build_index
chromsomes = {"chr1":random_exons(20, 4, 1, 200),"chr2":random_exons(50, 4, 1, 200)}
bins = [ ("chr1",i,i+5) for i in range(0,150,5)  ] + [("chr2",i,i+5) for i in range(0,200,5)]
bin_index = build_index(bins, 5, chromsomes)
print(bin_index["chr1"])
print(bin_index["chr2"])
'''
