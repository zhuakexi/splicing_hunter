from mid_search import mid_search
from standard import filter_search
from parser import con_parser
from parser import bed_parser
import sys
ref_name = sys.argv[1]
cell_name = sys.argv[2]
#sample_count = sys.argv[3]
try:
    out_name = sys.argv[4]
except IndexError:
    out_name = cell_name.split(".")[0] + ".hit"
# --------parsing files--------
## in: ref_name, cell_name, out_name
## out: *chromsomes, *ref_df, *cell
chromsomes, ref_df = bed_parser(ref_name)
cell = con_parser(cell_name, 1000)
# ------------sort exons according to mid-point------------
## in: chromsomes
## out: keys_of_all, *chromsomes
keys_of_all = {}
for i in chromsomes:
    chromsomes[i].sort(key=lambda exon : (int(exon[0])+int(exon[1])) /2)
    keys_of_all[i] = [ (int(exon[0]) + int(exon[1]))/2 for exon in chromsomes[i]]
# ------------do two search------------

result_mid = mid_search(cell, keys_of_all, chromsomes)
result_filt = filter_search(cell, chromsomes, ref_df)

# ------------test filter and two_way_search-----------
'''
from mid_sarch import two_way_search
from standard import filt_in_exon

for contact in cell:
    chr, leg1, leg2 = contact[0], contact[1], contact[2]
'''
    