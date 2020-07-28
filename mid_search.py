import time
import sys
import pandas as pd
import numpy as np
import bisect
import random
from parser import bed_parser

CHR = 0
LEG1 = 1
LEG2 = 2
# .contact parser
def con_parser(cellname):
    #list of entries: chromsome, leg1, leg2
    time_begin = time.time()
    con = pd.read_table(cell_name,header=None,sep='[,\t]')
    cell = list(zip(con[0],con[1],con[4]))
    #check parsing time
    #print("parsing cell:", time.time() - time_begin)
    return cell
def in_exon(locus:int):
    def working_func(exon:tuple):
        return exon[0] <= locus <= exon[1]
    return working_func 
def filt_in_exon(locus:int, exons:list) -> list:
    return list(filter(in_exon(locus), exons))
def two_way_search(exons, keys, locus):
    start = bisect.bisect( keys, float(locus) )
    #print(start)
    if start == 0 or start==len(exons):
        return []
    else:
        candidate = []
        this = start - 1
        while (exons[this][1] >= locus) and (this>=0):
            candidate.append(exons[this])
            this -= 1
        this = start
        while (exons[this][0] <= locus) and (this<=len(exons)):
            candidate.append((exons[this]))
            this += 1
        #print(candidate)
        return filt_in_exon(locus, candidate)
# --------parsing files--------
bed_name = sys.argv[1]
cell_name = sys.argv[2]
out_name = sys.argv[3] 
try:
    out_name = sys.argv[3]
except IndexError:
    out_name = cell_name.replace(".contacts.pairs.gz", ".rs_can")
time_begin = time.time()
chromsomes, ref_df = bed_parser(bed_name)
cell = con_parser(cell_name)
#check parsing time
print("parsing :", time.time() - time_begin)

# ------------sort exons according to mid-point------------
## in: chromsomes
## out: keys_of_all, *chromsomes
begin = time.time()
keys_of_all = {}
for i in chromsomes:
    chromsomes[i].sort(key=lambda exon : (int(exon[0])+int(exon[1])) /2)
    keys_of_all[i] = [ (int(exon[0]) + int(exon[1]))/2 for exon in chromsomes[i]]
'''
for i in range(0, 20):
    print(chromsomes["chr1"][i], keys_of_all["chr1"][i])
print("length of chr1 keys:", len(keys_of_all["chr1"]))
'''
print("sort & keying time:", time.time()-begin)

'''
sample = random.sample(cell, 100)
for contact in sample:
    position = bisect.bisect( keys_of_all[contact[CHR]], float(contact[LEG1]) ) 
    try:
        print(contact, contact[LEG1], keys_of_all[contact[CHR]][position],chromsomes[contact[CHR]][position-1], chromsomes[contact[CHR]][position], chromsomes[contact[CHR]][position+1])
    except IndexError:
        print("edge",contact)
'''

# ------------search leg1------------
##in: cell, keys_of_all, chromsomes
##out:
left_hit = []
right_hit = []
result = []
for contact in cell:
    chr = contact[CHR]
    leg1 = contact[LEG1]
    leg2 = contact[LEG2]
    left_hit = two_way_search(chromsomes[chr], keys_of_all[chr], leg1)
    right_hit = two_way_search(chromsomes[chr], keys_of_all[chr], leg2)
    if left_hit != [] and right_hit != [] :
        #print(contact,out)
        left_names = set([exon[2] for exon in left_hit])
        right_names = set([exon[2] for exon in right_hit])
        #print(out_names)
        if left_names.intersection(right_names) != set():
            #print(contact)
            result.append(contact)
sys.stderr.write("Total contacts:" + str( len(ref_df)) + "\n")
sys.stderr.write("Contacts hit: " + str( len(result) )+ " ratio: " + str( len(result)/len(ref_df)) +"\n")
# ------------write output file------------
title = cell_name + " in " + bed_name
content = map(str, result)
content = [i.strip("()") for i in content]
content = "\n".join(content)
with open(out_name,"w") as f:
    f.write(title+"\n"+content)
