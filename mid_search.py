import time
import sys
import pandas as pd
import numpy as np
import bisect
import random

CHR = 0
LEG1 = 1
LEG2 = 2

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



'''
for i in range(0, 20):
    print(chromsomes["chr1"][i], keys_of_all["chr1"][i])
print("length of chr1 keys:", len(keys_of_all["chr1"]))
'''
'''
sample = random.sample(cell, 100)
for contact in sample:
    position = bisect.bisect( keys_of_all[contact[CHR]], float(contact[LEG1]) ) 
    try:
        print(contact, contact[LEG1], keys_of_all[contact[CHR]][position],chromsomes[contact[CHR]][position-1], chromsomes[contact[CHR]][position], chromsomes[contact[CHR]][position+1])
    except IndexError:
        print("edge",contact)
'''

def mid_search(cell:list, keys_of_all:dict, chromsomes:dict) -> list:
    ##in: cell, keys_of_all, chromsomes
    ##out: *result
    left_hit = []
    right_hit = []
    result = []
    for contact in cell:
        chr = contact[CHR]
        leg1 = contact[LEG1]
        leg2 = contact[LEG2]
        left_hit = two_way_search(chromsomes[chr], keys_of_all[chr], leg1)
        right_hit = two_way_search(chromsomes[chr], keys_of_all[chr], leg2)
        if left_hit != []:
            print(contact, left_hit)
        if left_hit != [] and right_hit != [] :
            #print(contact,out)
            left_names = set([exon[2] for exon in left_hit])
            right_names = set([exon[2] for exon in right_hit])
            #print(out_names)
            if left_names.intersection(right_names) != set():
                #print(contact)
                result.append(contact)
    return result
