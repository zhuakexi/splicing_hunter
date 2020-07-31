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
    candidate = []
    #print(start)
    '''
    if start == 0 or start==len(exons):
        print("two_way_search: boundary")
        return []
    '''
    if start == 0:
        this = start
        while (exons[this][0] <= locus) and (this<=len(exons)):
            candidate.append((exons[this]))
            this += 1
    elif start == len(exons):
        this = start - 1
        while (exons[this][1] >= locus) and (this>=0):
            candidate.append(exons[this])
            this -= 1
    else:
        candidate = []
        this = start - 1
        while (exons[this][1] >= locus) and (this>=0):
            candidate.append(exons[this])
            this -= 1
        this = start
        while (exons[this][0] <= locus) and (this<len(exons)):
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
    begin = time.time()
    left_hit_contacts = [] #test mid-search method
    result = []
    for contact in cell:
        left_hit_exons, right_hit_exons = [],[]
        chr, leg1, leg2 = contact[CHR], contact[LEG1], contact[LEG2]
        left_hit_exons = two_way_search(chromsomes[chr], keys_of_all[chr], leg1)
        right_hit_exons = two_way_search(chromsomes[chr], keys_of_all[chr], leg2)
        if left_hit_exons != []: # test mid-search method
            left_hit_contacts.append(contact)
        if left_hit_exons != [] and right_hit_exons != [] :
            #print(contact,out)
            left_hit_genes = set([exon[2] for exon in left_hit_exons])
            right_hit_genes = set([exon[2] for exon in right_hit_exons])
            #print(out_names)
            if left_hit_genes.intersection(right_hit_genes) != set():
                #print(contact)
                result.append(contact)
    #sys.stderr.write("Total contacts:" + str( len(ref_df)) + "\n")
    #sys.stderr.write("Contacts hit: " + str( len(result) )+ " ratio: " + str( len(result)/len(ref_df)) +"\n")
    print("Mid_search left hit: " , str(len(left_hit_contacts))) # test mid-search method
    sys.stderr.write("Mid_search contacts hit: " + str( len(result) ) + "\n")
    sys.stderr.write("Mid_search time: " + str( time.time()-begin ) + "\n")
    return result
