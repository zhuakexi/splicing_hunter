import time
import sys
import pandas as pd
import numpy as np
BINSIZE=10000
def filt_in_exon(locus, exons):
    '''
    filter out exons envolope the locus
    '''
    return [exon for exon in exons if exon[0] <= locus <= exon[1]]
def filt_relate_exon(bin, exons):
    '''
    filter out exons relate to this bin
    note: for bed file, bp start = bin[0] + 1, bp end = bin[1]
    '''
    return [exon for exon in exons if not(bin[2]+1<exon[0] or bin[1]>exon[1]) ]
def set_list_value(target, index, value):
    '''
    set element value by index
    if index out of range, fill in with []
    '''
    try:
        target[index] = value
    except IndexError:
        target.extend( ([] for i in range(index-len(target)+1)) )
        target[-1] = value
def build_index(bins:"dict of bins", bin_size:int, chromsomes:"dict; list of exons by chromsome") -> "list of exon list by bin by chromsome in dict":
    '''
    build binning index from bins and reference
    '''
    begin_time = time.time()
    '''
    if bin_size != (bins[0][2] - bins[0][1]):
        sys.stderr.write("buil_index: Warning. Wrong bin_size.\n")
    '''
    bin_index = {}
    '''
    chromsome_names = set()
    for bin in bins:
        chromsome_names.add(bin[0])
    '''
    #print("bin chromsome numbers: ", len(chromsome_names))
    '''
    for name in chromsome_names:
        bin_index[name] = [[] for bin in bins if bin[0] == name]
    '''
    #print(bin_index)
    for chr_name, p_bins in bins.items():
        bin_index[chr_name] = [[] for bin in p_bins] #initiate bin_index with list of blank list
        for bin in p_bins:
            #print(bin)
            exons = chromsomes[chr_name]
            related = filt_relate_exon(bin, exons)
            the_index = bin_to_index(bin, bin_size)
            bin_index[chr_name][the_index] = related
    sys.stderr.write("build_index building time: " + str(time.time()-begin_time) + "\n")
    return bin_index
def build_chromsome_index(chr_name, p_bins:"list of bins of one chromsome", bin_size:int, chromsome:"list of exons of one chromsome") ->"exon list":
    t0 = time.time()
    p_bin_index = [[] for bin in p_bins] #initiate list of blank list same length as p_bins
    if bin_size != (p_bins[0][2] - p_bins[0][1]):
        sys.stderr.write("buil_index: Warning. Wrong bin_size.\n")
    for bin in p_bins:
        related = filt_relate_exon(p_bins, chromsome)
        the_index = bin_to_index(bin, bin_size)
        p_bin_index[the_index] = related
    sys.stderr.write("build_chromsome_index: build %s in %s \n" %( chr_name, str(time.time()-begin_time)) )
    return (chr_name, p_bin_index)
def key_to_index(locus, bin_size):
    return locus//bin_size
def bin_to_index(bin, bin_size):
    '''
    bin[1] can always been divided exactly by BINSIZE 
    '''
    return bin[1]//bin_size
def block_search(bin_index, binsize, cell):
    begin_time = time.time()
    result = []
    for contact in cell:
        #search left
        chr_name, left_index, right_index = contact[0], key_to_index(contact[1], binsize), key_to_index(contact[2], binsize)
        left_hit_exons = filt_in_exon(contact[1], bin_index[chr_name][left_index])
        if left_hit_exons != []:
            right_hit_exons = filt_in_exon(contact[1], bin_index[chr_name][right_index])
            left_hit_genes = set([exon[2] for exon in left_hit_exons])
            right_hit_genes = set([exon[2] for exon in right_hit_exons])
            #print(out_names)
            if left_hit_genes.intersection(right_hit_genes) != set():
                print(contact)
                result.append(contact)
    sys.stderr.write("block_search searching time: " + str(time.time()-begin_time) + "\n")
    return result