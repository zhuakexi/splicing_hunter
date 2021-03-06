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
    bin_index = {}
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

def build_chromsome_index(args):
    chr_name = args[0]
    p_bins = args[1]    
    bin_size = args[2]
    chromsome = args[3]
    '''
    print(chr_name)
    print(len(p_bins))
    print(bin_size)
    print(len(chromsome))
    '''
    return do_build_chromsome_index(chr_name, p_bins, bin_size, chromsome)
def do_build_chromsome_index(chr_name, p_bins:"list of bins of one chromsome", bin_size:int, chromsome:"list of exons of one chromsome") ->"exon list":
    t0 = time.time()
    p_bin_index = [[] for bin in p_bins] #initiate list of blank list same length as p_bins
    for bin in p_bins:
        exons = chromsome
        related = filt_relate_exon(bin, exons)
        the_index = bin_to_index(bin, bin_size)
        p_bin_index[the_index] = related
    sys.stderr.write("build_chromsome_index: build %s in %s \n" %( chr_name, str(time.time()-t0)) )
    return (chr_name, p_bin_index)
def key_to_index(locus, bin_size):
    return locus//bin_size
def bin_to_index(bin, bin_size):
    '''
    bin[1] can always been divided exactly by BINSIZE 
    '''
    return bin[1]//bin_size
def in_exon(contact:"line", bin_index:dict, binsize:int)->bool:
    if contact["chr1"] != contact["chr2"]:
        return False
    left_index, right_index = key_to_index(contact["pos1"], binsize), key_to_index(contact["pos2"], binsize)
    left_hit_exons = filt_in_exon(contact["pos1"], bin_index[contact["chr1"]][left_index])
    if left_hit_exons != []:
        right_hit_exons = filt_in_exon(contact["pos2"], bin_index[contact["chr2"]][right_index])
        left_hit_genes = set([exon[2] for exon in left_hit_exons])
        right_hit_genes = set([exon[2] for exon in right_hit_exons])
        return left_hit_genes.intersection(right_hit_genes) != set()
    else:
        return False
        
def block_search(bin_index:"dict of list", binsize:int, cell:"dataframe")->"data_frame":
    t0 = time.time()
    #vectorize using .pairs, target form
    mask = cell.apply(in_exon, axis=1, bin_index=bin_index, binsize=binsize)
    #print(cell[mask])
    cleaned_contacts = cell[~mask]
    sys.stderr.write("block_search searching time: %.2fs\n" % (time.time()-t0))
    return cleaned_contacts
    
