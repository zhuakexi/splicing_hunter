#v0.3 read file name, better out format,  new reference file
#v0.2 faster parser using hierarchical index
import time
import sys
import pandas as pd
import numpy as np

def filt_in_exon(locus, exons):
    '''
    get list of exons envelope locus
    '''
    return [exon for exon in exons if exon[0] <= locus <= exon[1]]
def back_search(gene_name, df):
    '''
    search exons from gene name
    '''
    #get target records
    short_df = df[df[3]==gene_name]
    left, right, names = short_df[1], short_df[2], short_df[3]
    return list(zip(left,right,names)) 

#--------calculate hit--------
def filter_search(cell, chromsomes, ref_df):
    time_begin = time.time()
    result = []
    records = []
    left_hit = [] # to test filter method
    for i in cell:
        chrom_name, locus1, locus2 = i[0], i[1], i[2] 
        locus1_in = filt_in_exon(locus1, chromsomes[chrom_name])
        if locus1_in != []:
            #print("left :",locus1_in)
            left_hit.append(i)
            gene_names = np.unique(np.array(locus1_in)[:,2])
            #print("gene names:", gene_names)
            for j in gene_names:
                j_exons = back_search(j, ref_df)
                hit_right = filt_in_exon(locus2, j_exons)
                if hit_right != []:
                    hit_left = filt_in_exon(locus1, j_exons)
                    hit_contact = str((chrom_name, locus1, locus2))
                    hit_text = ":".join(map(str, (hit_left, hit_right) ))
                    record = "----->".join( (hit_contact, hit_text) )
                    result.append(hit_contact)
                    records.append(record)
    print("Filter_search left hit",len(left_hit)) # to test filter method
    sys.stderr.write( "Filter_search hit: " + str( len(set(result)) ) + "\n"  )
    sys.stderr.write( "Filter_search time: " + str(time.time() - time_begin) +"\n" )

