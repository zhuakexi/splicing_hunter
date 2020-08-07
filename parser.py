#v0.2 new bed parser with hierach index
import pandas as pd
import time
import sys
# .bed parser
regular_chromsome_names = ["chr" + str(i) for i in range(1,23)]
regular_chromsome_names.extend(["chrX","chrY"])
def bed_parser(bed_name, regular="off"):
    # dict of chromsomes, for each entry: list of
    # ok with .gz thanks for pandas 
    begin = time.time()
    bed = pd.read_table(bed_name, header=None)
    #split gene id
    bed[3] = bed[3].str.split(".",expand=True)
    chromsomes = {}
    #pandas group dataframe by chr name
    grouped = bed.groupby(0)
    #chr name is the key; value format: (exon_start, exon_end, gene_id)
    if regular == "on":
        chromsomes = {name:list(zip(value[1], value[2], value[3])) for name, value in grouped if name in regular_chromsome_names}
    else:
        chromsomes = {name:list(zip(value[1], value[2], value[3])) for name, value in grouped }
    sys.stderr.write("bed_parser parsing time: " + str(time.time()-begin) + "\n")
    return chromsomes, bed
# .contact parser
def con_parser(cell_name,*sample_function):
    #list of entries: chromsome, leg1, leg2; only python engine support regex sep
    time_begin = time.time()
    con = pd.read_table(cell_name,header=None,sep='[,\t]',engine='python')
    if len(sample_function) != 0:
        con = con.sample(sample_function[0])
    cell = list(zip(con[0], con[1], con[3], con[4]))
    #check parsing time
    sys.stderr.write("con_parser parsing time: " + str(time.time() - time_begin) + "\n")
    return cell
def bin_parser(file_name, regular="off"):
    '''
    read bin.bed file
    return dictionary of bins by chromsome
    '''
    time_begin = time.time()
    bins = pd.read_table(file_name, header=None)
    print("bin_parser parsing time: " + str(time.time() - time_begin) + "\n")
    grouped = bins.groupby(0) #group by chromsome names
    if regular == "on":
        return {key:value.values for key,value in grouped if key in regular_chromsome_names}
        #return [bin for bin in bins.values if bin[0] in regular_chromsome_names]
    return {key:value.values for key,value in grouped}
def size_parser(file_name)-> dict:
    '''
    read chromsome size(start, end) from faid file
    '''
if __name__ == "main":
    #bed_name = "/share/Data/ychi/genome/hg38_RefSeq.bed"
    bed_name = "hg38_RefSeq.bed"
    begin = time.time()
    chromsomes, bed = bed_parser(bed_name)
    print("parsing time: ", time.time()-begin)
    print("chromsomes(with mito etc.):", len(chromsomes))
    print( "exon number of chr1:", len(chromsomes["chr1"]) )
    print("sample entry:", chromsomes['chr1'][0])


