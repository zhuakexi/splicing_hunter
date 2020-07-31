#v0.2 new bed parser with hierach index
import pandas as pd
import time
import sys
# .bed parser
def bed_parser(bed_name):
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
    cell = list(zip(con[0],con[1],con[4]))
    #check parsing time
    sys.stderr.write("con_parser parsing time: " + str(time.time() - time_begin) + "\n")
    return cell
if __name__ == "main":
    #bed_name = "/share/Data/ychi/genome/hg38_RefSeq.bed"
    bed_name = "hg38_RefSeq.bed"
    begin = time.time()
    chromsomes, bed = bed_parser(bed_name)
    print("parsing time: ", time.time()-begin)
    print("chromsomes(with mito etc.):", len(chromsomes))
    print( "exon number of chr1:", len(chromsomes["chr1"]) )
    print("sample entry:", chromsomes['chr1'][0])


