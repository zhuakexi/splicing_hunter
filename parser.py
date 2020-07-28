#v0.2 new bed parser with hierach index
import pandas as pd
import time
'''
#bed parser
ref = pd.read_table("/share/Data/ychi/genome/hg38_RefSeq.bed",header=None)
bed = {}
getlist = lambda i, df: list(zip(df[df[0]==i][1],df[df[0]==i][2]))
for i in bed[0].unique():
  chromsomes[i] = getlist(i, bed)
#con parser
cell_name = "GM001.impute.con.gz"
con = pd.read_table(cell_name,header=None,sep='[,\t]')
cell = list(zip(con[0],con[1],con[4]))
'''
# bed file parser
def bed_parser(bed_name, chr_names=[]):
    # dict of chromsomes, for each entry: list of 
    # two level index: chromsome, left end
    bed = pd.read_table(bed_name, header=None, index_col=[0,1])
    #split gene id
    bed[3] = bed[3].str.split(".",expand=True)
    chromsomes = {}
    #create keys from indices
    if chr_names == []:
        all_names = dict(bed.index).keys()
    else:
        all_names = chr_names
    for chr_name in all_names: 
        #key: chr_name, values: list of entries with each entry: exon start, exon end, gene id
        #bed.loc[chr_name]: choose all rows of chr_name, .index is left leg since level 1(chromsome_name) has been used.
        chromsomes[chr_name] = list(zip(bed.loc[chr_name].index ,bed.loc[chr_name][2], bed.loc[chr_name][3]))
    return chromsomes, bed
# .contact parser
def con_parser(cell_name):
    #list of entries: chromsome, leg1, leg2
    time_begin = time.time()
    con = pd.read_table(cell_name,header=None,sep='[,\t]')
    cell = list(zip(con[0],con[1],con[4]))
    #check parsing time
    #print("parsing cell:", time.time() - time_begin)
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


