import sys
from parser import bed_parser, con_parser 
from mid_search import mid_search
ref_name = sys.argv[1]
cell_name = sys.argv[2]
out_name = sys.argv[3]

try:
    out_name = sys.argv[3]
except IndexError:
    out_name = cell_name.replace(".contacts.pairs.gz", ".rs_can")

# --------parsing files--------
## in: ref_name, cell_name, out_name
## out: *chromsomes, *ref_df, *cell
chromsomes, ref_df = bed_parser(ref_name)
cell = con_parser(cell_name)
# ------------sort exons according to mid-point------------
## in: chromsomes
## out: keys_of_all, *chromsomes
keys_of_all = {}
for i in chromsomes:
    chromsomes[i].sort(key=lambda exon : (int(exon[0])+int(exon[1])) /2)
    keys_of_all[i] = [ (int(exon[0]) + int(exon[1]))/2 for exon in chromsomes[i]]
# ------------do the search------------
result = mid_search(cell, keys_of_all, chromsomes)
# ------------write output file------------
## in: result
## out: fd
title = "splicint_hunter v0.1 " + cell_name + " in " + ref_name
content = map(str, result)
content = [i.strip("()") for i in content]
content = "\n".join(content)
with open(out_name,"w") as f:
    f.write(title+"\n"+content)
sys.stderr.write("Total contacts:" + str( len(ref_df)) + "\n")
sys.stderr.write("Contacts hit: " + str( len(result) )+ " ratio: " + str( len(result)/len(ref_df)) +"\n")
