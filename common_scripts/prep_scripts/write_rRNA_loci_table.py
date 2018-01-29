#180117 MKT
#make table of rRNA-mapping transcripts to rRNA loci
import sys
import os
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
fasta_file = '/Users/maryk.thompson/Desktop/Davislab/C2.16/test_rRNA_txts.fa'

def write_table(fasta_file, outname):
    #rdict = defaultdict(int)
    r_dict = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        atts = record.description.split('; ')
        k = [i.split('=')[0] for i in atts]
        v = [i.split('=')[1] for i in atts]
        d = dict(zip(k,v))
        locus = d['name'].split(':')[0]
        r_dict[record.id] = locus

    df = pd.DataFrame(list(r_dict.items()), columns = ['txt', 'locus'])
    df.to_csv('%s.csv' % outname, index = False)
    
def main(arglist):
    
    fasta_file, outname = arglist
    write_table(fasta_file, outname)
    
if __name__ == '__main__':
    main(sys.argv[1:])