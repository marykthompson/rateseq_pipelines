#180123 MKT
#remove really long unspliced genomic regions and corresponding spliced transcripts
#necessary to keep total mouse genome size < 2^31 and allow Salmon index to build


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def remove_long(infasta, outfasta, max_len):
    '''Remove the spliced and unspliced versions of the same transcript if >200k long'''
    record_dict = SeqIO.index(infasta, "fasta")
    ids = set(list(record_dict))

    #1) find transcripts that are too long. Many of these will be the _unspliced versions
    long_txts = set()
    for i in ids:
        length = len(record_dict[i].seq)
        if length > max_len:
            long_txts.add(i.split('_unspliced')[0])

    #2) make new set that contains both the _unspliced and normal names,
    unspliced = set([i+'_unspliced' for i in long_txts])
    to_remove = long_txts | unspliced
    to_keep = ids.difference(to_remove)

    #3) Write transcripts that aren't too long to new file
    with open(outfasta, 'w') as g:
        for i in to_keep:
            SeqIO.write(record_dict[i], g, "fasta")

def main(args):
    infasta = snakemake.input['infasta']
    outfasta = snakemake.output['outfasta']
    max_len = snakemake.params['max_length']
    #infasta, outfasta = args
    remove_long(infasta, outfasta, max_len)
if __name__ == '__main__':
    main(sys.argv[1:])
    #main()