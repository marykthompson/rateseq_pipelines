#MKT 18019
#find rRNA transcripts from the miscRNA.fasta file on Flybase

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip

def is_rRNA(record):
    fields = record.description.split(' ')
    txt_id = fields[0]
    atts = [i.rstrip(';') for i in fields[1:]]
    att_dict = dict(zip([i.split('=')[0] for i in atts], [i.split('=')[1] for i in atts]))
    if att_dict['type'] == 'rRNA':
        return True
    else:
        return False

def main():
    infasta = snakemake.input['infasta']
    outfasta = snakemake.output['outfasta']
    
    with gzip.open(infasta, 'rt') as f:
        records = list(SeqIO.parse(f, "fasta"))

    rRNA_records = (record for record in records if is_rRNA(record))
    SeqIO.write(rRNA_records, outfasta, "fasta")

if __name__ == '__main__':
    main()