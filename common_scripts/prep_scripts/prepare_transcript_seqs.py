#180119 MKT
#build salmon index for transcript mapping
#allow option of whether to include introns (e.g. include transcriptx_unspliced) or not
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import gffutils

def get_transcript_boundaries(db = None, featuretypes = None):
    '''
    Given transcript features of specified type (i.e. mRNA, ncRNA),
    get their genomic coordinates (to use for construction of _unspliced features
    in the Salmon index
    featuretypes = ['mRNA', 'ncRNA']
    '''
    genes_2_coord_dict = {}
    for i in featuretypes:
        for region in db.features_of_type(i):
            region_id=region.attributes['transcript_id'][0]
            genes_2_coord_dict[region_id] = region

    return genes_2_coord_dict

def get_unspliced_seqs(genome_fasta = None, genomic_seqs_file = None, genes_2_coord_dict = None):
    '''Get the genomic sequences for transcripts, corresponding to unspliced transcripts'''
    
    #1) Get unspliced sequences for these features
    f = gzip.open(genome_fasta, 'rt')
    #gzip open 'rt' = text mode
    #I couldn't get the .index() function to work with the file handle method
    chr_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    f.close()
    #4) Get the genomic region sequences, convert to SeqRecord objects and write to outfile
    records = []
    for i in genes_2_coord_dict:
        region = genes_2_coord_dict[i]
        if region.strand == '+':
            record = SeqRecord(chr_dict[region.seqid].seq[region.start-1:region.end], id = '%s_unspliced' % i)
        else:
            record = SeqRecord(chr_dict[region.seqid].seq[region.start-1:region.end].reverse_complement(), id = '%s_unspliced' % i)
        records.append(record)

    with open(genomic_seqs_file, 'w') as g:
        SeqIO.write(records, g, 'fasta')

def main():
    include_unspliced = snakemake.params['include_unspliced']
    if include_unspliced == True:
        db = gffutils.FeatureDB(snakemake.input['db_file'])
        genes_2_coord_dict = get_transcript_boundaries(db = db, featuretypes = snakemake.params['included_featuretypes'])
        get_unspliced_seqs(genome_fasta = snakemake.input['genome_fasta'], genomic_seqs_file = snakemake.output['genomic_seqs_file'], genes_2_coord_dict = genes_2_coord_dict)
        #otherwise just create empty file to stick with the same workflow
    else:
        with open(snakemake.output['genomic_seqs_file'], 'w') as g:
            g.write('')

if __name__ == '__main__':
    main()
    