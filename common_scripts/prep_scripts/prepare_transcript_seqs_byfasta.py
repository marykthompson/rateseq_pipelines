#180119 MKT
#build salmon index for transcript mapping
#allow option of whether to include introns (e.g. include transcriptx_unspliced) or not
#the _byfasta version is going to match the seqs from a fasta file containing the transcripts
#instead of finding them by featuretye

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import gffutils

def get_transcript_boundaries(db = None, included_txt_files = None):
    '''
    Given transcript features of specified type (i.e. mRNA, ncRNA),
    get their genomic coordinates (to use for construction of _unspliced features
    in the Salmon index
    featuretypes = ['mRNA', 'ncRNA']
    '''
    genes_2_coord_dict = {}
    all_txt_ids = []
    for i in included_txt_files:
        included_txts = list(SeqIO.to_dict(SeqIO.parse(i, "fasta")))
        all_txt_ids.extend(included_txts)
        
    for i in all_txt_ids:
        #remove the transcript version number as this is not stored in the db
        txt_name = i.split('.')[0]
        region = db[txt_name]
        genes_2_coord_dict[i] = region

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

def write_intronless_file(db = None, sirv_db_file = None, include_sirvs = None, intronless_table = None, included_erccs = None, transcript_ids = None):
    '''
    Write file containing the ids of all features without introns. 
    To be used after mapping to combine the _unspliced and _spliced counts for these genes
    '''
    intronless = set()
    #get all the included transcripts from the included transcripts files:
    for this_id in transcript_ids:
        exons = list(db.children(this_id, featuretype = 'exon'))
        if len(exons) == 1:
            intronless.add(this_id)

    if include_sirvs == True:
        sirv_db = gffutils.FeatureDB(sirv_db_file)
        for region in sirv_db.features_of_type('transcript'):
            this_id = region.attributes['transcript_id'][0]
            intronless.add(this_id)

    included_erccs = set(included_erccs)
    all_intronless = intronless | included_erccs

    print('num intronless', len(all_intronless))

    with open(intronless_table, 'w') as g:
        for i in all_intronless:
            g.write('%s\n' % i)

def main():
    sirv_db_file = snakemake.input['sirv_db_file']
    db_file = snakemake.input['db_file']
    genome_fasta = snakemake.input['genome_fasta']
    
    included_transcript_files = snakemake.params['transcript_fastas']
    include_unspliced = snakemake.params['include_unspliced']
    included_erccs = snakemake.params['included_erccs']
    include_sirvs = snakemake.params['include_sirvs']

    intronless_table = snakemake.output['intronless_table']
    genomic_seqs_file = snakemake.output['genomic_seqs_file']
 
    print('hi running!')
    if include_unspliced == True:
        db = gffutils.FeatureDB(db_file)
        genes_2_coord_dict = get_transcript_boundaries(db = db, included_txt_files = included_transcript_files)
        get_unspliced_seqs(genome_fasta = genome_fasta, genomic_seqs_file = genomic_seqs_file, genes_2_coord_dict = genes_2_coord_dict)
        included_transcript_ids =  set(list(genes_2_coord_dict))
        write_intronless_file(db = db, sirv_db_file = sirv_db_file, include_sirvs = include_sirvs, included_erccs = included_erccs, intronless_table = intronless_table, transcript_ids = included_transcript_ids)
        #otherwise just create empty file to stick with the same workflow
    else:
        with open(genomic_seqs_file, 'w') as g:
            g.write('')
        with open(intronless_table, 'w') as g:
            g.write('')

if __name__ == '__main__':
    main()
    