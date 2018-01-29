#180119 MKT
#get the spike-in seqs and convert into fasta for incoporation into mapping index
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils

#Convert provided SIRV file to simple multifasta record, remove the extra sequence that's not included in the gtf file

sirv_infasta = snakemake.input['infasta']
sirv_gtf = snakemake.input['gtf']
include_sirvs = snakemake.params['include_sirvs']
db_file = snakemake.output['sirv_db_file']
sirv_outfasta = snakemake.output['outfasta']

if include_sirvs == True:

    #make db to get SIRV annotations
    #here, leaving infer_transcripts and infer_genes is useful because only exons are listed in the file
    #it outputs 74/136 (54%), don't know what this means
    db = gffutils.create_db(sirv_gtf, db_file)
    sirv_dict = SeqIO.to_dict(SeqIO.parse(sirv_infasta, "fasta"))

    #Get sequences for spliced or unspliced SIRVs
    records = []
    for txt in db.features_of_type('transcript'):
        exons = list(db.children(txt, featuretype = 'exon', order_by = 'start'))

        spliced_seq = Seq(''.join([str(sirv_dict[txt.seqid].seq[i.start-1:i.end]) for i in exons]))
        unspliced_seq = sirv_dict[txt.seqid].seq[exons[0].start-1:exons[-1].end]

        if txt.strand == '+':
            spliced_record = SeqRecord(spliced_seq, id = txt.id)
            unspliced_record = SeqRecord(unspliced_seq, id = '%s_unspliced' % txt.id)

        else:
            spliced_record = SeqRecord(spliced_seq.reverse_complement(), id = txt.id)
            unspliced_record = SeqRecord(unspliced_seq.reverse_complement(), id = '%s_unspliced' % txt.id)

        records.append(spliced_record)
        records.append(unspliced_record)

    with open(sirv_outfasta, 'w') as g:
        SeqIO.write(records, g, 'fasta')

else:
    with open(sirv_outfasta, 'w') as g:
        g.write('')
    with open(db_file, 'w') as g:
        g.write('')