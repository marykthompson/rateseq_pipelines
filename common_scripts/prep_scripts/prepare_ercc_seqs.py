#180119 MKT
#get the spike-in seqs and convert into fasta for incoporation into mapping index

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools

#1) Select the used ERCC transcripts and write to output file, also write another version called ercc-xxx_unspliced to match the endogenous transcript file
ercc_infile = snakemake.input['infasta']
ercc_outfile = snakemake.output['outfasta']

#Just set the included_erccs to [] if no erccs in the experiment, then should write empty file
included_erccs = snakemake.params['included_erccs']

records = list(SeqIO.parse(ercc_infile, "fasta"))
spliced_records = (record for record in records if record.id in included_erccs)
unspliced_records = (SeqRecord(record.seq, id = '%s_unspliced' % record.id) for record in records if record.id in included_erccs)
all_ercc_records = itertools.chain(unspliced_records, spliced_records)
SeqIO.write(all_ercc_records, ercc_outfile, "fasta")