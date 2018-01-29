#180119 MKT
#build gffutils database

import gffutils

if snakemake.params['gtf_format'] == 'ensembl':
    db = gffutils.create_db(snakemake.input['gtf_file'], snakemake.output['db_file'], disable_infer_transcripts = True, disable_infer_genes = True)

#flybase gtf file needs id_spec in order to have mRNAs as ids in the db
elif snakemake.params['gtf_format'] == 'flybase':
    key_dict = {'gene':'gene_id', 'mRNA': 'transcript_id', 'ncRNA': 'transcript_id'}
    db = gffutils.create_db(snakemake.input['gtf_file'], snakemake.output['db_file'], id_spec = key_dict, disable_infer_transcripts = True, disable_infer_genes = True)