'''
180119 MKT
Write file containing the intronless genes, so that this can be used after mapping to combine the
_unspliced and spliced counts for these genes
'''

import gffutils

##Figure out how to exclude SIRVs from this if not in the experiment
#I think you'd need to turn sirv_db into an optional parameter in snakemake
sirv_db_file = snakemake.input['sirv_db']
db_file = snakemake.input['db']
included_erccs = snakemake.params['included_erccs']
included_featuretypes = snakemake.params['included_featuretypes']
intronless_table = snakemake.output['intronless_table']
include_sirvs = snakemake.params['include_sirvs']

db = gffutils.FeatureDB(db_file)

intronless = set()

for i in included_featuretypes:
    for region in db.features_of_type(i):
        #seems like you can pass the gene object to get the children but you have to pass the mRNA/ncRNA id, weird
        this_id = region.attributes['transcript_id'][0]
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
