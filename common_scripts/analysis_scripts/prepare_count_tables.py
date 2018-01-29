
# coding: utf-8

# Goals:
# - take in the quant.sf values from Salmon mapping
# - convert to RPKM and make intron and exon counts tables to feed into INSPEcT
# - write INSPEcT script for each library, next step will run these from within R
#
# Version 3:
# - combine reads from intron and exon mapping for purposes of calculating read cutoffs
#
# Version 4:
# - implement combining _spliced and _unspliced mapping counts for intronless genes
# - make different count tables either including introns or not to run the analysis separately and see what happens
#
# Version 5:
# - make an analysis version that combines reads mapping to different transcripts by gene and repeats all the analysis by gene
#
# Version 6:
# - Add the spike_in seqs to the txt-> gene dict so that they can still be incorporated downstream

#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

#2) Filter by count and make tables for input into INSPEcT

#version3,
#get time pt & labelling parameters from snakemake
#write the processing rates as an outfile

input_df = pd.read_csv(snakemake.input[0], sep = '\t', index_col = 'Name')
pd_df = pd.read_csv(snakemake.input[1], sep = '\t', index_col = 'Name')
intronless = set(pd.read_csv(snakemake.params['intronless_file'], header = None, index_col = 0).index.values)

db = gffutils.FeatureDB(snakemake.params['db_file'])
print('getting all mRNAs')
all_mRNAs = set([i.id for i in db.all_features()])
txt_2_name = {i: db[i].attributes['gene_id'][0] for i in all_mRNAs}
print('mRNAs found')
#note which rows to use for the intron counts
#then remove '_unspliced' from the name as it will be written to a new file for INSPEcT and needs to be matched
#note: do the gene names need to be in the same order for INSPEcT?
exp_dict = {'total': input_df, 'foursu': pd_df}
for exp in exp_dict:
    print('exp', exp)
    df = exp_dict[exp]
    df['intron'] = df.index.map(lambda x: x.split('_')[-1] == 'unspliced')
    df['Name'] = df.index.map(lambda x: x.split('_')[0])

    #split dfs into intron and exon mapping:
    exon_df = df[df['intron'] == False][['Name', 'TPM']].copy().set_index('Name')
    intron_df = df[df['intron'] == True][['Name', 'TPM']].copy().set_index('Name')

    print('getting exon and intron reads...')
    #set TPM_total to the sum of intron and exon reads
    df = pd.merge(exon_df, intron_df, left_index = True, right_index = True, suffixes = ('_exon', '_intron'))
    #now get the IDs not in the txt_2_name dict (like the spike-ins) and add to the dict
    all_ids = set(df.index.tolist())
    old_keys = set(txt_2_name)
    missing_genes = all_ids.difference(old_keys)
    for i in missing_genes:
        txt_2_name[i] = i

    df['TPM_total'] = df['TPM_exon'] + df['TPM_intron']

    #set TPM_exon = TPM spliced or TPM total depending on whether the gene is intronless or not
    #set TPM_intron = 0 if gene is intronless
    df['TPM_exon'] = df.apply(lambda row: row['TPM_total'] if row.name in intronless else row['TPM_exon'], axis = 1)
    df['TPM_intron'] = df.apply(lambda row: 0 if row.name in intronless else row['TPM_intron'], axis = 1)
    print('counting by gene1')
    #write output files by the summing the TPM from individual transcripts
    df['transcript'] = df.index
    df['gene_id'] = df['transcript'].map(txt_2_name)

    print('writing bygene plusI out')
    #groupby gene and sum:
    sum_by_gene = df.groupby('gene_id').sum()
    sum_by_gene[['TPM_exon']].to_csv(snakemake.output['byGene_plusI_%s_exons_file' % exp], header = False, sep = ' ')
    sum_by_gene[['TPM_intron']].to_csv(snakemake.output['byGene_plusI_%s_introns_file' % exp], header = False, sep = ' ')
    print('writing bytxt plusI out')
    #this are the values to use in the + intron analysis, write to output file:
    df[['TPM_exon']].to_csv(snakemake.output['byTxt_plusI_%s_exons_file' % exp], header = False, sep = ' ')
    df[['TPM_intron']].to_csv(snakemake.output['byTxt_plusI_%s_introns_file' % exp], header = False, sep = ' ')

    #now set all TPM_intron to 0 to exlcude intron reads for the minus intron analysis
    df['TPM_intron'] = 0
    print('counting by gene2')
    print('writing bygene minusI out')

    sum_by_gene = df.groupby('gene_id').sum()
    sum_by_gene[['TPM_exon']].to_csv(snakemake.output['byGene_minusI_%s_exons_file' % exp], header = False, sep = ' ')
    sum_by_gene[['TPM_intron']].to_csv(snakemake.output['byGene_minusI_%s_introns_file' % exp], header = False, sep = ' ')
    print('writing bytxt minusI out')

    df[['TPM_exon']].to_csv(snakemake.output['byTxt_minusI_%s_exons_file' % exp], header = False, sep = ' ')
    df[['TPM_intron']].to_csv(snakemake.output['byTxt_minusI_%s_introns_file' % exp], header = False, sep = ' ')


for method in ['byTxt', 'byGene']:
    for itype in ['minusI', 'plusI']:
        with open(snakemake.output['%s_%s_rscript' % (method, itype)], 'w') as f:
            f.write('library(INSPEcT)\n')
            #each replicate will get it's own INSPEcT analysis and output files
            f.write('su_exons= read.table("%s", row.names=1)\n' % snakemake.output['%s_%s_foursu_exons_file' % (method, itype)])
            f.write('su_introns= read.table("%s", row.names=1)\n' % snakemake.output['%s_%s_foursu_introns_file' % (method, itype)])
            f.write('t_exons= read.table("%s", row.names=1)\n' % snakemake.output['%s_%s_total_exons_file' % (method, itype)])
            f.write('t_introns= read.table("%s", row.names=1)\n' % snakemake.output['%s_%s_total_introns_file' % (method, itype)])
            f.write('mycounts<-list(foursu_exons= su_exons, foursu_introns= su_introns, total_exons= t_exons, total_introns= t_introns)\n')
            f.write('tpts<-c(0)\n')
            #f.write('tL<-0.33\n')
            f.write('tL<-%s\n' % snakemake.params['labelling_hrs'])
            f.write('myAnalysis<-newINSPEcT(tpts, tL, mycounts$foursu_exons, mycounts$total_exons, mycounts$foursu_introns, mycounts$total_introns, BPPARAM = SerialParam(), degDuringPulse=%s)\n' % snakemake.params['deg_during_pulse'])
            f.write('write.csv(ratesFirstGuess(myAnalysis, "synthesis"), "%s")\n' % snakemake.params['%s_%s_synthesis_file' % (method, itype)])
            f.write('write.csv(ratesFirstGuess(myAnalysis, "degradation"), "%s")\n' % snakemake.params['%s_%s_decay_file' % (method, itype)])
            f.write('write.csv(ratesFirstGuess(myAnalysis, "processing"), "%s")\n' % snakemake.params['%s_%s_processing_file' % (method, itype)])
