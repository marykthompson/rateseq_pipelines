
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
#
# Version 7:
# - prepare count tables that keep both the counts and the TPM output for each experiment. Change Rscript to handle this.

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

with open(snakemake.params['txt_2_gene_pickle'], 'rb') as f:
    txt_2_name = pickle.load(f)

#note which rows to use for the intron counts
#then remove '_unspliced' from the name as it will be written to a new file for INSPEcT and needs to be matched
#note: do the gene names need to be in the same order for INSPEcT?
exp_dict = {'total': {'infile': input_df, 'out': None}, 'foursu': {'infile': pd_df, 'out': None}}
for exp in exp_dict:
    print('exp', exp)
    df = exp_dict[exp]['infile']
    df['intron'] = df.index.map(lambda x: x.split('_')[-1] == 'unspliced')
    df['Name'] = df.index.map(lambda x: x.split('_')[0])

    #split dfs into intron and exon mapping:
    exon_df = df[df['intron'] == False][['Name', 'TPM', 'NumReads']].copy().set_index('Name')
    intron_df = df[df['intron'] == True][['Name', 'TPM', 'NumReads']].copy().set_index('Name')

    print('getting exon and intron reads...')
    #merge intron and exon reads into df
    df = pd.merge(exon_df, intron_df, left_index = True, right_index = True, suffixes = ('_exon', '_intron'))
    
    #get the transcript IDs not in the txt_2_name dict (like the spike-ins) and add to the name dict
    all_ids = set(df.index.tolist())
    old_keys = set(txt_2_name)
    missing_genes = all_ids.difference(old_keys)
    for i in missing_genes:
        txt_2_name[i] = i

    df['TPM_total'] = df['TPM_exon'] + df['TPM_intron']
    df['counts_total'] = df['NumReads_exon'] + df['NumReads_intron']

    #set TPM_exon = TPM spliced or TPM total depending on whether the gene is intronless or not
    #set TPM_intron = 0 if gene is intronless
    df['TPM_exon'] = df.apply(lambda row: row['TPM_total'] if row.name in intronless else row['TPM_exon'], axis = 1)
    df['TPM_intron'] = df.apply(lambda row: 0 if row.name in intronless else row['TPM_intron'], axis = 1)
    
    #save new df after counting intron and exon reads
    exp_dict[exp]['out'] = df
    
    
    '''
    #add gene ids from transcript IDs
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

    '''

#combine foursu and total dfs:
df = pd.merge(exp_dict['foursu']['out'], exp_dict['total']['out'], left_index = True, right_index = True, suffixes = ('_foursu', '_total'))

#get gene names by transcript ID
df['transcript'] = df.index
df['gene_id'] = df['transcript'].map(txt_2_name)

#groupby gene and sum for the by gene analyses
#write one output file 'by gene' and one output file 'by transcript'
sum_by_gene_df = df.groupby('gene_id').sum()
sum_by_gene_df.index.names = ['Name']
cols2write = sum_by_gene_df.columns.values
df[cols2write].to_csv(snakemake.output['byTxt_plusI_file'], sep = ' ')
sum_by_gene_df[cols2write].to_csv(snakemake.output['byGene_plusI_file'], sep = ' ')

#now set all TPM_intron and NumReads_intron to 0 to exlcude intron reads for the minus intron analysis
df[['TPM_intron_foursu', 'TPM_intron_total', 'NumReads_intron_foursu', 'NumReads_intron_total']] = 0

sum_by_gene_df = df.groupby('gene_id').sum()
sum_by_gene_df.index.names = ['Name']
df[cols2write].to_csv(snakemake.output['byTxt_minusI_file'], sep = ' ')
sum_by_gene_df[cols2write].to_csv(snakemake.output['byGene_minusI_file'], sep = ' ')

#write Rscript to run INSPEcT on the data
for method in ['byTxt', 'byGene']:
    for itype in ['minusI', 'plusI']:
        with open(snakemake.output['%s_%s_rscript' % (method, itype)], 'w') as f:
            f.write('library(INSPEcT)\n')
            f.write('df<-read.table("%s", header = TRUE, row.names=1)\n' % snakemake.output['%s_%s_file' % (method, itype)])
            f.write('su_exons<-df["TPM_exon_foursu"]\n')
            f.write('su_introns<-df["TPM_intron_foursu"]\n')
            f.write('t_exons<-df["TPM_exon_total"]\n')
            f.write('t_introns<-df["TPM_intron_total"]\n')
            f.write('mycounts<-list(foursu_exons= su_exons, foursu_introns= su_introns, total_exons= t_exons, total_introns= t_introns)\n')
            f.write('tpts<-c(0)\n')
            f.write('tL<-%s\n' % snakemake.params['labelling_hrs'])
            f.write('myAnalysis<-newINSPEcT(tpts, tL, mycounts$foursu_exons, mycounts$total_exons, mycounts$foursu_introns, mycounts$total_introns, BPPARAM = SerialParam(), degDuringPulse=%s)\n' % snakemake.params['deg_during_pulse'])
            f.write('syn<-ratesFirstGuess(myAnalysis, "synthesis")\n')
            f.write('deg<-ratesFirstGuess(myAnalysis, "degradation")\n')
            f.write('proc<-ratesFirstGuess(myAnalysis, "processing")\n')
            f.write('sum_df<-cbind(syn, deg, proc)\n')
            f.write('write.csv(sum_df, "%s")\n' % snakemake.params['%s_%s_rate_file' % (method, itype)])