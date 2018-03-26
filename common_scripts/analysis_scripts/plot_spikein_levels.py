#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 15:27:12 2018

@author: maryk.thompson

Designed to test the correspondence between the SIRV spike-in reads and the known concentrations.
Can be used to look at the effect of different read depths or analysis methods on isoform quantification
"""

#1) Imports, this relies on utils keeping same relative path
import sys
util_dir = '../../common_scripts/pipe_utils/'
#util_dir = '../pipe_utils/' #this is path relative to /common_scripts
#change this later to get from snakemake config file maybe
sys.path.append(util_dir)
from import_file import *
import re

def plot_spikeins(big_df, genes_dict = None, outdir = None):
    '''Count number of reads mapping to spike-in sets.
    Write summary to file and also plot bar graphs
    genes_dict = {'gene_set1':{}, 'gene_set2':{}}
    Note: in future, maybe move this code to module, and import to do the same thing for rRNA mapping'''

    big_df['txt'] = big_df.index
    gene_sets = list(genes_dict)
    all_genes = set().union(*list(genes_dict.values()))

    this_outdir = os.path.join(outdir, 'spike_ins')
    if not os.path.exists(this_outdir):
        os.makedirs(this_outdir)

    #1) get component column names for each replicate, which are split into exon and intron
    colnames = [i for i in big_df.columns.values if 'counts' in i]
    d = defaultdict(list)
    for i in colnames:
        coldata = i.split('_')
        exp, rep = coldata[0], coldata[-1]
        newname = '%s_%s' % (exp, rep)
        d[newname].append(i)
    for k in d:
        big_df[k] = big_df[d[k]].sum(axis = 1)

    #2) Make dfs for each gene set and get percentages of mapping counts
    #make a dataframe for each set of genes
    df_dict = {}
    #for gene set (i.e. erccs)
    for s in genes_dict:
        df_dict[s] = big_df.loc[big_df.index.isin(genes_dict[s])].copy()
        #make summary df of spike-ins
        count_dfs = []
        percent_dfs = []
        #adjust the width of graph based on number of samples plotted:
        #for experiment (e.g. total_1)
        for exp in d:
            #get the total of all counts for each experiment
            total = big_df[exp].sum()
            outname = os.path.join(this_outdir, '%s_%s_bar') % (s, exp)
            #print('this df', df_dict[s])
            if not df_dict[s].empty:
                new_df = pipeline_aux.quick_barplot(data = df_dict[s], y_col = exp, x_col = 'txt', divide_by = total, percent = True, outname = outname, x_label = 'transcript', y_label = '% of mapping transcripts')
            else:
                new_df = df_dict[s]
                new_df['percent'] = np.nan
                new_df['frac'] = np.nan
            count_dfs.append(new_df[[exp]])
            new_df.rename(columns = {'percent':'percent_%s' % exp}, inplace = True)
            percent_dfs.append(new_df[['percent_%s' % exp]])

        #write summary output files
        summary_df = pd.concat(count_dfs + percent_dfs, axis = 1)
        summary_df.to_csv(os.path.join(this_outdir, '%s_spike_counts.csv' % s))

    #Remove ERCC and SIRV genes from dataset as we don't want to plot these for gene reproducibility
    big_df.drop(labels = all_genes, inplace = True)
    return big_df

def main():
    
    SIRV_quant_file = "/Users/maryk.thompson/Desktop/Davislab/C2.18.spike_in_stability/datasources/spike_ins/SIRV_spike_counts.csv"
    SIRV_abundance_file = "/Users/maryk.thompson/Desktop/Davislab/C2.18.spike_in_stability/notebooks/SIRV_spikein_molarities.csv"
    #datafile = snakemake.input['SIRV_quant_file']
    #abundances = snakemake.input['SIRV_abundance_file']
    
    #big_df = pd.read_csv(datafile, index_col = 'transcript')
    #get set of ercc and sirv transcripts by name
    #all_genes = big_df.indexl
    
    SIRV_quant_df = pd.read_csv(SIRV_quant_file, index_col = 'transcript')
    SIRV_ab_df = pd.read_csv(SIRV_abundance_file, index_col = 'Name')
    
    quant_cols = SIRV_quant_df.columns.values
    
    #get the columns that match total_1, etc. or foursu_1, etc.
    tot_pattern = re.compile("total_[0-9]")
    foursu_pattern = re.compile("foursu_[0-9]")

    total_cols = [i for i in quant_cols if tot_pattern.match(i) != None]
    foursu_cols = [i for i in quant_cols if foursu_pattern.match(i) != None]
    
    exps = total_cols + foursu_cols
    for i in exps:
        merged_df = pd.merge(SIRV_quant_df[[i]], SIRV_ab_df, left_index = True, right_index = True)
        pipeline_aux.scatter_plot(merged_df, x_col = 'E2', y_col = i, x_name = 'conc. (fmol/ul)', y_name = 'read counts', filename = 'test_spike_%s' % i, logbase = 10)
        #def scatter_plot(df, x_col = None, y_col = None, x_name = None, y_name = None, filename = None, logbase = None, title = None):

    #check reproducibility between foursu and total for each replicate
    for i in range(0, len(total_cols)):
        merged_df = pd.merge(SIRV_quant_df[[total_cols[i]]], SIRV_quant_df[[foursu_cols[i]]], left_index = True, right_index = True)
        pipeline_aux.scatter_plot(merged_df, x_col = total_cols[i], y_col = foursu_cols[i], x_name = 'total RNA levels (TPM)', y_name = 'foursu RNA levels (TPM)', filename = 'test_compare_%s' % i, logbase = 10)

    #print(SIRV_ab_df.head())
    #print(SIRV_quant_df.head())
    
    #print(merged_df.head())
    
    #total_cols = [i for i in quant_cols if ('total' in i) and ('scaled' not in i)]
    #foursu_cols = [i for i in quant_cols if ('foursu' in i) and ('scaled' not in i)]
    
    #print("total_cols", total_cols)
    #print("foursu_cols", foursu_cols)
    #print('quant_cols', quant_cols)
    #print('foursu_cols', foursu_cols)
    #print(SIRV_quant_df.head())
    #print(SIRV_ab_df.head())
    '''
    erccs = set([i for i in all_genes if i.startswith('ERCC')])
    sirvs = set([i for i in all_genes if i.startswith('SIRV')])
    big_df = plot_spikeins(big_df, genes_dict = {'ERCC':erccs, 'SIRV':sirvs}, outdir = snakemake.params['plot_dir'])
    plot_reproducibility(big_df, outdir = snakemake.params['plot_dir'], gene_list_dir = snakemake.params['genelist_dir'], min_reads = snakemake.params['min_reads'], min_read_type = snakemake.params['min_read_type'])
    '''
if __name__ == '__main__':
    main()