
#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

#version5 will take the list of genes passing filter
#make:
#1) summary file with correlations and number of genes plotted
#2) plots containing each set of filtered genes
#3) plots and summary file with spike-in mapping info --> now move the barplot code into module file
#version 7: plot tpm reproducibility, not just rate reproducibility
#version 8: fix it to use the rep hierarchical indexing sheet

def plot_spikeins(df, genes_dict = None, outdir = None):
    '''
    Count number of reads mapping to spike-in sets.
    Write summary to file and also plot bar graphs
    genes_dict = {'gene_set1':{}, 'gene_set2':{}}
    Return the df with the spike-ins removed so that they won't be included in the gene reproducibility plots.
    '''

    #get all spikein gene names
    gene_sets = list(genes_dict)
    all_genes = set().union(*list(genes_dict.values()))

    this_outdir = os.path.join(outdir, 'spike_ins')
    os.makedirs(this_outdir, exist_ok = True)

    #get rep hierarchical index
    reps = [i for i  in df.columns.levels[0] if i.startswith('rep')]
    to_plot = [i for i in df.columns.levels[1] if 'counts' in i]

    #make plot and summary table for each set of spikeins (ERCC, SIRV)
    for s in genes_dict:
        count_dfs = []
        percent_dfs = []
        for cat in to_plot:
            for rep in reps:
                outname = os.path.join(this_outdir, '_'.join([s, cat, rep, 'bar']))
                total = df[(rep, cat)].sum()
                data = df.loc[df.index.isin(genes_dict[s]), rep]
                #add index to column for easy passing to seaborn
                data['txt'] = data.index
                expname = '%s_%s' % (cat, rep)
                #if there are reads mapped to spikeins, then make plots and output files
                if not data.empty:
                    new_df = pipeline_aux.quick_barplot(data = data, y_col = cat, x_col = 'txt', divide_by = total, percent = True, outname = outname, x_label = 'transcript', y_label = '% of mapping transcripts')
                else:
                    new_df = data
                    new_df['percent'] = np.nan
                    new_df['frac'] = np.nan

                new_df.rename(columns = {'percent': 'percent_%s' % expname, 'frac': 'frac_%s' % expname}, inplace = True)

                count_dfs.append(new_df[['frac_%s' % expname]])
                percent_dfs.append(new_df[['percent_%s' % expname]])

        #write summary output files
        summary_df = pd.concat(count_dfs + percent_dfs, axis = 1)
        summary_df.to_csv(os.path.join(this_outdir, '%s_spike_counts.csv' % s))

    #Remove ERCC and SIRV genes from dataset as we don't want to plot these for gene reproducibility
    df.drop(labels = all_genes, inplace = True)

def plot_reproducibility(df, outdir = None, gene_list_dir = None, min_reads = None, min_read_type = None):
    '''Make scatterplots showing reproducibility between replicates and how many genes passed filter'''

    #plot reproducibility at different filtering levels for each of the following
    to_plot = ['synthesis', 'degradation', 'processing', 'TPM_exon_foursu', 'TPM_intron_foursu', 'TPM_exon_total', 'TPM_intron_total']
    reps = [i for i in df.columns.levels[0] if i.startswith('rep')]
    gene_files = [os.path.join(gene_list_dir, i) for i in os.listdir(gene_list_dir) if i.endswith('.csv')]

    this_outdir = os.path.join(outdir, 'scatterplots')
    os.makedirs(this_outdir, exist_ok = True)
    for rate in to_plot:
        os.makedirs(os.path.join(this_outdir, rate), exist_ok = True)

    corr_results_dict = {}
    n = 0
    for i in gene_files:
        read_type, co = os.path.basename(i).split('.csv')[0].split('_')[-2:]
        #get dataframe with only the subset of genes that passed filter
        passed_genes = set(pd.read_csv(i, names = ['txt'], index_col = 'txt').index.tolist())
        filtered_df = df.loc[df.index.isin(passed_genes)]
        for rate in to_plot:
            #get location of each replicate from the hierachical index
            these_cols = list(itertools.product(*[reps, [rate]]))
            exp_name = '%s_%s_%s' % (rate, read_type, co)
            filename = os.path.join(this_outdir, rate, '%s' % exp_name)

            title = exp_name
            #give index level 0 so that it will write 'rep_1', etc. in the axis titles
            label_index_level = 0
            axis_title_suffix = '%s (log10)' % rate
            results = pipeline_aux.plot(filtered_df, cols = these_cols, plottype = 'multiscatter', logbase = 10, title = exp_name, label_index_level = 0, axis_title_suffix = axis_title_suffix, filename = filename)
            exp_d = {'rate': rate, 'filtering_method': read_type, 'min_reads': co, 'num_genes': results['num_plotted']}
            exp_d.update(results['corr_dict'])
            corr_results_dict[n] = exp_d
            n += 1

    #write summary file with number of genes passes and repllicate correlations
    filter_summary_df = pd.DataFrame.from_dict(corr_results_dict, orient = 'index')
    filter_summary_df = filter_summary_df.round(4)
    filter_summary_df.to_csv(os.path.join(outdir, 'filtering_summary.csv'), index = False)

def main():
    datafile = snakemake.input['all_data_file']

    df = pd.read_csv(datafile, header = [0,1])
    df.set_index(('gene_info', 'transcript'), inplace = True, drop = False)
    df.index.rename('transcript', inplace = True)

    #get set of ercc and sirv transcripts by name
    all_genes = df.index
    erccs = set([i for i in all_genes if i.startswith('ERCC')])
    sirvs = set([i for i in all_genes if i.startswith('SIRV')])
    #running plot_spikeins() first removes the spikein genes from df, so they don't end up in the reproducibility scatter plots
    plot_spikeins(df, genes_dict = {'ERCC':erccs, 'SIRV':sirvs}, outdir = snakemake.params['plot_dir'])
    plot_reproducibility(df, outdir = snakemake.params['plot_dir'], gene_list_dir = snakemake.params['genelist_dir'], min_reads = snakemake.params['min_reads'], min_read_type = snakemake.params['min_read_type'])

if __name__ == '__main__':
    main()
