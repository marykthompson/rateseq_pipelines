
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

def plot_reproducibility(df, outdir = None, gene_list_dir = None, min_reads = None, min_read_type = None):
    '''Make scatterplots showing reproducibility between replicates and how many genes passed filter'''

    #i) log transform synthesis and decay rates, make df = subset of columns with values
    colnames = df.columns.values
    synth_cols = [i for i in colnames if 'synthesis' in i]
    deg_cols = [i for i in colnames if 'decay' in i]
    proc_cols = [i for i in colnames if 'processing' in i]
    foursu_exon_cols = [i for i in colnames if ('foursu_exon' in i) and ('scaled' not in i)]
    foursu_intron_cols = [i for i in colnames if ('foursu_intron' in i) and ('scaled' not in i)]
    total_exon_cols = [i for i in colnames if ('total_exon' in i) and ('scaled' not in i)]
    total_intron_cols = [i for i in colnames if ('total_intron' in i) and ('scaled' not in i)]
    synth_est_cols = [i for i in colnames if 'synth_est' in i]
    decay_est_cols = [i for i in colnames if 'decay_est' in i]
    decay_from_total_est_cols = [i for i in colnames if 'decay_simple_est' in i]
    proc_est_cols = [i for i in colnames if 'proc_est' in i]

    for cols in [synth_cols, deg_cols, proc_cols, foursu_exon_cols, foursu_intron_cols, total_exon_cols, total_intron_cols, synth_est_cols, decay_est_cols, decay_from_total_est_cols, proc_est_cols]:
        for j in cols:
            df['{0}_log10_{1}'.format('_'.join(j.split('_')[:-1]), j.split('_')[-1])] = df[j].apply(np.log10)

    newcolnames = df.columns.values
    synth_cols = [i for i in newcolnames if 'synthesis_rate_log10' in i]
    deg_cols = [i for i in newcolnames if 'decay_rate_log10' in i]
    proc_cols = [i for i in newcolnames if 'processing_rate_log10' in i]
    foursu_exon_cols = [i for i in newcolnames if 'foursu_exon_counts_log10' in i]
    foursu_intron_cols = [i for i in newcolnames if 'foursu_intron_counts_log10' in i]
    total_exon_cols = [i for i in newcolnames if 'total_exon_counts_log10' in i]
    total_intron_cols = [i for i in newcolnames if 'total_intron_counts_log10' in i]
    synth_est_cols = [i for i in newcolnames if 'synth_est_log10' in i]
    decay_est_cols = [i for i in newcolnames if 'decay_est_log10' in i]
    decay_simple_est_cols = [i for i in newcolnames if 'decay_simple_est_log10' in i]
    proc_est_cols = [i for i in newcolnames if 'proc_est_log10' in i]

    df = df[synth_cols + deg_cols + proc_cols + foursu_exon_cols + foursu_intron_cols + total_exon_cols + total_intron_cols + synth_est_cols + decay_est_cols + decay_simple_est_cols + proc_est_cols]
    
    to_plot_dict = {'synthesis': synth_cols, 'decay': deg_cols, 'processing': proc_cols, 'foursu_exon': foursu_exon_cols, 'foursu_intron': foursu_intron_cols, 'total_exon': total_exon_cols, 'total_intron': total_intron_cols, 'synthesis_est': synth_est_cols, 'decay_est': decay_est_cols, 'decay_simple_est': decay_simple_est_cols, 'processing_est': proc_est_cols}
  
    #ii) make mini_dfs for each filtering combination and plot
    num_exps = len(synth_cols)
    pairs = [pair for pair in itertools.combinations(range(num_exps), 2)]
    rep_comb = ['pearson_r2_rep%s_v_rep%s' % (i[1]+1, i[0]+1) for i in pairs]
    header = ['rate', 'filtering_method', 'min_reads', 'num_genes'] + rep_comb
    filter_summary_df = pd.DataFrame(columns = header)
    gene_files = [os.path.join(gene_list_dir, i) for i in os.listdir(gene_list_dir) if i.endswith('.csv')]

    this_outdir = os.path.join(outdir, 'scatterplots')
    if not os.path.exists(this_outdir):
        os.makedirs(this_outdir)
        for rate in to_plot_dict:
            os.makedirs(os.path.join(this_outdir, rate))

    for i in gene_files:
        read_type, co = os.path.basename(i).split('.csv')[0].split('_')[-2:]
        #get empty dataframe with genes that passed the cutoff
        passed_genes = pd.read_csv(i, names = ['txt'], index_col = 'txt')
        #do inner join instead of left, because we'll already have gotten rid of the ERCC stds from df but these show up in passed
        filtered_df = pd.merge(passed_genes, df, how = 'inner', left_index = True, right_index = True)
        
        #despite filtering for read count cutoff, still won't be able to get synthesis or decay rate for a lot of genes
        #iii) Plot rate reprocibility
        for rate in to_plot_dict:
            fig = plt.figure(figsize = (8,8))
            n = 1
            corr_vals = []
            for pair in pairs:
                x_name = to_plot_dict[rate][pair[0]]
                y_name = to_plot_dict[rate][pair[1]]
               
                #because np.log(0) = -inf instead of np.nan, need to convert these values
                valid_df = filtered_df.replace([np.inf, -np.inf], np.nan, inplace = False)
                valid_df.dropna(axis = 0, how = 'any', subset = [x_name, y_name], inplace = True)

                num_plotted = len(valid_df)
                corr = valid_df[x_name].corr(valid_df[y_name])
                r2_val = corr**2
                corr_vals.append(r2_val)
                ax = fig.add_subplot(num_exps - 1, num_exps - 1, n)
                ax.scatter(valid_df[x_name], valid_df[y_name], color = 'k', s = 10)
                ax.text(0.1, 0.9, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)
                ax.set_xlabel(x_name)
                ax.set_ylabel(y_name)
                n += 1

            plt.savefig(os.path.join(this_outdir, rate, '%s_%s_%s.png' % (rate, read_type, co)))
            plt.close(fig)
            result_dict = {}
            result_dict['rate'] = rate
            result_dict['filtering_method'] = read_type
            result_dict['min_reads'] = co
            result_dict['num_genes'] = num_plotted
            for j in range(0, len(rep_comb)):
                result_dict[rep_comb[j]] = corr_vals[j]
            filter_summary_df = filter_summary_df.append(result_dict, ignore_index = True)

    filter_summary_df = filter_summary_df.round(4)
    filter_summary_df.to_csv(os.path.join(outdir, 'filtering_summary.csv'), index = False)

def main():
    datafile = snakemake.input['all_data_file']
    big_df = pd.read_csv(datafile, index_col = 'transcript')
    #get set of ercc and sirv transcripts by name
    all_genes = big_df.index
    erccs = set([i for i in all_genes if i.startswith('ERCC')])
    sirvs = set([i for i in all_genes if i.startswith('SIRV')])
    big_df = plot_spikeins(big_df, genes_dict = {'ERCC':erccs, 'SIRV':sirvs}, outdir = snakemake.params['plot_dir'])
    plot_reproducibility(big_df, outdir = snakemake.params['plot_dir'], gene_list_dir = snakemake.params['genelist_dir'], min_reads = snakemake.params['min_reads'], min_read_type = snakemake.params['min_read_type'])

if __name__ == '__main__':
    main()
