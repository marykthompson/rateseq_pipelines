
#1) Imports, this relies on utils keeping same relative path
import sys
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *
import argparse

#1) summary file with correlations and number of genes plotted
#2) plots containing each set of filtered genes
#3) plots and summary file with spike-in mapping info
#4) fraction recovery of endogenous mRNAs and ERCCs compared to input
## TODO:
# Get mpl to autoscale axes to include all points but also be square and the same across replicates

def plot_recovery(df, outdir = None, spike_sets = ['SIRV', 'ERCC'], intron_plot = True):
    '''
    1) Plot fraction of transcripts unspliced vs. spliced for input and foursu libraries
    2) Plot scaled foursu TPM vs. scaled total TPM to look at recovery by gene
    '''
    this_outdir = os.path.join(outdir, 'foursu_recovery')
    os.makedirs(this_outdir, exist_ok = True)

    reps = [i for i  in df.columns.levels[0] if i.startswith('rep')]
    cols = df.columns.levels[1]

    spike_names = {'SIRV': 'sirv_scaled', 'ERCC': 'ercc_scaled'}
    spike_strings = [spike_names[k] for k in spike_sets]

    for rep in reps:
        data = df[(rep)].copy() #this selects the data for the replicate
        #using TPM for now, but could change to counts if desired
        m = 'TPM'
        for s in spike_strings:
            l = [i for i in cols if i.startswith(m) and i.endswith(s)]
            foursu_exon = [i for i in l if 'exon' in i and 'foursu' in i][0]
            input_exon = [i for i in l if 'exon' in i and 'foursu' not in i][0]
            foursu_intron = [i for i in l if 'intron' in i and 'foursu' in i][0]
            input_intron = [i for i in l if 'intron' in i and 'foursu' not in i][0]

            data['%s_inputall_%s' % (m, s)] = data[input_exon] + data[input_intron]
            data['%s_foursuall_%s' % (m, s)] = data[foursu_exon] + data[foursu_intron]

            foursu_all_sum = data['%s_foursuall_%s' % (m, s)].sum()
            input_all_sum = data['%s_inputall_%s' % (m, s)].sum()

            frac_intron_foursu = data[foursu_intron].sum()/foursu_all_sum
            frac_intron_input = data[input_intron].sum()/input_all_sum

            percent_foursu = (foursu_all_sum/input_all_sum)*100

            #plot recovery: scaled TPM of foursu vs. scaled TPM of total library
            outname = os.path.join(this_outdir, '_'.join(['recovery', s, rep]))
            results = pipeline_aux.plot(data, cols = ['%s_inputall_%s' % (m, s), '%s_foursuall_%s' % (m, s)], plottype = 'regplot', logbase = 10, labels = {'ylabel':'scaled foursu RNA levels (log10 TPM)', 'xlabel':'scaled total RNA levels (log10 TPM)'}, text = ['% input recovered in foursu: {:02.2f}'.format(percent_foursu), 0.1, 0.8], ax_loc = 1.0, limits = {'x':[-11, -1], 'y':[-11, -1]}, xy_line = True, filename = outname)

            if intron_plot == True:
                sum_dict = {'foursu':{'intron': frac_intron_foursu, 'total': 1}, 'input': {'intron': frac_intron_input, 'total': 1}}

                sum_df = pd.DataFrame.from_dict(sum_dict, orient = 'index')
                sum_df['library'] = sum_df.index

                #order the plot to have input, then foursu
                new_index = ['input', 'foursu']
                sum_df = sum_df.reindex(new_index)

                #plot fraction of unspliced vs. spliced transcripts
                outname = os.path.join(this_outdir, '_'.join(['intron_stackedbar', s, rep]))
                results = pipeline_aux.plot(sum_df, cols = ['library', 'total', 'intron'], cat_labels = ['spliced', 'unspliced'], plottype = 'stacked_bar', filename = outname, labels = {'ylabel':'fraction TPM'})

def plot_spikeins(df, genes_dict = None, outdir = None, sirv_mol_file = None, in_pd_ERCC = None, pd_in_RNA = None):
    '''
    Count number of reads mapping to spike-in sets.
    Write summary to file and also plot bar graphs
    genes_dict = {'gene_set1':{}, 'gene_set2':{}}
    Return the df with the spike-ins removed so that they won't be included in the gene reproducibility plots.
    '''
    SIRV_ab_df = pd.read_csv(sirv_mol_file, index_col = 'Name')

    #get all spikein gene names
    gene_sets = list(genes_dict)
    all_genes = set().union(*list(genes_dict.values()))

    this_outdir = os.path.join(outdir, 'spike_ins')
    os.makedirs(this_outdir, exist_ok = True)

    #get rep hierarchical index
    reps = [i for i  in df.columns.levels[0] if i.startswith('rep')]
    #get the colnames corresponding to the unscaled counts and TPM vals
    count_cols = [i for i in df.columns.levels[1] if 'counts' in i]
    tpm_cols = [i for i in df.columns.levels[1] if 'TPM_total' in i]

    #because I'm using the word total to equal both total RNA and intron+exon, need to do this
    #perhaps later change to eq instead of total?
    total_count_cols = [i for i in count_cols if 'foursu' not in i]
    foursu_count_cols = [i for i in count_cols if 'foursu' in i]

    total_tpm_cols = [i for i in tpm_cols if 'foursu' not in i]
    foursu_tpm_cols = [i for i in tpm_cols if 'foursu' in i]

    count_exps = total_count_cols + foursu_count_cols
    tpm_exps = total_tpm_cols + foursu_tpm_cols

    #make plot and summary table for each set of spikeins (ERCC, SIRV)
    for s in genes_dict:
        count_dfs = []
        percent_dfs = []
        for rep in reps:
            data = df.loc[df.index.isin(genes_dict[s]), rep]
            #prevent breaking if no spikein data present
            if data.empty:
                continue
            #1) and #2) For each library==>
            #1) plot counts vs. known molar spikein concentration for each lib.
            for i in count_exps:
                expname = '%s_%s' % (i, rep)
                #only make this plot for SIRVs, don't know [ERCC]
                if s == 'SIRV':
                    outname = os.path.join(this_outdir, '_'.join([s, i, rep, 'byMol']))
                    merged_df = pd.merge(data, SIRV_ab_df, left_index = True, right_index = True)
                    pipeline_aux.plot(merged_df, cols = ['E2', i], plottype = 'scatter', logbase = 10, labels = {'ylabel':'read counts', 'xlabel':'conc. (fmol/ul)'}, filename = outname)

            #2) Plot fraction spikeins in each library
            for i in tpm_exps:
                outname = os.path.join(this_outdir, '_'.join([s, i, rep, 'tpm_bar']))
                total = df[(rep, i)].sum()
                pipeline_aux.plot(data.reset_index(), cols = ['transcript', i], plottype = 'quickbar', filename = outname, labels = {'ylabel':'% of TPM', 'xlabel':'transcript'}, divide_by = total, percent = True)
            for i in count_exps:
                outname = os.path.join(this_outdir, '_'.join([s, i, rep, 'bar']))
                total = df[(rep, i)].sum()
                #if there are reads mapped to spikeins, then make plots and output files
                new_df = pipeline_aux.plot(data.reset_index(), cols = ['transcript', i], plottype = 'quickbar', filename = outname, labels = {'ylabel':'% of mapping reads', 'xlabel':'transcript'}, divide_by = total, percent = True)['df']
                #because of the hack to get it to use the index as xlabel, need to set here
                new_df.set_index('transcript', inplace = True)
                new_df.rename(columns = {'percent': 'percent_%s' % expname, 'frac': 'frac_%s' % expname}, inplace = True)

                count_dfs.append(new_df[['frac_%s' % expname]])
                percent_dfs.append(new_df[['percent_%s' % expname]])

            #3) For each matched set of foursu and total ==>
            #3) check reproducibility between foursu and total for each replicate
            for i in range(0, len(total_tpm_cols)):
                #will generally only be one graph per rep unless we subdivide more
                outname = os.path.join(this_outdir, '_'.join([s, str(i), rep, 'compare']))
                pipeline_aux.plot(data, cols = [total_tpm_cols[i], foursu_tpm_cols[i]], plottype = 'scatter', logbase = 10, labels = {'ylabel':'foursu RNA levels (log10 TPM)', 'xlabel':'total RNA levels (log10 TPM)'}, filename = outname)

        #write summary output files
        summary_df = pd.concat(count_dfs + percent_dfs, axis = 1)
        summary_df.to_csv(os.path.join(this_outdir, '%s_spike_counts.csv' % s))

    #4) Plot ERCC recovery, scaled by SIRVs
    ercc_df = df.loc[df.index.isin(genes_dict['ERCC'])].copy()
    adj_ratio = in_pd_ERCC * pd_in_RNA
    for rep in reps:
        ercc_df.loc[:, (rep, 'TPM_exon_foursu_sirv_scaled')] = df.loc[:, (rep, 'TPM_exon_foursu_sirv_scaled')] * adj_ratio
        ercc_df.loc[:, (rep, 'TPM_intron_foursu_sirv_scaled')] = df.loc[:, (rep, 'TPM_intron_foursu_sirv_scaled')] * adj_ratio

        plot_recovery(ercc_df, outdir = this_outdir, spike_sets = ['ERCC'], intron_plot = False)

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

def main(arglist):

    if not 'snakemake' in globals():
        parser = argparse.ArgumentParser()
        parser.add_argument('-datafile', help = 'csv file with results from 3 reps')
        parser.add_argument('-plot_dir', default = 'test_plots')
        parser.add_argument('-genelist_dir', help = 'dir containing list of genes passing filter')
        parser.add_argument('-min_reads', default = 100)
        parser.add_argument('-min_read_type', default = 'perlib')
        parser.add_argument('-SIRV_mols', default = "/Users/maryk.thompson/Desktop/Davislab/C2.18.spike_in_stability/notebooks/SIRV_spikein_molarities.csv")
        parser.add_argument('-in_pd_ERCC', default = 9)
        parser.add_argument('-pd_in_RNA', default = 48)
        args = parser.parse_args(args = arglist)
        datafile, plot_dir, gene_list_dir, min_reads, min_read_type, SIRV_mols, in_pd_ERCC, pd_in_RNA = args.datafile, args.plot_dir, args.gene_list_dir, args.min_reads, args.min_read_type, args.SIRV_mols, args.in_pd_ERCC, args.pd_in_RNA
        #can we make more compact? https://docs.python.org/3/library/argparse.html
    else:
        datafile = snakemake.input['all_data_file']
        plot_dir = snakemake.params['plot_dir']
        genelist_dir = snakemake.params['genelist_dir']
        min_reads = snakemake.params['min_reads']
        min_read_type = snakemake.params['min_read_type']
        #ADD THESE AS SNAKEMAKE PARAMS!!!
        SIRV_mols = snakemake.params['SIRV_molarity_file']
        in_pd_ERCC = snakemake.params['in_pd_spike_ratio_ERCC']
        pd_in_RNA = snakemake.params['pd_in_RNA_ratio']

    df = pd.read_csv(datafile, header = [0,1])
    df.set_index(('gene_info', 'transcript'), inplace = True, drop = False)
    df.index.rename('transcript', inplace = True)

    #get set of ercc and sirv transcripts by name
    all_genes = df.index
    erccs = set([i for i in all_genes if i.startswith('ERCC')])
    sirvs = set([i for i in all_genes if i.startswith('SIRV')])
    #running plot_spikeins() first removes the spikein genes from df, so they don't end up in the reproducibility scatter plots
    plot_spikeins(df, genes_dict = {'ERCC':erccs, 'SIRV':sirvs}, outdir = plot_dir, sirv_mol_file = SIRV_mols, in_pd_ERCC = in_pd_ERCC, pd_in_RNA = pd_in_RNA)
    plot_recovery(df, outdir = plot_dir)
    plot_reproducibility(df, outdir = plot_dir, gene_list_dir = genelist_dir, min_reads = min_reads, min_read_type = min_read_type)

'''
#local test:
python plot_reproducibility.py -genelist_dir /Users/maryk.thompson/Desktop/Davislab/comp_labbook_backup/data/computational_projects/C2/C2.15.SnakeMakeSetup/pipeline_testing/byGene/plusI/filtered_genes/ -datafile /Users/maryk.thompson/Desktop/Davislab/comp_labbook_backup/data/computational_projects/C2/C2.15.SnakeMakeSetup/pipeline_testing/byGene/plusI/inspect_data_full.csv
'''
if __name__ == '__main__':
    main(sys.argv[1:])
