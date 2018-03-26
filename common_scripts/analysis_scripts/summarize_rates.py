
# coding: utf-8

# ## summarize_rates.py
#
# Goals:
# - take the output of the INSPEcT pipeline and create summary datasheets and graphs
# - allow filtering of included genes using either the sum of the input and pd experiments or per library
#
# Version 1:
# potential future developments:
# 1) fix unnecessary reparsing of the quant.sf files -- can you save the count_dfs from the snakemake run
# For now, I am actually reading in the output of the prepare_count_table intron and and exon tables
# This has the added advantage, that using actual written files helps snakemake construct the DAG tree and I'm not
# sure if this would be made to work with python variable names. I doubt it.
# version 3 incorporate:
# processing rates
# version 5:
# incorporate rate estimates using simple versions of formulas, +/- intron counts
# make script compatible with both intron and no intron version of pipeline
# this version has now excluded the spike_ins! Fix tis.
# version 6:
# changing to take the rates summary files that contains synth, deg, and proc.
# changing to take the counts summary file that contains intron, exon for foursu and total for each experiment

#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

def estimate_rates_simple(df, tL = None, in_pd_spike_ratio_ERCC = None, in_pd_spike_ratio_SIRV = None):
    '''Get an estimate of rates using simple formulas from dePretis_GR_2015
    Also generate scaled estimates of the read counts (not TPM, which is already scaled) using the spike-in counts
    Explanation of the spike adjustment scaling factors:
    pd_in_RNA_ratio = the ratio of total RNA used in the PD to the amount used in the input
    (e.g. for 'SPLICE' experiments, 500 ng was used for input vs. 24 ug for pulldown, therefore PD/IN = 48)
    in_pd_spike_ratio = the pipetting ratio of input/pd (i.e. ratio of spike-in amount added to tube for input vs. pd)
    if we multiple the spike ratio by the pd_in_RNA_ratio, then this will make the spike-ins to the scale
    that they would have been if we'd added at the same scale relative to the input [RNA].
    This will make the spike-in values a lot higher and therefore the genecounts a lot lower.
    For now try scaling by either ERCC or SIRV or total counts
    '''

    all_txts = set(df.index.values)

    ERCCs = set([i for i in all_txts if i.startswith('ERCC')])
    SIRVs = set([i for i in all_txts if i.startswith('SIRV')])
    allspike = ERCCs | SIRVs

    ERCC_df = df[df.index.isin(ERCCs)].copy()
    SIRV_df = df[df.index.isin(SIRVs)].copy()
    allspike_df = df[df.index.isin(allspike)].copy()

    ERCC_totals = ERCC_df.sum()
    SIRV_totals = SIRV_df.sum()
    allspike_totals = allspike_df.sum()


    #1) get hierarchical index-- equiv. to the replicates
    reps = df.columns.levels[0]

    #2) adjust pd spike-in levels by dilution factor to match input spike-in levels
    #counts total = intron + exon mapping reads, this will be used for calculating the spike-in levels
    scale_factors = {'SIRV':{}, 'ERCC':{}, 'spike':{}}
    for rep in reps:
        scale_factors['SIRV'][rep] = {}
        scale_factors['ERCC'][rep] = {}
        scale_factors['spike'][rep] = {}

        #calculate the scaling correction factor if using both ERCC and SIRV spike_ins
        #I'm using the weighted average (weighting from total SIRV or ERCC counts), is this the best way?
        allERCC_counts = ERCC_totals[rep]['counts_total_total'] + ERCC_totals[rep]['counts_total_foursu']
        allSIRV_counts = SIRV_totals[rep]['counts_total_total'] + SIRV_totals[rep]['counts_total_foursu']
        frac_ERCC = allERCC_counts/(allERCC_counts + allSIRV_counts)
        frac_SIRV = 1 - frac_ERCC
        in_pd_spike_ratio_all = (frac_ERCC*in_pd_spike_ratio_ERCC + frac_SIRV*in_pd_spike_ratio_SIRV)/2

        scale_factors['ERCC'][rep]['total'] = ERCC_totals[rep]['counts_total_total']
        scale_factors['ERCC'][rep]['foursu'] = ERCC_totals[rep]['counts_total_foursu']*in_pd_spike_ratio_ERCC

        scale_factors['SIRV'][rep]['total'] = SIRV_totals[rep]['counts_total_total']
        scale_factors['SIRV'][rep]['foursu'] = SIRV_totals[rep]['counts_total_foursu']*in_pd_spike_ratio_SIRV

        scale_factors['spike'][rep]['total'] = ERCC_totals[rep]['counts_total_total'] + SIRV_totals[rep]['counts_total_total']
        scale_factors['spike'][rep]['foursu'] = ERCC_totals[rep]['counts_total_foursu']*in_pd_spike_ratio_all + SIRV_totals[rep]['counts_total_foursu']*in_pd_spike_ratio_all

    #3) normalize counts for each gene to the spike-in counts for that library
    #note that for now I'm doing this by dividing by the number of spikein counts for the library, not by applying scaling factor around 1
    count_classes = {'foursu': {'intron': 'NumReads_intron_foursu', 'exon': 'NumReads_exon_foursu'}, 'total': {'intron': 'NumReads_intron_total', 'exon': 'NumReads_exon_total'}}
    for rep in reps:
        for fraction in count_classes:
            for region in count_classes[fraction]:
                df.loc[:, (rep, '%s_ercc_scaled' % count_classes[fraction][region])] = df.loc[:, (rep, count_classes[fraction][region])]/scale_factors['ERCC'][rep][fraction]
                df.loc[:, (rep, '%s_sirv_scaled' % count_classes[fraction][region])] = df.loc[:, (rep, count_classes[fraction][region])]/scale_factors['SIRV'][rep][fraction]
                df.loc[:, (rep, '%s_spike_scaled' % count_classes[fraction][region])] = df.loc[:, (rep, count_classes[fraction][region])]/scale_factors['spike'][rep][fraction]

    #4) calculate simple-estimated rates using the scaled counts
    for rep in reps:
        df.loc[:, (rep, 'synth_est')] = df.loc[:, (rep, '%s_sirv_scaled' % count_classes['foursu']['exon'])]/tL
        df.loc[:, (rep, 'decay_est')] = df.loc[:, (rep, 'synth_est')]/(df.loc[:, (rep, '%s_sirv_scaled' % count_classes['total']['exon'])] - df.loc[:, (rep, '%s_sirv_scaled' % count_classes['total']['intron'])])
        df.loc[:, (rep, 'decay_simple_est')] = df.loc[:, (rep, 'synth_est')]/df.loc[:, (rep, '%s_sirv_scaled' % count_classes['total']['exon'])]
        df.loc[:, (rep, 'proc_est')] = df.loc[:, (rep, 'synth_est')]/df.loc[:, (rep, '%s_sirv_scaled' % count_classes['total']['intron'])]

    return df

def write_filtered_genes(pass_filter_dict, filtering_outdir = None, min_read_type = None, min_reads = None):
    '''Given the filtered dfs for each specification, 1) concat them together by replicate,
    2) make a set of genes that pass filtering for all replicates, 3) write set to file,
    4) return list of pass filter dfs for the specified cutoffs--this will be used to make
    the main output spreadsheet'''

    for read_type in pass_filter_dict:
        for co in pass_filter_dict[read_type]:
            df = pd.concat(pass_filter_dict[read_type][co], axis = 1)
            df['all_passed_filter'] = df['pass_filter'].all(axis = 1) #double brackets preserve column names in output
            df_passed = df[df['all_passed_filter'] == True].copy()
            df_passed['txt'] = df_passed.index
            df_passed['txt'].to_csv(os.path.join(filtering_outdir, 'filtered_transcripts_%s_%s.csv' % (read_type, co)), index = False)

    return pass_filter_dict[min_read_type][min_reads] #this is a list of dfs

def reorder_columns(df):
    '''
    hacky way to reorder columns because pandas sort_index to sort level 0 (replicates)
    sorts both the first level and the second level, which I don't want
    '''
    first_level = [x.tolist() for x in df.columns.levels[:1]]
    myorder = df.columns.levels[1].tolist()
    #set pass filter as the fourth column after the rate estimates
    myorder.remove('pass_filter')
    myorder.insert(3, 'pass_filter')
    first_level.append(myorder)
    new_index = list(itertools.product(*first_level))
    df = df[new_index]

    return df

def build_df(count_files = None, rate_files = None, min_read_type = 'perlib', min_reads = 100, filtering_outdir = None, labelling_hrs = None, in_pd_spike_ratio_ERCC = None, in_pd_spike_ratio_SIRV = None):
    '''Filter data according to specifications and combine into one big spreadsheet
    Also try filtering with some different specifications and output these for plotting specifically'''

    count_dfs = [pd.read_csv(i, sep = ' ', index_col = 'Name') for i in count_files]
    rate_dfs = [pd.read_csv(i, header = 0, names = ['Name', 'synthesis', 'degradation', 'processing'], skiprows = 1, index_col = 'Name') for i in rate_files]
    #this could be made less ugly if we got R to write an index name!

    #1) Go through by replicate and apply filtering
    print('filtering genes by expression level')
    pass_filter_dict = {}
    min_read_types = ['sum', 'perlib']
    read_cos = [0.1, 1, 10, 100]
    for i in min_read_types:
        pass_filter_dict[i] = {}
        for j in read_cos:
            pass_filter_dict[i][j] = [] #list to collect all the filtered results

    #apply each of the filters by replicate
    for i in range(0, len(count_dfs)):
        count_df = count_dfs[i]
        for read_type in pass_filter_dict:
            for co in pass_filter_dict[read_type]:
                filter_df = pd.DataFrame(index = count_df.index)
                if min_read_type == 'perlib':
                    filter_df['pass_filter'] = (count_df['counts_total_total'] >= min_reads) & (count_df['counts_total_foursu'] >= min_reads)
                elif min_read_type == 'sum':
                    filter_df['pass_filter'] = (count_df['counts_total_total'] + count_df['counts_total_foursu'] >= min_reads)
                pass_filter_dict[read_type][co].append(filter_df)

    print('done filtering...')
    #write files containing list of genes that pass filter by different specifications
    #because Im controlling the file names from Python and not snakemake, need to make filtering outdir
    if not os.path.exists(filtering_outdir):
        os.makedirs(filtering_outdir)

    #return filtered_dfs for the selected filtering method (min_read_type) and the min required reads (min_reads)
    filtered_dfs = write_filtered_genes(pass_filter_dict, filtering_outdir = filtering_outdir, min_read_type = min_read_type, min_reads = min_reads)

    #2) Make a df for each replicate containing the rates, counts and filtering
    rep_dfs = []
    for i in range(0, len(filtered_dfs)):
        rep_dfs.append(pd.concat([rate_dfs[i], count_dfs[i], filtered_dfs[i]], axis = 1))

    #3) concatenate the replicate dfs into one large df with a hierarchical index containing the info about the replicate
    replist = ['rep_%s' % (i + 1) for i in range(len(filtered_dfs))]
    big_df = pd.concat(rep_dfs, axis = 1, keys = replist)

    #4) calculate the rough estimate rates using simple scaling and formulas from dePretis_GR_2015
    print('calculating simple rate estimates')
    big_df = estimate_rates_simple(big_df, tL = labelling_hrs, in_pd_spike_ratio_ERCC = in_pd_spike_ratio_ERCC, in_pd_spike_ratio_SIRV = in_pd_spike_ratio_SIRV)

    #5) reorder columns by replicate for output sheet
    big_df = reorder_columns(big_df)

    return big_df

def summarize_df(df):
    '''Take averages and write to output file for quick summary'''
    mean_df = pd.DataFrame(index = df.index)
    reps = df.columns.levels[0]
    #haven't been able to find a good built-in way of doing these averages with multiIndex, so just using list of tuples to select rep data
    mean_df['all_passed_filter'] = df[list(itertools.product(*[reps, ['pass_filter']]))].all(axis = 1)

    cols_to_average = ['synthesis', 'degradation', 'processing', 'synth_est', 'decay_est', 'decay_simple_est', 'proc_est', 'TPM_exon_foursu', 'TPM_intron_foursu', 'TPM_exon_total', 'TPM_intron_total',\
                       'NumReads_exon_foursu', 'NumReads_intron_foursu', 'NumReads_exon_total', 'NumReads_intron_total']

    for col in cols_to_average:
        data_index = list(itertools.product(*[reps, [col]]))
        mean_df['%s_mean' % col] = df[data_index].mean(axis = 1)

    mean_df['mean_halflife_INSPEcT_(min)'] = (math.log(2)/mean_df['degradation_mean'])*60
    mean_df['mean_halflife_est_(min)'] = (math.log(2)/mean_df['decay_est_mean'])*60

    return mean_df

def main():
    db = gffutils.FeatureDB(snakemake.params['annotation_db'])
    #take in a list of counts files and rates files, ordered by replicate
    count_files = snakemake.input['count_files']
    rate_files = snakemake.input['rate_files']

    min_read_type = snakemake.params['min_read_type']
    min_reads = snakemake.params['min_reads']
    filtering_outdir = snakemake.params['filtering_outdir']
    all_data_file = snakemake.output['all_data_file']
    summary_data_file = snakemake.output['summary_data_file']
    labelling_hrs = snakemake.params['labelling_hrs']
    in_pd_spike_ratio_ERCC = snakemake.params['in_pd_spike_ratio_ERCC']
    in_pd_spike_ratio_SIRV = snakemake.params['in_pd_spike_ratio_SIRV']

    #combine replicate data into one large dataframe with level 0 index for replicate
    df = build_df(count_files = count_files, rate_files = rate_files, min_read_type = min_read_type, min_reads = min_reads, filtering_outdir = filtering_outdir, labelling_hrs = labelling_hrs, in_pd_spike_ratio_ERCC = in_pd_spike_ratio_ERCC, in_pd_spike_ratio_SIRV = in_pd_spike_ratio_SIRV)

    #make summary sheet with the averages of all replicate data
    summary_df = summarize_df(df)

    #add the common gene names to each spreadsheet
    name_df = pd.DataFrame(df.index, index = df.index, columns = ['transcript'])
    all_mRNAs = set([i.id for i in db.all_features()])
    name_df['gene_id'] = name_df['transcript'].apply(lambda x: db[x].attributes['gene_id'][0] if x in all_mRNAs else np.nan)
    name_df['gene_symbol'] = name_df['transcript'].apply(lambda x: db[x].attributes['gene_symbol'][0] if x in all_mRNAs else np.nan)
    #give the name df a hierarchical index so that it will write nicely to spreadsheet
    name_df_layered = pd.concat([name_df], axis = 1, keys = ['gene_info'])

    alldata_df = pd.concat([name_df_layered, df], axis = 1)
    summary_df = pd.concat([name_df, summary_df], axis = 1)
    filtered_summary_df = summary_df[summary_df['all_passed_filter'] == True].copy()

    alldata_df.to_csv(all_data_file, index = False)
    summary_df.to_csv(summary_data_file, index = False)
    filtered_summary_df.to_csv('%s_filtered.csv' % summary_data_file.split('.csv')[0], index = False)

if __name__ == '__main__':
    main()
