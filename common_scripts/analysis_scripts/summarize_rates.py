
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

#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

def relabel_rep_df(df, base = None):
    '''Rename the columns after concatenation'''
    old_labels = df.columns.values
    #if no base string provided, infer from the old labels
    if base == None:
        base = old_labels[0]
    new_labels = ['%s_%s' % (base, i+1) for i in range(len(old_labels))]
    df.columns = new_labels
    return df

def filter_df(input_df, pd_df, min_read_type, min_reads):
    '''Filter the df according to specifications
    Several filters utlimately will be applied and list of genes passing
    output to file'''

    #store boolean values of filtering operation
    #I'm making a df instead of a series to use the df methods, but then I need to pass the 'pass_filter' column specifically
    filter_df = pd.DataFrame(input_df.index, index = input_df.index)
    
    if min_read_type == 'perlib':
        filter_df['pass_filter'] = (input_df['total'] >= min_reads) & (pd_df['total'] >= min_reads)
    elif min_read_type == 'sum':
        filter_df['pass_filter'] = (input_df['total'] + pd_df['total'] >= min_reads)

    return filter_df[['pass_filter']].copy()

def estimate_rates_simple(df, tL = None):
    '''Get an estimate of rates using simple formulas from dePretis_GR_2015'''

    foursu_exon_cols = ['foursu_exon_counts_1', 'foursu_exon_counts_2', 'foursu_exon_counts_3']
    foursu_intron_cols = ['foursu_intron_counts_1', 'foursu_intron_counts_2', 'foursu_intron_counts_3']
    total_exon_cols = ['total_exon_counts_1', 'total_exon_counts_2', 'total_exon_counts_3']
    total_intron_cols = ['total_intron_counts_1', 'total_intron_counts_2', 'total_intron_counts_3']

    num = len(foursu_exon_cols)
    lib_sizes = {}
    for i in range(0, num):
        foursu_lib_size = df[foursu_exon_cols[i]].sum() + df[foursu_intron_cols[i]].sum()
        total_lib_size = df[total_exon_cols[i]].sum() + df[total_intron_cols[i]].sum()
        lib_sizes['foursu_%s' % i] = foursu_lib_size
        lib_sizes['total_%s' % i] = total_lib_size

    av_size = np.average(list(lib_sizes.values()))
    scale_factors = {k: av_size/lib_sizes[k] for k in lib_sizes}
    #scale counts by average lib size
    for i in range(0, num):
        df['%s_scaled' % foursu_exon_cols[i]] = df[foursu_exon_cols[i]]*scale_factors['foursu_%s' % i]
        df['%s_scaled' % foursu_intron_cols[i]] = df[foursu_intron_cols[i]]*scale_factors['foursu_%s' % i]
        df['%s_scaled' % total_exon_cols[i]] = df[total_exon_cols[i]]*scale_factors['total_%s' % i]
        df['%s_scaled' % total_intron_cols[i]] = df[total_intron_cols[i]]*scale_factors['total_%s' % i]
        
    cols = df.columns.values
    #3) calculate rates
    foursu_exon_cols_scaled = [i for i in cols if (i.startswith('foursu_exon_counts') and i.endswith('scaled'))]
    foursu_intron_cols_scaled = [i for i in cols if (i.startswith('foursu_intron_counts') and i.endswith('scaled'))]
    total_exon_cols_scaled = [i for i in cols if (i.startswith('total_exon_counts') and i.endswith('scaled'))]
    total_intron_cols_scaled = [i for i in cols if (i.startswith('total_intron_counts') and i.endswith('scaled'))]

    #labeling time in hr, this is what I have been giving INSPEcT previously
    #1) calculate synthesis rates
    num = len(foursu_exon_cols)
    synth_cols = []
    decay_cols = []
    decay_fromtotal_cols = []
    proc_cols = []
    for i in range(0, num):
        df['synth_est_%s' % i] = df[foursu_exon_cols_scaled[i]]/tL
        df['decay_est_%s' % i] = df['synth_est_%s' % i]/(df[total_exon_cols_scaled[i]] - df[total_intron_cols_scaled[i]])
        df['decay_simple_est_%s' % i] = df['synth_est_%s' % i]/df[total_exon_cols_scaled[i]]
        df['proc_est_%s' % i] = df['synth_est_%s' % i]/df[total_intron_cols_scaled[i]]

    #add 1 to each of the est rep names
    est_cols = [i for i in df.columns.values if 'est' in i]
    new_names = {}
    for col in est_cols:
        thisname = '_'.join(col.split('_')[0:-1])
        num = int(col.split('_')[-1])
        new_names[col] = '%s_%s' % (thisname, num + 1)
        
    df.rename(columns = new_names, inplace = True)        
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

def build_df(foursu_exons = None, foursu_introns = None, total_exons = None, total_introns = None, synthesis_files = None, decay_files = None, processing_files = None, ann_db = None, min_read_type = 'perlib', min_reads = 100, filtering_outdir = None, labelling_hrs = None):
    '''Filter data according to specifications and combine into one big spreadsheet
    Also try filtering with some different specifications and output these for plotting specifically'''

    print('label', labelling_hrs)
    #convert datalist into strings for using as pandas column names
    datatype_dict = {foursu_exons: 'foursu_exons', foursu_introns: 'foursu_introns', total_exons: 'total_exons', total_introns: 'total_introns', synthesis_files: 'synthesis_files', decay_files: 'decay_files', processing_files: 'processing_files'}

    all_count_files = [foursu_exons, foursu_introns, total_exons, total_introns]
    all_summary_files = [synthesis_files, decay_files, processing_files]

    #0) Read all files in and convert to pandas dfs:
    all_dfs = {}
    for file_list in all_count_files:
        these_dfs = [pd.read_csv(i, sep = ' ', header=None, names = ['transcript', 'counts'], index_col = 'transcript')\
                     for i in file_list]
        all_dfs[datatype_dict[file_list]] = these_dfs
    for file_list in all_summary_files:
        these_dfs = [pd.read_csv(i, header = 0, names = ['transcript', 'rate'], skiprows = 1, index_col = 'transcript')\
                    for i in file_list]
        all_dfs[datatype_dict[file_list]] = these_dfs

    #1) Go through by replicate and apply filtering
    print('filtering genes by expression level')
    pass_filter_dict = {}
    min_read_types = ['sum', 'perlib']
    read_cos = [50, 100, 200, 400, 800]
    for i in min_read_types:
        pass_filter_dict[i] = {}
        for j in read_cos:
            pass_filter_dict[i][j] = [] #list to collect all the filtered results

    #apply each of the filters by replicate
    for i in range(0, len(foursu_exons)):
        foursu_exon_df = all_dfs['foursu_exons'][i]
        foursu_intron_df = all_dfs['foursu_introns'][i]
        total_exon_df = all_dfs['total_exons'][i]
        total_intron_df = all_dfs['total_introns'][i]

        #combine intron and exon reads
        input_df = pd.merge(total_exon_df, total_intron_df, left_index = True, right_index = True, suffixes = ('_exon', '_intron'))
        pd_df = pd.merge(foursu_exon_df, foursu_intron_df, left_index = True, right_index = True, suffixes = ('_exon', '_intron'))
        input_df['total'] = input_df['counts_exon'] + input_df['counts_intron']
        pd_df['total'] = pd_df['counts_exon'] + pd_df['counts_intron']

        #for now I will not allow only input or pd to pass filter, all or nothing for each replicate
        #therefore, temporarily store filtering info in input dict
        for read_type in pass_filter_dict:
            for co in pass_filter_dict[read_type]:
                filtered_df = filter_df(input_df, pd_df, read_type, co)
                pass_filter_dict[read_type][co].append(filtered_df)

    print('done filtering...')
    #write files containing list of genes that pass filter by different specifications
    #because Im controlling the file names from Python and not snakemake, need to make filtering outdir
    if not os.path.exists(filtering_outdir):
        os.makedirs(filtering_outdir)

    these_filtered_dfs = write_filtered_genes(pass_filter_dict, filtering_outdir = filtering_outdir, min_read_type = min_read_type, min_reads = min_reads)

    all_dfs['pass_filter'] = these_filtered_dfs
        #get set of genes which pass filter
        #passed_features = set(input_df[input_df['pass_filter']==True].index.tolist())

    #2) Make large dataframe, including summary columns
    #2i) Make each individual by concating all the replicates
    data_order = ['synthesis_files', 'decay_files', 'processing_files', 'pass_filter', 'total_exons', 'total_introns', 'foursu_exons', 'foursu_introns']
    base_dict = {'synthesis_files': 'synthesis_rate', 'decay_files': 'decay_rate', 'processing_files': 'processing_rate', 'total_exons': 'total_exon_counts',\
                 'total_introns': 'total_intron_counts', 'foursu_exons': 'foursu_exon_counts',\
                 'foursu_introns': 'foursu_intron_counts', 'pass_filter': 'pass_filter'}
    data_to_combine = []
    for i in data_order:
        summary_df = pd.concat(all_dfs[i], axis =1)
        summary_df = relabel_rep_df(summary_df, base = base_dict[i])
        data_to_combine.append(summary_df)

    big_df = pd.concat(data_to_combine, axis = 1)

    #calculate the rough estimate rates using simple scaling and formulas from dePretis_GR_2015
    big_df = estimate_rates_simple(big_df, tL = labelling_hrs)

    #get the corresponding gene names and symbols to add to the table
    name_df = pd.DataFrame(big_df.index, index = big_df.index, columns = ['transcript'])
    all_mRNAs = set([i.id for i in ann_db.all_features()])
    name_df['gene_id'] = name_df['transcript'].apply(lambda x: ann_db[x].attributes['gene_id'][0] if x in all_mRNAs else np.nan)
    name_df['gene_symbol'] = name_df['transcript'].apply(lambda x: ann_db[x].attributes['gene_symbol'][0] if x in all_mRNAs else np.nan)

    final_df = pd.concat([name_df, big_df], axis = 1)
    return final_df

def summarize_df(final_df):
    '''Take avearges and write to output file for quick summary'''

    #get column names for ones we want to take average of
    colnames = final_df.columns.values
    #print('original col names', colnames)
    synth = [i for i in colnames if 'synthesis' in i]
    #need to keep as 'decay_rate' otherwise will pick up decay_est, etc.
    deg = [i for i in colnames if 'decay_rate' in i]
    proc = [i for i in colnames if 'processing' in i]
    pass_filter = [i for i in colnames if 'filter' in i]
    total_exon = [i for i in colnames if 'total_exon' in i]
    total_intron = [i for i in colnames if 'total_intron' in i]
    foursu_exon = [i for i in colnames if 'foursu_exon' in i]
    foursu_intron = [i for i in colnames if 'foursu_intron' in i]
    synth_est = [i for i in colnames if 'synth_est' in i]
    decay_est = [i for i in colnames if 'decay_est' in i]
    decay_simple_est = [i for i in colnames if 'decay_simple_est' in i]
    proc_est = [i for i in colnames if 'proc_est' in i]

    #print('decay est cols', decay_est)
    #print('decay simple est cols', decay_simple_est)
    
    num_reps = len(total_exon)

    #caclulate whether each replicate passes the filter. Also get means
    #TO FIX: I THINK THIS IS CURRENTLY TAKING THE AVERAGE EVEN IF NO DATA FOR SOME EXPERIMENTS
    #PASS_FILTER IS NOT A STRINGENT WAY TO CHECK AS THERE CAN BE GENES THAT PASS BUT DON'T HAVE RATES
    final_df['all_passed_filter'] = final_df[pass_filter].all(axis = 1)
    final_df['mean_synthesis_rate'] = final_df[synth].mean(axis = 1)
    final_df['mean_decay_rate'] = final_df[deg].mean(axis = 1)
    final_df['mean_processing_rate'] = final_df[proc].mean(axis = 1)
    final_df['mean_halflife_INSPEcT_(min)'] = (math.log(2)/final_df['mean_decay_rate'])*60
    final_df['mean_total_reads'] = final_df[total_exon + total_intron].sum(axis = 1)/num_reps
    final_df['mean_foursu_reads'] = final_df[foursu_exon + foursu_intron].sum(axis = 1)/num_reps
    final_df['mean_synth_est'] = final_df[synth_est].mean(axis = 1)
    final_df['mean_decay_est'] = final_df[decay_est].mean(axis = 1)
    final_df['mean_halflife_est_(min)'] = (math.log(2)/final_df['mean_decay_est'])*60
    final_df['mean_decay_simple_est'] = final_df[decay_simple_est].mean(axis = 1)
    final_df['mean_proc_est'] = final_df[proc_est].mean(axis = 1)

    #reorder columns
    reordered_cols = ['transcript', 'gene_id', 'gene_symbol', 'all_passed_filter', 'mean_synthesis_rate', 'mean_decay_rate',\
                 'mean_synth_est', 'mean_decay_est', 'mean_decay_simple_est', 'mean_proc_est', 'mean_halflife_INSPEcT_(min)',  'mean_halflife_est_(min)', 'mean_total_reads', 'mean_foursu_reads']\
                      + synth + deg + proc + synth_est + decay_est + decay_simple_est + proc_est
    #print('reordered_cols', reordered_cols)
    summary_df = final_df[reordered_cols]
    #print('summarydf cols', summary_df.columns.values)
    return summary_df

def main():

    db = gffutils.FeatureDB(snakemake.params['annotation_db'])

    foursu_exons = snakemake.input["foursu_exons_files"]
    foursu_introns = snakemake.input["foursu_introns_files"]
    total_exons = snakemake.input["total_exons_files"]
    total_introns = snakemake.input["total_introns_files"]

    synthesis_files = snakemake.input['synthesis_files']
    decay_files = snakemake.input['decay_files']
    processing_files = snakemake.input['processing_files']
    min_read_type = snakemake.params['min_read_type']
    min_reads = snakemake.params['min_reads']
    filtering_outdir = snakemake.params['filtering_outdir']
    all_data_file = snakemake.output['all_data_file']
    summary_data_file = snakemake.output['summary_data_file']
    labelling_hrs = snakemake.params['labelling_hrs']

    #PUT THIS CHECK IN TEMPORARILY SO THAT DON'T NEED TO REBUILD FILE DURING DEBUG
    #I SHOULD PROBABLY SOON BREAK THIS UP INTO 2 SCRIPTS: 1) FILTER AND 2) SUMMARIZE
    if not os.path.exists(all_data_file):
        print('file does not exist, building...')
        print('outfile path', all_data_file)
        big_df = build_df(foursu_exons = foursu_exons, foursu_introns = foursu_introns, total_exons = total_exons, total_introns = total_introns, synthesis_files = synthesis_files, decay_files = decay_files, processing_files = processing_files, min_read_type = min_read_type, min_reads = min_reads, ann_db = db, filtering_outdir = filtering_outdir, labelling_hrs = labelling_hrs)
        big_df.to_csv(all_data_file, index = False)
        print(big_df.columns.values)
    else:
        big_df = pd.read_csv(all_data_file)
        print('file exisits')
    
    if not os.path.exists(summary_data_file):
        summary_df = summarize_df(big_df)
        summary_df.to_csv(summary_data_file, index = False)
        filtered_summary_df = summary_df[summary_df['all_passed_filter'] == True].copy()
        filtered_file = '%s_filtered.csv' % summary_data_file.split('.csv')[0]
        filtered_summary_df.to_csv(filtered_file, index = False)
    else:
        summary_df = pd.read_csv(summary_data_file)
    #plot_reproducibility(big_df, outdir = snakemake.params['plot_dir'])
    

if __name__ == '__main__':
    main()
