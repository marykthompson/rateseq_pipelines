#module for storing RNA-Seq pipeline auxiliary functions
import sys
#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

def quick_barplot(data = None, y_col = None, x_col = None, divide_by = 1, percent = False, outname = None, x_label = None, y_label = None):
    '''
    Make quick barplot to summarize count data,
    for example from the spike-in or rRNA-mapping reads
    Input:
    a dataframe (df) containing y_col values to plot, x_col ids
    divide_by = a number which all values in y_col will be divided by
    percent = convert fraction to percent for plotting
    '''
    #normalize values, this will auto change the col_to_plot
    df = data
    col_to_plot = y_col
    if divide_by != 1:
        df['frac'] = df[y_col]/divide_by
        col_to_plot = 'frac'

    if percent == True:
        df['percent'] = df['frac']*100
        col_to_plot = 'percent'

    height = 4
    width_ratio = len(df)/16
    width = height*width_ratio
    if width < 4:
        width = 4
    fig = plt.figure(figsize = (width, height))
    ax = fig.add_subplot(111)
    ax = sns.barplot(data = data, x = x_col, y = col_to_plot, ax = ax)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    #just using text.set_rotation(45) doesn't allow you to align wrt the axis
    xlabels = df[x_col]
    #lining the right side of text box up with tick looks best of options, still not great
    ax.set_xticklabels(xlabels, rotation = 45, ha = 'right')
    plt.tight_layout()
    plt.savefig('%s.png' % outname)
    plt.close(fig)

    return df

def plot_genomic_region(coverage, chrom, start, end, strand, positions = None):
    '''
    Given a genomic region and coverage HTSeq GA, plot reads mapping to that region
    [start, end) 0-based
    Also and option to pass a postion file which will then mark regions that overlap with those positions
    '''
    window = HTSeq.GenomicInterval(chrom, start, end, strand)
    wincvg = np.fromiter(coverage[window], dtype='i')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(range(start, end), wincvg)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_ylim(bottom=0)
    ax.set_xlim(start, end)
    plt.xticks(rotation=45)
    ax.set_ylabel('read counts')
    ax.set_xlabel('position')

    return ax

def limit_df(df, colnames = None, min_val = None, max_val = None):
    '''
    Filter out data that is outside given min and max values
    '''

    for col in colnames:
        df = df[df[col].between(min_val, max_val, inclusive = True)].copy()

    df.dropna(subset = colnames, inplace = True)

    return df

def clean_df(df, cols = None):
    '''
    Replace inf values with nan & drop nan-containing rows over given columns
    This needs to be done before analysis and/or after log-transform
    '''
    df.replace([np.inf, -np.inf], np.nan, inplace = True)
    df.dropna(axis = 0, how = 'any', subset = cols, inplace = True)

def cdf_plot(df, x_col = None, bg_group = 'all', group_col = None, nbins = 100, x_label = None, y_label = None, filename = None, logbase = None, title = None):
    '''
    Plot data in given column and a group_by variable as CDF
    x_col = the name of the column containing the data
    group_col = the name of the column containing the classification variable, e.g. short, long
    nbins = number of bins to use
    '''
    fig = plt.figure(figsize = (8, 8))
    ax = fig.add_subplot(111)

    groups = df[group_col].unique()
    groups = sorted(groups)

    #bring bg group to front if we are not using a specific bg group
    if bg_group != 'all':
        groups.insert(0, groups.pop(groups.index(bg_group)))
        first_df = None
        #build first df from first group
    #otherwise the first group will be the whole dataset, named 'all'
    else:
        groups.pop(groups.index('bg'))
        groups.insert(0, 'all')
        first_df = df

    all_bins = []
    all_cdfs = []
    handles = []
    labels = []

    for i in range(0, len(groups)):
        if i == 0:
            if first_df is not None:
                sub_df = first_df
            else:
                sub_df = df[df[group_col] == groups[i]].copy()
        else:
            sub_df = df[df[group_col] == groups[i]].copy()

        data = sub_df[x_col]
        counts, bin_edges = np.histogram(data, bins = nbins, normed = True)
        cdf = np.cumsum(counts)
        all_cdfs.append(cdf)
        all_bins.append(bin_edges)
        l, = plt.plot(bin_edges[1:], cdf/cdf[-1])
        handles.append(l)

        if i == 0:
            if bg_group != 'all':
                bg_data = data
            else:
                bg_data = df[x_col]
            labels.append('{:}, n={:}'.format(groups[i], len(data)))
        else:
            stat, pval = sp.stats.ks_2samp(bg_data, data)
            labels.append('{:}, n={:}, p={:.2E}'.format(groups[i], len(data), Decimal(pval)))

    plt.legend(handles, labels)

    #ax.text(0.1, 0.9, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes, fontsize = 12)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    plt.title(title)
    if filename != None:
        plt.savefig('%s.png' % filename)
        plt.close(fig)
    else:
        return ax

def seaborn_box():
	'''Seaborn style box plot'''
	#boxplot(df, xname, yname)
	sns.set(style='ticks')
	fig=plt.figure(figsize=(8,4))
	ax=sns.boxplot(data=data, showfliers=False)

	#sns.boxplot(x=xname, y=yname, data=df, showfliers=False)

	#sns.boxplot(x=xname, y=yname, hue=colname, data=df, showfliers=False)
	#sns.despine(offset=10, trim=True)
	sns.despine(offset=10)
	return ax

def scatter_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', title = ''):
    '''
    Make a multiscatter plot of all the combinations in given columns
    Also retun correlations for each
    '''
    #store results of each replicate as correlation dict
    corr_dict = {}
    num_plotted = len(df)
    fig = plt.figure(figsize = (8,8))
    xname = cols[0]
    yname = cols[1]
    corr = df[xname].corr(df[yname])
    r2_val = corr**2

    ax = fig.add_subplot(111)
    ax.scatter(df[xname], df[yname], color = 'k', s = 10)
    ax.text(0.1, 0.8, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)

    label_axes(ax, xname = xname, yname = yname, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix)

    return {'num_plotted': num_plotted, 'fig': fig}

def multiscatter_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', filename = None, title = ''):
    '''
    Make a multiscatter plot of all the combinations in given columns
    Also retun correlations for each
    '''
    #store results of each replicate as correlation dict
    corr_dict = {}
    num_plotted = len(df)
    pairs =  [pair for pair in itertools.combinations(range(len(cols)), 2)]
    fig = plt.figure(figsize = (8,8))
    n = 1
    for pair in pairs:
        xi = pair[0]
        yi = pair[1]
        xname = cols[pair[0]]
        yname = cols[pair[1]]
        corr = df[xname].corr(df[yname])
        r2_val = corr**2
        corr_dict['%s_v_%s' % (yname[label_index_level], xname[label_index_level])] = r2_val
        ax = fig.add_subplot(len(pairs) - 1, len(pairs) - 1, n)
        ax.scatter(df[xname], df[yname], color = 'k', s = 10)
        ax.text(0.1, 0.8, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)
        ax.set_xlabel('%s %s' % (xname[label_index_level], axis_title_suffix))
        ax.set_ylabel('%s %s' % (yname[label_index_level], axis_title_suffix))
        n += 1

    return {'ax': ax, 'corr_dict': corr_dict, 'num_plotted': num_plotted, 'fig': fig}

def log_transform(df, cols = None, logbase = None, label_index_level = None):
    '''Transform columns with given logbase and return new df and column names'''
    if logbase in np_log_transforms:
        logcols = []
        for col in cols:
            #if col is a tuple, we're dealing with a hierarchical index
            #for now I'm assuming the next level (label_index_level + 1) will have the name of the data to transform

            #if label_index_level = 0, indicates that we are plotting different reps against each other
            #if label_index_level = 1, indicates that we are plotting different columns from same rep against each other
            #maybe this should be made more flexible in the future
            if type(col) == tuple:
                newcol = (col[0], '%s_log' % col[1])
            else:
                newcol = '%s_log' % col

            df[newcol] = df[col].apply(np_log_transforms[logbase])
            logcols.append(newcol)
    else:
        raise NotImplementedError('logbase %s not supported' % logbase)

    #shouldn't need to retun df, as it should be modified here
    return logcols

def label_axes(ax, xname = None, yname = None, label_index_level = None, axis_title_suffix = ''):
    '''
    label x and y-axes of plot
    '''
    #if it's a multiIndex df, we'd like to specify which level to use,
    #but if it's a single index, then this is going to slice the string which is not what we want
    if type(xname) == tuple:
        ax.set_xlabel('%s %s' % (xname[label_index_level], axis_title_suffix))
        ax.set_ylabel('%s %s' % (yname[label_index_level], axis_title_suffix))
    else:
        ax.set_xlabel('%s %s' % (xname, axis_title_suffix))
        ax.set_ylabel('%s %s' % (yname, axis_title_suffix))

def save(fig, filename = None, title = None, figformat = 'png'):
    '''
    Save figure
    '''
    plt.tight_layout()
    plt.suptitle(title)
    #plt.subplots_adjust(top = 0.92)
    plt.savefig('%s.%s' % (filename, figformat))
    plt.close(fig)

def plot(df, cols = None, plottype = None, logbase = None, title = '', label_index_level = 0, axis_title_suffix = '', filename = None, figformat = 'png'):
    '''
    Given a df and plottype, clean up data and then send to plotting function
    df = pandas dataframe, cols = list of columns with data, plottype = 'scatter', etc.
    '''
    #make a copy of the df here. All subsequent operations will modify this copy
    df = df[cols].copy()
    clean_df(df, cols = cols)

    #get new log-transformed columns if required
    if logbase != None:
        cols = log_transform(df, cols = cols, logbase = logbase, label_index_level = label_index_level)
        clean_df(df, cols = cols)

    #send to plotting fxn, which will return plot-specific analyses in the results dict
    results = plot_fxn_dict[plottype](df, cols = cols, title = title, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix)
    save(results['fig'], filename = filename, title = title, figformat = figformat)
    return results

plot_fxn_dict = {'multiscatter':multiscatter_plot, 'scatter':scatter_plot}
np_log_transforms = {2:  np.log2, 10: np.log10}
