#module for storing RNA-Seq pipeline auxiliary functions
import sys
#1) Imports, this relies on utils keeping same relative path
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *
import matplotlib.ticker as plticker

def get_fasta(infasta, outfasta, write_all = True, get_chr = None):
    '''
    Given a multifasta, split into individual files.
    Can be useful for testing tools.
    If get_chr != None, will write a matching fasta. Otherwise will write all.
    #also see notebook C2.19c
    '''
    records = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

    for k in records:
        records[k].name = ''
        records[k].description = ''

    with open(outfasta, 'w') as g:
        for k in records:
            SeqIO.write(records[k], g, 'fasta')

def test_plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(range(0, 10), range(0, 10), color = 'k', s = 10)
    #ax.text(0.1, 0.8, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)
    plt.savefig('testplot.png')
    plt.close(fig)

    #label_axes(ax, xname = xname, yname = yname, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix)
    return fig
    #return {'num_plotted': num_plotted, 'fig': fig, 'ax': ax}


def remove_spikeins(df, spikenames = ['ERCC', 'SIRV']):
    '''
    Remove all rows with index starting with spikein names, e.g. ERCC and SIRV
    '''
    all_txts = set(df.index.values)
    allspike = set()
    for s in spikenames:
        geneset =  set([i for i in all_txts if i.startswith(s)])
        allspike = allspike | geneset

    #Remove ERCC and SIRV genes from dataset as we don't want to plot these for gene reproducibility
    df.drop(labels = allspike, inplace = True)

def quick_barplot(df, cols = None, label_index_level = None, axis_title_suffix = '', title = '', limits = None, ticklabels = None, **kwargs):
    '''
    Make quick barplot to summarize count data,
    for example from the spike-in or rRNA-mapping reads
    Input:
    a dataframe (df) containing y_col values to plot, x_col ids
    divide_by = a number which all values in y_col will be divided by
    percent = convert fraction to percent for plotting
    '''
    #normalize values, this will change the col_to_plot
    x_col, y_col = cols[0:]
    col_to_plot = y_col
    if kwargs['divide_by'] != 1:
        df['frac'] = df[y_col]/kwargs['divide_by']
        col_to_plot = 'frac'

    if kwargs['percent'] == True:
        df['percent'] = df['frac']*100
        col_to_plot = 'percent'

    height = 4
    width_ratio = len(df)/16
    width = height*width_ratio
    if width < 4:
        width = 4
    fig = plt.figure(figsize = (width, height))
    ax = fig.add_subplot(111)
    ax = sns.barplot(data = df, x = x_col, y = col_to_plot, ax = ax)
    #just using text.set_rotation(45) doesn't allow you to align wrt the axis
    xlabels = df[x_col]
    #lining the right side of text box up with tick looks best of options, still not great
    ax.set_xticklabels(xlabels, rotation = 45, ha = 'right')
    plt.tight_layout()

    return {'df':df, 'fig': fig, 'ax': ax}

def stacked_bar(df, cols = None, label_index_level = None, axis_title_suffix = '', title = '', limits = None, ticklabels = None, **kwargs):
    #this scales the width of the bars but the legend is still plotted on top of the bars
    #normalize values, this will change the col_to_plot
    x_col = cols[0]
    y_cols = cols[1:]

    width = 4
    height = 4
    fig = plt.figure(figsize = (width, height))
    ax = fig.add_subplot(111)

    sns.barplot(data = df, x = x_col, y = y_cols[0], color = 'red', ax = ax)
    #need to pass 'cat_labels' to kwargs in order of cols
    bars = [plt.Rectangle((0,0),1,1,fc="red", edgecolor = 'none')]
    for i in range(1, len(y_cols)):
        bottom_plot = sns.barplot(data = df, x = x_col, y = y_cols[i], color = 'blue', ax = ax)
        bars.append(plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor = 'none'))

    l = plt.legend(bars, kwargs['cat_labels'], loc = (1.04, 0.75), prop={'size':16})
    l.draw_frame(False)
    return {'fig': fig, 'ax': ax, 'bars':bars, 'extra_artists': [l]}

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
    #print('cleaning')
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

def barhist_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', filename = None, title = '', limits = None, ticklabels = None):
    '''
    Make histogram with bars
    I think this is currently only compatible with plotting a single histogram on the axis
    '''
    #mu = np.median(to_plot_x)
    #sigma = np.std(to_plot_x)
    normed = False

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(111)
    # the histogram of the data
    ##n, bins, patches = plt.hist(df[cols], 100, normed=normed, facecolor='blue', alpha=0.8, histtype = 'bar')

    # add a 'best fit' line
    #y = mlab.normpdf( bins, mu, sigma)
    nbins = 100
    ax.hist(df[cols], bins = nbins)
    ##ax = plt.plot(bins, patches, 'r--', linewidth=1)
    return {'ax': ax}

def barhist_plot2(df, cols = None, label_index_level = None, axis_title_suffix = '', title = ''):
    '''
    Make histogram with bars
    I think this is currently only compatible with plotting a single histogram on the axis
    '''
    #mu = np.median(to_plot_x)
    #sigma = np.std(to_plot_x)

    # the histogram of the data
    n, bins, patches = plt.hist(df[cols], 100, normed=normed, facecolor='blue', alpha=0.8, histtype = 'bar')

    # add a 'best fit' line
    #y = mlab.normpdf( bins, mu, sigma)

    ax = plt.plot(bins, y, 'r--', linewidth=1)
    return {'ax': ax}

def seaborn_box(df, cols = None, label_index_level = None, axis_title_suffix = '', title = '', limits = None, ticklabels = None):
    '''Seaborn style box plot'''
    sns.set(style = 'ticks')
    fig = plt.figure(figsize = (8,8))
    ax = sns.boxplot(data = df[cols])
    fig = ax.get_figure()
    sns.despine(offset = 10)

    return {'ax':ax, 'fig': fig}

def scatter_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', title = '', limits = None, ticklabels = None, **kwargs):
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
    ax.text(0.1, 0.9, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)

    #label_axes(ax, xname = xname, yname = yname, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix)

    return {'num_plotted': num_plotted, 'fig': fig, 'ax': ax}

def reg_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', title = '', limits = None, ticklabels = None, **kwargs):
    '''
    Make a scatterplot using seaborn's regplot
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
    ax = sns.regplot(data = df, x = xname, y = yname, fit_reg = False)
    #ax.scatter(df[xname], df[yname], color = 'k', s = 10)
    ax.text(0.1, 0.9, 'r2 = %1.3f\nn = %s' % (r2_val, num_plotted), transform = ax.transAxes)
    set_lim(ax, limits = limits)

    #label_axes(ax, xname = xname, yname = yname, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix)

    return {'num_plotted': num_plotted, 'fig': fig, 'ax': ax}


def multiscatter_plot(df, cols = None, label_index_level = None, axis_title_suffix = '', filename = None, title = '', limits = None, ticklabels = None):
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
        set_lim(ax, limits = limits)
        #put this in temporarily to see if it will draw it now:
        #fig.canvas.draw()
        set_ticklabels(ax, fig, ticklabels =  ticklabels)

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
    If passed with only xname or yname, probably a 1D plot, only label that axis
    '''
    #if it's a multiIndex df, we'd like to specify which level to use,
    #but if it's a single index, then this is going to slice the string which is not what we want
    if type(xname) == tuple:
        ax.set_xlabel('%s %s' % (xname[label_index_level], axis_title_suffix))
        ax.set_ylabel('%s %s' % (yname[label_index_level], axis_title_suffix))
    else:
        ax.set_xlabel('%s %s' % (xname, axis_title_suffix))
        ax.set_ylabel('%s %s' % (yname, axis_title_suffix))

def set_lim(ax, limits = None):
    if limits != None:
        if 'x' in limits:
            ax.set_xlim(limits['x'][0], limits['x'][1])
        if 'y' in limits:
            ax.set_ylim(limits['y'][0], limits['y'][1])

def set_ticklabels(ax, fig, ticklabels = None):
    '''
    Swap out default tick labels for custom
    '''
    if ticklabels != None:
        #without calling canvas.draw(), ticklabels may all be set to ''
        fig.canvas.draw()
        if 'xlabel' in ticklabels:
            for tick in ax.get_xticklabels():
                newlabel = ticklabels['xlabel'].get(tick.get_text(), '')
                tick.set_text(newlabel)
            ax.set_xticklabels(ax.get_xticklabels())
        if 'ylabel' in ticklabels:
            for tick in ax.get_yticklabels():
                newlabel = ticklabels['ylabel'].get(tick.get_text(), '')
                tick.set_text(newlabel)
            ax.set_yticklabels(ax.get_yticklabels())

def save(fig, filename = None, title = '', figformat = 'png', extra_artists = None):
    '''
    Save figure
    '''
    #extra_artists = extra_artists[0]
    #plt.tight_layout()
    plt.suptitle(title)
    #plt.subplots_adjust(top = 0.9)
    plt.savefig('%s.%s' % (filename, figformat), bbox_extra_artists = (extra_artists), bbox_inches = 'tight')
    ##plt.savefig('%s.%s' % (filename, figformat))
    ##plt.close(fig)

def filter_df(df, filter_col = None):
    '''
    Remove rows that are not = True in this column
    If filter_col  is actually a list, e.g. [('rep1', 'filter'), ('rep2', 'filter'),...]
    Then this is probably from a multiIndex and test if they all match filter
    '''

    #test if there are one or two levels in df:
    if type(filter_col) == list:
        filter_mask = df.loc[:, filter_col].all(axis = 1)
        df = df[filter_mask].copy()

    #only 1 level of indexing
    else:
        df = df[df[filter_col] == True].copy()

    ##if you don't return the copy made here, then it won't update the one we're working on
    #is there a pandas command to directly modify the copy inplace rather than creating yet another copy?
    return df

def add_text(ax, s, x, y):
    '''
    note: Adding this because doing it after return seems to put in different place
    '''
    ax.text(x, y, s, transform = ax.transAxes)

def plot(df, cols = None, plottype = None, logbase = None, title = '', label_index_level = 0, axis_title_suffix = '', filter_col = None, filename = None, limits = None, figformat = 'png', ticklabels = None, labels = None, **kwargs):
    '''
    Given a df and plottype, clean up data and then send to plotting function
    df = pandas dataframe, cols = list of columns with data, plottype = 'scatter', etc.
    filter_column = name of column whose values need to be True to include in the analysis
    limits = {'x':[-1, 1], 'y':[-1, 1]}, for a 1D plot, only use x
    labels = {'ylabel':ylabel, 'xlabel':xlabel} #will overwrite any inferred labels from the column names
    ticklabels = {'xticks':[tick1, tick2,...], 'yticks':[tick1, tick2, ...]} #will overwrite any existing ticklabels
    '''
    #make a copy of the df here. All subsequent operations will modify this copy
    df = df.copy()

    if filter_col != None:
        df = filter_df(df, filter_col =  filter_col)

    df = df[cols].copy()
    clean_df(df, cols = cols)

    #get new log-transformed columns if required
    if logbase != None:
        cols = log_transform(df, cols = cols, logbase = logbase, label_index_level = label_index_level)
        clean_df(df, cols = cols)

    #send to plotting fxn, which will return plot-specific analyses in the results dict
    results = plot_fxn_dict[plottype](df, cols = cols, title = title, label_index_level = label_index_level, axis_title_suffix = axis_title_suffix, limits = limits, ticklabels = ticklabels, **kwargs)

    #formatting, because this works by axis, this is called separately in the multiaxis functions
    if plottype not in ['multiscatter', 'multiscatter2']:
        if limits != None:
            set_lim(results['ax'], limits = limits)

        if ticklabels != None:
            #without calling canvas.draw(), ticklabels may all be set to ''
            set_ticklabels(results['ax'], results['fig'], ticklabels =  ticklabels)


    #this doesn't work:(, get_xticklabels not defined. Dang... How can you accomplish this then?

    if labels != None:
        if 'xlabel' in labels:
            results['ax'].set_xlabel(labels['xlabel'])
        if 'ylabel' in labels:
            results['ax'].set_ylabel(labels['ylabel'])

    #todo: maybe set an option to print the ticklabels
    #as with different ranges will be really difficult to tell what the  actual labels will be in advance
    if 'text' in kwargs:
        add_text(results['ax'], *kwargs['text'])

    if 'ax_loc' in kwargs:
        loc = plticker.MultipleLocator(base = float(kwargs['ax_loc'])) # this locator puts ticks at regular intervals
        results['ax'].yaxis.set_major_locator(loc)
        results['ax'].xaxis.set_major_locator(loc)

    if 'xy_line' in kwargs:
        lims = [
            np.min([results['ax'].get_xlim(), results['ax'].get_ylim()]),
            np.max([results['ax'].get_xlim(), results['ax'].get_ylim()])
        ]
        # min of both axes  # max of both axes
        results['ax'].plot(lims, lims, 'k-', linestyle = '--')

    #hacky fix for now, but should make a default results_dict with these keys
    if 'extra_artists' not in results:
        results['extra_artists'] = None

    save(results['fig'], filename = filename, title = title, figformat = figformat, extra_artists = results['extra_artists'])
    #avoid cropping legend on file save: https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box
    #results['fig'].savefig(filename, bbox_extra_artists=(l,), bbox_inches='tight')

    #Uncomment for normal use but while debugging, if you close, will not display in Jupyter
    plt.close(results['fig'])

    #print('saved fig')
    return results

plot_fxn_dict = {'multiscatter':multiscatter_plot, 'scatter':scatter_plot, 'box':seaborn_box, 'hist':barhist_plot, 'quickbar':quick_barplot, 'regplot': reg_plot, 'stacked_bar': stacked_bar}
np_log_transforms = {2:  np.log2, 10: np.log10}

#https://stackoverflow.com/questions/41122923/getting-empty-tick-labels-before-showing-a-plot-in-matplotlib
