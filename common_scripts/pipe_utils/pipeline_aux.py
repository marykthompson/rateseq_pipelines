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
    plt.savefig('%s.png' % outname)
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
