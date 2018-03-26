#171106 Specialized plotting for specific genomic regions and or meta genomic plots
import copy
import numpy as np
import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import HTSeq

#plt.style.use('presentation')
def read_sites_from_moods(pos_file):
    '''
    Given the output file from MOODS motif finder,
    convert to a dictionary of HTSeq genomic intervals.
    '''
    pos_dict = {}

    with open(pos_file, 'r') as f:
        n = 0
        pos_dict = {}
        for line in f:
            fields = line.strip('\n').split(',')
            chrom, start, strand, length = [fields[0], int(fields[2]), fields[3], len(fields[5])]
            pos_dict[n] = HTSeq.GenomicInterval(chrom, start, start + length, strand)
            n += 1
    return pos_dict

def read_sites(pos_file):
    '''
    Given file of positions (tab-separated: chr, 0-based postion, strand), convert to dictionary of HTSeq.GenomicIntervals
    Infer from the first line of the file whether there are transcript names associated with the positions.
    Otherwise, assign a number as the ID
    '''
    with open(pos_file, 'r') as f:
        first_line = f.readline()
    num_fields = len(first_line.split('\t'))
    if num_fields == 4:
        names_included = True
    elif num_fields == 3:
        names_included = False
    else:
        raise NotImplementedError('unsupported input file format')

    with open(pos_file, 'r') as f:
        n = 0
        pos_dict = {}
        for line in f:
            fields = line.strip('\n').split('\t')
            if names_included == True:
                name, chrom, start, strand = [fields[0], fields[1], int(fields[2]), fields[3]]

            else:
                chrom, start, strand = [fields[0], int(fields[1]), fields[2]]
                name = n

            pos_dict[name]=HTSeq.GenomicInterval(chrom, start, start, strand).start_d_as_pos
            n += 1

    return pos_dict

def write_positions(pos_dict, outfile, write_names = False):
    '''
    Given position dict (name->HTSeq.GenomicPosition), write text file of positions
    If you set write_names = True, then the output file will have the names and 4 fields (otherwise only 3 fields)
    '''
    with open(outfile, 'w') as g:
        for i in pos_dict:
            if write_names == True:
                g.write('%s\t%s\t%s\t%s\n' % (i, pos_dict[i].chrom, pos_dict[i].start, pos_dict[i].strand))
            else:
                g.write('%s\t%s\t%s\n' % (pos_dict[i].chrom, pos_dict[i].start, pos_dict[i].strand))

def generate_profile_multi(pos_dict, coverage_files, halfwinwidth = 500, norm_by = 'window'):
    '''Given a dictionary of HTSeq genomic intervals with start positions, plot coverage around those positions.
    pos_dict = {txt:HTSeq.GenomicInterval}, coverage = HTSeq.GenomicArray
    norm_by = how reads will be normalized for each individual start/end site. choose between (downstream, upstream, window).
    Window = the whole window. Downstream/upstream does not include reads at the middle position. In each case, the mean is used.
    _multi() version will average all the coverage vectors from different libraries
    '''
    all_profiles = []
    for c in coverage_files:
        num_genes = 0
        profile = np.zeros(2*halfwinwidth + 1)
        for txt in pos_dict:
            p = pos_dict[txt]
            try:
                window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth + 1, p.strand)
                wincvg = np.fromiter(c[window], dtype='i', count=2*halfwinwidth + 1)

                if norm_by == 'window':
                    denominator = wincvg.mean()

                elif norm_by == 'upstream':
                    denominator = wincvg[:len(wincvg)//2].mean()

                elif norm_by == 'downstream':
                    denominator = wincvg[len(wincvg)//2+1:].mean()

                else:
                    raise NotImplementedError('please choose one of [window, upstream, downstream]')

                #if 0 reads in the window, don't count it
                if denominator == 0.0:
                    continue

                wincvg = wincvg/denominator

                if p.strand == '+':
                    profile+= wincvg

                else:
                    profile+= wincvg[::-1]

                num_genes+=1

            #this will prevent the couple of ones that go off the chromosome from throwing an error.
            except IndexError:
                continue

            #key error can occur if there are no reads mapping to the feature from which the window is derived
            except KeyError:
                continue

        profile = profile/num_genes
        all_profiles.append(profile)

    final_profile = np.mean(np.array(all_profiles), axis = 0)
    return final_profile, num_genes

def write_positions(pos_dict, outfile, write_names = False):
    '''
    Given position dict (name->HTSeq.GenomicPosition), write text file of positions
    If you set write_names = True, then the output file will have the names and 4 fields (otherwise only 3 fields)
    '''
    with open(outfile, 'w') as g:
        for i in pos_dict:
            if write_names == True:
                g.write('%s\t%s\t%s\t%s\n' % (i, pos_dict[i].chrom, pos_dict[i].start, pos_dict[i].strand))
            else:
                g.write('%s\t%s\t%s\n' % (pos_dict[i].chrom, pos_dict[i].start, pos_dict[i].strand))

def generate_profile(pos_dict, coverage, halfwinwidth = 500, norm_by = 'window'):
    '''Given a dictionary of HTSeq genomic intervals with start positions, plot coverage around those positions.
    pos_dict = {txt:HTSeq.GenomicInterval}, coverage = HTSeq.GenomicArray
    norm_by = how reads will be normalized for each individual start/end site. choose between (downstream, upstream, window).
    Window = the whole window. Downstream/upstream does not include reads at the middle position. In each case, the mean is used.
    '''

    num_genes = 0
    profile = np.zeros(2*halfwinwidth + 1)
    for txt in pos_dict:
        p = pos_dict[txt]
        try:
            window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth + 1, p.strand)
            wincvg = np.fromiter(coverage[window], dtype='i', count=2*halfwinwidth + 1)

            if norm_by == 'window':
                denominator = wincvg.mean()

            elif norm_by == 'upstream':
                denominator = wincvg[:len(wincvg)//2].mean()

            elif norm_by == 'downstream':
                denominator = wincvg[len(wincvg)//2+1:].mean()

            else:
                raise NotImplementedError('please choose one of [window, upstream, downstream]')

            #if 0 reads in the window, don't count it
            if denominator == 0.0:
                continue

            wincvg = wincvg/denominator

            if p.strand == '+':
                profile+= wincvg

            else:
                profile+= wincvg[::-1]

            num_genes+=1

        #this will prevent the couple of ones that go off the chromosome from throwing an error.
        except IndexError:
            continue

        #key error can occur if there are no reads mapping to the feature from which the window is derived
        except KeyError:
            continue

    profile = profile/num_genes
    return profile, num_genes

def plot_profiles(profile_dict, plot_order = None, xlabel = 'relative position (nt)', ylabel = 'mean normalized reads', replacement_labels = None, legend_below = True):
    '''
    Plot read profiles in the profile_dict
    profile_dict = {sample_name:{'profile':profile, 'num_genes':num_genes}}
    Note that all profiles have to be the same length window
    replacement_labels = {really_long_annoying_genelist_label:new label}
    plot_order = really_long_annoying_genelist_label in a list order for which ones to plot firs
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    handles = []
    names = []
    if plot_order == None:
        plot_order = list(profile_dict.keys())
    for sample in plot_order:
        halfwinwidth = len(profile_dict[sample]['profile'])//2
        this_range = np.arange(-halfwinwidth, halfwinwidth + 1)
        d, = ax.plot(np.arange(-halfwinwidth, halfwinwidth + 1), profile_dict[sample]['profile'])
        handles.append(d)
        if replacement_labels == None:
            names.append('%s (n = %s )' % (sample, profile_dict[sample]['num_genes']))
        else:
            names.append('%s (n = %s )' % (replacement_labels[sample], profile_dict[sample]['num_genes']))

    #how to get the legend below the axis:
    #https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    #ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), fancybox = True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #hard code limits for now
    ax.set_xlim(-500, 500)
    plt.legend(handles, names, bbox_to_anchor = (0.55, 1))
    plt.tight_layout()

    #plt.legend(handles, names)
    #plt.legend(handles, names, loc = 'upper center', bbox_to_anchor = (0.5, -0.05), fancybox = True, mode = 'expand')

    return ax

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
