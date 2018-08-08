
#1) Imports, this relies on utils keeping same relative path
import sys
util_dir = '../../common_scripts/pipe_utils/'
sys.path.append(util_dir)
from import_file import *

#version1 will take the sam file output by bowtie2 with mapping to rRNA transcripts
#make:
#1) plots containing rRNA contam percentage of different rRNA loci
#2) table summary file with counts and percentages

def plot_rRNA_percent(samfiles, txt_df = None, outdir = None):
    '''
    Count and summarize the rRNA contaminants from library
    input = sam file, e.g. produced by bowtie2 mapping to rRNA transcripts
    '''

    sample_dfs = []
    for i in samfiles:
        print('i', i)
        #name = os.path.basename(i).split('-')[-1]
        name = os.path.basename(i).split('.')[0]
        with open(i, 'r') as f:
            mapping_dict = defaultdict(int)
            #initialize on first record
            thisref = ''
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip('\n').split('\t')
                refname, chrom = fields[0], fields[2]
                if thisref != refname: #we're at a new pair
                    #add the old pair
                    mapping_dict[chrom] += 1
                    #re-assign:
                    thisref = refname

        #get total number mapping reads from sam file:
        total = sum(mapping_dict.values())
        print('total mapping', total)
        del(mapping_dict['*'])

        count_df = pd.DataFrame(list(mapping_dict.items()), columns = ['txt', 'counts'])
        count_df.set_index('txt', drop = True, inplace = True)
        contam_df = pd.merge(count_df, txt_df, how = 'left', left_index = True, right_index = True)
        rRNA_counts = contam_df.groupby(['locus'])['counts'].sum()
        rRNA_df = rRNA_counts.to_frame()
        rRNA_df['rRNA'] = rRNA_df.index

        outname = os.path.join(outdir, '%s_rRNA_contam' % name)
        #if there's not rRNA contam, df could be empty
        if not rRNA_df.empty:
            #new_df = pipeline_aux.quick_barplot(data = rRNA_df, y_col = 'counts', x_col = 'rRNA', divide_by = total, percent = True, outname = outname, x_label = 'transcript', y_label = '% of reads')
            new_df = pipeline_aux.quick_barplot(rRNA_df, cols = ['rRNA', 'counts'], divide_by = total, percent = True, outname = outname, x_label = 'transcript', y_label = '% of reads')['df']
            df2 = new_df[['percent']].copy()
            df2.rename(columns = {'percent':'percent_%s' % name}, inplace = True)
        else:
            rRNA_df['percent_%s' % name] = 0
            df2 = rRNA_df['percent_%s' % name].copy()
        sample_dfs.append(df2)

    #print('sample_dfs', sample_dfs)
    #make large df and write as summary outfile
    summary_df =  pd.concat(sample_dfs, axis = 1)
    summary_df.loc['total_contam'] = summary_df.sum()
    summary_df.to_csv(os.path.join(outdir, 'rRNA_contamination_summary.csv'))

def plot_contam_positions(samfiles, txt_df = None, outdir = None):
    '''
    Go through the rRNA contam. mapping sam file
    and plot the start position of the reads along the transcripts
    This is not currently very correct because it works by binning since there
    are multiple versions of each transcript and some of them are very different lengths
    '''

    sample_dfs = []
    for i in samfiles:
        print('i', i)
        #name = os.path.basename(i).split('-')[-1]
        name = os.path.basename(i).split('.')[0]
        aln_file = HTSeq.SAM_Reader(i)

        nreads = 0
        big_dict = defaultdict(list)

        #add the positions the hits occur at to the list of transcripts
        #because we mapped to transcripts, the chrom should be the transcript
        for aln in aln_file:
            if aln.aligned:
                #print(aln)
                big_dict[aln.iv.chrom].append(aln.iv.start)
                nreads += 1

        txt_2_family = txt_df['locus'].to_dict()

        #go through binned transcript and add to transcript family vector
        count_by_pos_bin = {}
        binned_pos = {}
        txt_family_arrays = defaultdict(list)

        for txt in big_dict:
            big_dict[txt] = np.array(big_dict[txt])
            hist = np.histogram(big_dict[txt], bins = 100, density = False)
            count_by_pos_bin[txt] = hist[0]
            binned_pos[txt_2_family[txt]] = hist[1]
            txt_family_arrays[txt_2_family[txt]].append(hist[0])

        #get sum of txt_family arrays
        summed_arrays = {}
        for txt in txt_family_arrays:
            arrays = np.vstack(txt_family_arrays[txt])
            summed = np.sum(arrays, axis = 0)
            summed_arrays[txt] = summed

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(binned_pos[txt][:-1], summed)
            plt.title('%s contamination' % txt)
            ax.set_ylabel('read counts')
            ax.set_xlabel('position')
            plt.savefig(os.path.join(outdir, '%s_%s_contam.png' % (name, txt)))
            plt.close(fig) #CHECK IF CORRECT SYNTAX


def plot_contam_positions_by_txt(samfiles, txt_df = None, outdir = None, plot_rRNA_gene_positions = False):
    '''
    Go through the rRNA contam. mapping sam file
    and plot the start position of the reads along the transcripts.
    This function will specifically plot reads mapping to each transcript separately,
    rather than pooling transcripts from the same locus as this doesn't match up for some of them well.
    We need to know the length of the rRNA transcript in order to make the array
    '''

    #overlay specific gene regions on top of the transcript
    #genes_to_plot_dict = {}

    txt_2_family = txt_df['locus'].to_dict()

    ##NOTE: THIS IS CURRENTLY HARDCODED FOR THE FB 6.18 VERSION
    prerRNA_dict = {'FBtr0346881': {'prerRNA':[66807, 74924], '18S': [67668, 69662], '5.8S': [70389, 70511], '2S': [70540, 70569], '28S': [70955, 74924]},\
    'FBtr0346873': {'prerRNA': [42623, 49485], '18S': [43484, 45478], '5.8S': [46205, 46327], '2S': [46356, 46385], '28S': [46771, 49485]},\
    'FBtr0346877': {'prerRNA': [55104, 59294], '18S': [55965, 57959], '5.8S': [58654, 58776], '2S': [58805, 58834]}}

    #adjust each to the coordinates of the pre-RNA
    adj_pre_dict = {}

    for pre in prerRNA_dict:
        adj_pre_dict[pre] = {}
        start, end = prerRNA_dict[pre]['prerRNA'][0:]
        for txt in prerRNA_dict[pre]:
            if txt == 'prerRNA':
                continue
            else:
                adj_pre_dict[pre][txt] = [prerRNA_dict[pre][txt][0] - start, prerRNA_dict[pre][txt][1] - start]
    #'FBtr0346881' #the largest preRNA rDNA:66807..74924
    #'FBtr0346882' #the 18S associated rDNA:67668..69662
    #chrom, start, strand, length = [fields[0], int(fields[2]), fields[3], len(fields[5])]
    #pos_dict[n] = HTSeq.GenomicInterval(chrom, start, start + length, strand)

    for i in samfiles:
        length_dict = {}
        this_sam = HTSeq.FileOrSequence(i)
        for line in this_sam:
            if line.startswith('@SQ'):
                fields = line.strip('\n').split('\t')
                length_dict[fields[1].split(':')[1]] = int(fields[2].split(':')[1])
            elif line.startswith('@PG'):
                break

        name = os.path.basename(i).split('.')[0]

        cvg = genome_functions.get_coverage_vector(i)

        for txt in length_dict:
            if txt in cvg.chrom_vectors:
                start = 0
                end = length_dict[txt] - 1
                ax = genome_functions.plot_genomic_region(cvg, txt, start, end, '+')
                if plot_rRNA_gene_positions == True:
                    if txt in adj_pre_dict:
                        for subtxt in adj_pre_dict[txt]:
                            ax.axvline(x = adj_pre_dict[txt][subtxt][0], color = 'r', linestyle = '--', linewidth = 0.2)
                            ax.axvline(x = adj_pre_dict[txt][subtxt][1], color = 'r', linestyle = '--', linewidth = 0.2)

                plt.savefig(os.path.join(outdir, '%s_%s_%s_contam.png' % (name, txt_2_family[txt], txt)))


def main(arglist):

    if not 'snakemake' in globals():
        parser = argparse.ArgumentParser()
        parser.add_argument('-rRNA_txt_table', help = 'file containing the names of the rRNA transcripts')
        parser.add_argument('-rRNA_mapping_files', nargs = '+', help = 'list of sam files from reads mapped to rRNA transcripts')
        parser.add_argument('-rRNA_outdir', help = 'directory to output the plots and tables with rRNA contamination analysis')
        parser.add_argument('--plot_rRNA_gene_positions', action = 'store_true', help = 'plot the positions of the genes on the pre-rRNA transcripts')
        args = parser.parse_args(args = arglist)
        rRNA_txt_table, rRNA_mapping_files, rRNA_outdir, plot_rRNA_gene_positions = args.rRNA_txt_table, args.rRNA_mapping_files, args.rRNA_outdir, args.plot_rRNA_gene_positions

    else:
        rRNA_txt_table = snakemake.input['rRNA_txt_table']
        rRNA_mapping_files = snakemake.input['rRNA_mapping_files']
        rRNA_outdir = snakemake.params['rRNA_outdir']
        plot_rRNA_gene_positions = False

    os.makedirs(rRNA_outdir, exist_ok = True)
    txt_df = pd.read_csv(rRNA_txt_table, index_col = 'txt')
    plot_rRNA_percent(rRNA_mapping_files, txt_df = txt_df, outdir = rRNA_outdir)
    #plot_contam_positions(rRNA_mapping_files, txt_df = txt_df, outdir = rRNA_outdir)
    plot_contam_positions_by_txt(rRNA_mapping_files, txt_df = txt_df, outdir = rRNA_outdir, plot_rRNA_gene_positions = plot_rRNA_gene_positions)

'''
#local test:
Marys-MacBook-Pro-5:analysis_scripts maryk.thompson$ python plot_rRNA_contam.py -rRNA_txt_table ~/Desktop/Davislab/comp_labbook_backup/genomes/dmel_r6.18_FB/rRNA_txt_table.csv -rRNA_mapping_files ~/Desktop/Davislab/comp_labbook_backup/data/computational_projects/C2/C2.15.SnakeMakeSetup/pipeline_testing/byTxt/plusI/intermediates/rRNA/input1_head.sam ~/Desktop/Davislab/comp_labbook_backup/data/computational_projects/C2/C2.15.SnakeMakeSetup/pipeline_testing/byTxt/plusI/intermediates/rRNA/pd1_head.sam -rRNA_outdir ~/Desktop/Davislab/2.25_rRNA_sub_2/rRNA_contam_spliceseq/rRNA_contam9
'''
if __name__ == '__main__':
    main(sys.argv[1:])
