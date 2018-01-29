
#1) Imports, this relies on utils keeping same relative path
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
            new_df = pipeline_aux.quick_barplot(data = rRNA_df, y_col = 'counts', x_col = 'rRNA', divide_by = total, percent = True, outname = outname, x_label = 'transcript', y_label = '% of reads')
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


def main():
    txt_file = snakemake.input['rRNA_txt_table']
    txt_df = pd.read_csv(txt_file, index_col = 'txt')
    plot_rRNA_percent(snakemake.input['rRNA_mapping_files'], txt_df = txt_df, outdir = snakemake.params['rRNA_outdir'])
    plot_contam_positions(snakemake.input['rRNA_mapping_files'], txt_df = txt_df, outdir = snakemake.params['rRNA_outdir'])
if __name__ == '__main__':
    main()
