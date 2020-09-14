#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Capture-C Pipeline
==================

This script processes data from the capture-c or NG capture-c sequencing
protocols designed to identify 3d interactions in the genome.

It takes Illumina paired-end sequencing reads in fastq format (gz compressed)
as input and performs the following steps:

1: Identifies all restriction fragments in the genome (digest_genome.py)
2: Quality control of raw reads (fastqc, multiqc)
3: Removal of PCR duplicates (based on exact sequence matches; deduplicate_fastq.py)
4: Trimming of reads to remove adaptor sequence (trim_galore)
5: Combining overlapping read pairs (flash)
6: Splits fastq(s) into smaller files (enables fast parallel processing; split_fastq.py)
7: In silico digestion of reads (digest_fastq.py)
8: Alignment of fastq files with a user specified aligner (i.e. bowtie/bowtie2; BWA is not supported)
9: Analysis of alignment statistics (picard CollectAlignmentSummaryMetrics, multiqc)
10: Intersection of mapped reads with: capture probes, exclusion regions, blacklist, restriction fragments
11: Classification of mapped read slices (ccanalyser.py)
12: Generation of bedgraphs/BigWigs (convert_tsv_to_bedgraph.py and make_bigwig)
13: Generation of a UCSC track hub for visualisation
14: Collation of run statistics and generation of a run report


@authors: asmith, dsims
"""
# import packages
from ccanalysis.ccanalyser import SliceFilter
import sys
import os
#import gzip
import seaborn as sns
import matplotlib.colors
import pandas as pd
#from pybedtools import BedTool
import itertools
from pysam import FastxFile
from cgatcore import pipeline as P
from ruffus import (
    mkdir,
    follows,
    transform,
    merge,
    originate,
    collate,
    split,
    regex,
    add_inputs,
    suffix,
    active_if,
)


# Global vars
ON_VALS = ['true', 't', 'on', 'yes', 'y', '1']

# Script location
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)

# Read in parameter file
params = 'capturec_pipeline.yml'
P.get_parameters(params)

# Sort pipeline global params
P.PARAMS['SCRIPT_DIR'] = SCRIPT_DIR
P.PARAMS['PACKAGE_DIR'] = PACKAGE_DIR
P.PARAMS['HUB_DIR'] = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
P.PARAMS['ASSEMBLY_DIR'] = os.path.join(P.PARAMS['HUB_DIR'], P.PARAMS['genome_name'])


def main(argv=None):
    if argv is None:
        argv = sys.argv

    P.main(argv)


@follows(mkdir('ccanalysis/restriction_enzyme_map/'))
@transform(
    P.PARAMS['genome_fasta'],
    regex(r'.*/(.*).fa.*'),
    r'ccanalysis/restriction_enzyme_map/genome.digest.bed.gz',
)
def digest_genome(infile, outfile):
    '''Digest genome using restriction enzyme and output fragments in bed file'''

    tmp = outfile.replace('.gz', '')
    restriction_enzyme = (
        f'-r {P.PARAMS["ccanalyser_re"]}' if P.PARAMS["ccanalyser_re"] else ''
    )
    restriction_site = (
        f'-s {P.PARAMS["ccanalyser_re_site"]}' if P.PARAMS["ccanalyser_re_site"] else ''
    )

    statement = '''ccanalyser utils digest_genome
                   -i %(infile)s -o %(tmp)s %(restriction_enzyme)s %(restriction_site)s
                   -l %(tmp)s.log
                   && gzip %(tmp)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir('fastq_pre-processing'), mkdir('fastq_pre-processing/fastqc'))
@transform(
    '*.fastq*', regex(r'(.*).fastq.*'), r'fastq_pre-processing/fastqc/\1_fastqc.zip'
)
def qc_reads(infile, outfile):
    '''Quality control of raw sequencing reads'''
    outdir = os.path.dirname(outfile)
    statement = '''fastqc
                   -q
                   -t %(run_options_threads)s
                   --nogroup %(infile)s
                   --outdir %(outdir)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=P.PARAMS['run_options_threads'],
    )


@follows(mkdir('run_statistics'))
@merge(qc_reads, 'run_statistics/fastqc_report.html')
def multiqc_reads(infile, outfile):
    '''Collate fastqc reports into single report using multiqc'''

    bn = os.path.basename(outfile)
    dn = os.path.dirname(outfile)

    statement = '''rm -f %(outfile)s &&
                   export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc
                   fastq_pre-processing/fastqc/
                   -o %(dn)s
                   -n %(bn)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'], job_memory='2G')


@follows(mkdir('fastq_pre-processing/deduplicated'))
@collate(
    '*.fastq*',
    regex(r'(.*)_R*[12].fastq.*'),
    r'fastq_pre-processing/deduplicated/\1_1.fastq.gz',
)
def deduplicate_reads(infiles, outfile):

    '''Checks for duplicate read1/read2 pairs in a pair of fastq files
       any duplicates are discarded'''

    fq1, fq2 = infiles
    out1, out2 = outfile, outfile.replace('_1.fastq.gz', '_2.fastq.gz')
    stats_file = out1.replace('_1.fastq.gz', '.stats')

    if str(P.PARAMS['deduplication_pre-dedup']).lower() in ON_VALS:

        statement = '''ccanalyser utils deduplicate_fastq
                               -1 %(fq1)s -2 %(fq2)s
                               --out1 %(out1)s --out2 %(out2)s
                               --stats_file %(stats_file)s
                               -c %(run_options_compression_level)s
                              '''
    elif '.fastq.gz' in fq1:
        # If deduplication turned off sylink input files and count number of reads
        statement = '''ln -s $(pwd)/%(fq1)s %(out1)s &&
                       ln -s $(pwd)/%(fq2)s %(out2)s &&
                       lc=$(zcat %(fq1)s | wc -l);
                       logfile=%(stats_file)s;
                       echo -e "Read_pairs_processed\\t$(($lc / 4))" > $logfile;
                       echo -e "Read_pairs_unique\\t$(($lc / 4))" >> $logfile;
                       echo -e "Read_pairs_removed\\t0" >> $logfile'''

    else:
        # If deduplication turned off and not gzipped, count number of reads and gzip fastq
        statement = '''cat %(fq1)s | pigz -p 6 > %(out1)s &
                       cat %(fq2)s | pigz -p 6  > %(out2)s &
                       lc=$(cat %(fq1)s | wc -l);
                       logfile=%(stats_file)s;
                       echo -e "Read_pairs_processed\\t$(($lc / 4))" > $logfile;
                       echo -e "Read_pairs_unique\\t$(($lc / 4))" >> $logfile;
                       echo -e "Read_pairs_removed\\t0" >> $logfile'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'], job_memory='8G')


@follows(mkdir('fastq_pre-processing/trimmed'), deduplicate_reads)
@collate(
    r'fastq_pre-processing/deduplicated/*.fastq.gz',
    regex(r'fastq_pre-processing/deduplicated/(.*)_[12].fastq.gz'),
    r'fastq_pre-processing/trimmed/\1_1_val_1.fq.gz',
)
def trim_reads(infiles, outfile):
    '''Trim adaptor sequences using Trim-galore'''

    fastq1, fastq2 = infiles
    outdir = os.path.dirname(outfile)
    statement = '''trim_galore --cores %(run_options_threads)s --paired %(trim_options)s -o %(outdir)s
                     %(fastq1)s %(fastq2)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=P.PARAMS['run_options_threads'],
    )


@follows(trim_reads, mkdir('fastq_pre-processing/flashed'))
@collate(
    'fastq_pre-processing/trimmed/*.fq.gz',
    regex(r'fastq_pre-processing/trimmed/(.*)_[12]_.*.fq.gz'),
    r'fastq_pre-processing/flashed/\1.extendedFrags.fastq.gz',
)
def combine_reads(infiles, outfile):
    '''Combine overlapping paired-end reads using flash'''
    fastq1, fastq2 = infiles
    output_prefix = outfile.replace('.extendedFrags.fastq.gz', '')
    statement = '''flash
                   -z
                   -t %(run_options_threads)s
                   -o %(output_prefix)s
                    %(fastq1)s
                    %(fastq2)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=P.PARAMS['run_options_threads'],
    )


@follows(mkdir('fastq_pre-processing/split'), combine_reads)
@transform(
    'fastq_pre-processing/flashed/*.fastq.gz',
    regex(r'fastq_pre-processing/flashed/(.*).fastq.gz'),
    r'fastq_pre-processing/split/\1_0.fastq.gz',
)
def split_fastq(infile, outfile):
    '''Splits the combined (flashed) fastq files into chunks for parallel processing'''

    # Small error in function as only processes chunksize - 1 reads
    output_prefix = outfile.replace('_0.fastq.gz', '')
    statement = '''ccanalyser utils split_fastq
                  -i %(infile)s
                  -n %(output_prefix)s
                  --chunksize %(split_n_reads)s
                  -c %(run_options_compression_level)s '''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir('fastq_pre-processing/digested'), split_fastq)
@transform(
    'fastq_pre-processing/split/*.fastq.gz',
    regex(r'fastq_pre-processing/split/(.*).extendedFrags_(\d+).fastq.gz'),
    r'fastq_pre-processing/digested/\1.flashed_\2.fastq.gz',
)
def digest_flashed_reads(infile, outfile):
    '''In silico restriction enzyme digest of combined (flashed) read pairs'''

    restriction_enzyme = (
        f'-r {P.PARAMS["ccanalyser_re"]}' if P.PARAMS["ccanalyser_re"] else ''
    )
    restriction_site = (
        f'-s {P.PARAMS["ccanalyser_re_site"]}' if P.PARAMS["ccanalyser_re_site"] else ''
    )

    statement = '''ccanalyser utils digest_fastq
                   flashed
                   -o %(outfile)s
                   --stats_file %(outfile)s.stats
                   %(restriction_enzyme)s
                   %(restriction_site)s
                   -m 18
                   -c %(run_options_compression_level)s
                   -p 4
                   -i %(infile)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=4,
    )


@follows(split_fastq)
@collate(
    'fastq_pre-processing/split/*.fastq.gz',
    regex(r'fastq_pre-processing/split/(.*).notCombined_[12]_(\d+).fastq.gz'),
    r'fastq_pre-processing/digested/\1.pe_\2.fastq.gz',
)
def digest_pe_reads(infiles, outfile):
    '''In silico restriction enzyme digest of non-combined (non-flashed) read pairs'''

    fq1, fq2 = infiles
    restriction_enzyme = (
        f'-r {P.PARAMS["ccanalyser_re"]}' if P.PARAMS["ccanalyser_re"] else ''
    )
    restriction_site = (
        f'-s {P.PARAMS["ccanalyser_re_site"]}' if P.PARAMS["ccanalyser_re_site"] else ''
    )

    statement = '''ccanalyser utils digest_fastq
                   unflashed
                   --stats_file %(outfile)s.stats
                    %(restriction_enzyme)s
                    %(restriction_site)s
                   -o %(outfile)s
                   -c %(run_options_compression_level)s
                   -m 18
                   -p 4
                   -1 %(fq1)s
                   -2 %(fq2)s
                   '''

    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=4,
    )


@follows(mkdir('aligned'))
@transform(
    [digest_flashed_reads, digest_pe_reads],
    regex(r'fastq_pre-processing/digested/(.*).fastq.gz'),
    r'aligned/\1.bam',
)
def align_reads(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    aligner = P.PARAMS['align_aligner']
    index_flag = P.PARAMS['align_index_flag'] if P.PARAMS['align_index_flag'] else ''
    options = P.PARAMS['align_options'] if P.PARAMS['align_options'] else ''

    statement = '''%(aligner)s %(options)s %(index_flag)s %(align_index)s %(infile)s
                    | samtools view -bS > %(outfile)s
                    && samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(run_options_threads)s
                    && mv %(outfile)s.sorted.bam %(outfile)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_threads=P.PARAMS['run_options_threads'],
        job_memory='4G',
    )


@collate(align_reads, regex(r'aligned/(.*)_(\d+).bam'), r'aligned/\1.bam')
def merge_bam_files(infiles, outfile):
    '''Combines bam files (by flashed/non-flashed status and sample)'''
    fnames = ' '.join(infiles)

    statement = '''samtools merge %(outfile)s %(fnames)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@transform(merge_bam_files, regex('(.*)'), r'\1.bai')
def index_bam(infile, outfile):
    statement = '''samtools index %(infile)s'''
    P.run(statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_memory='1G')


@follows(mkdir('aligned/mapping_statistics'), index_bam)
@transform(
    merge_bam_files,
    regex(r'aligned/(.*).bam'),
    r'aligned/mapping_statistics/\1.picard.metrics',
)
def mapping_qc(infile, outfile):
    '''Uses picard CollectAlignmentSummaryMetrics to get mapping information.'''

    cmd = [
        'picard',
        'CollectAlignmentSummaryMetrics',
        'R=%(genome_fasta)s',
        'I=%(infile)s',
        'O=%(outfile)s',
        '&> %(outfile)s.log',
    ]

    statement = ' '.join(cmd)

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@merge(mapping_qc, 'run_statistics/mapping_report.html')
def mapping_multiqc(infiles, outfile):
    '''Combines mapping metrics using multiqc'''

    indir = os.path.dirname(infiles[0])
    out_fn = os.path.basename(outfile)
    out_dn = os.path.dirname(outfile)
    statement = '''rm -f %(outfile)s &&
                   export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc
                   %(indir)s
                   -o %(out_dn)s
                   -n %(out_fn)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'], job_memory='8G')


@follows(mkdir('ccanalysis/annotations'))
@transform(
    align_reads, regex(r'aligned/(.*).bam'), r'ccanalysis/annotations/\1.bam.bed.gz'
)
def bam_to_bed(infile, outfile):
    '''Converts bam files to bed for faster intersection'''
    tmp = outfile.replace('.gz', '')
    statement = '''bedtools bamtobed
                    -i %(infile)s | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@transform(
    P.PARAMS["ccanalyser_capture"],
    regex(P.PARAMS["ccanalyser_capture"]),
    r'ccanalysis/annotations/exclude.bed.gz',
)
def build_exclusion_bed(infile, outfile):
    '''Generates exclusion window around each capture site'''

    statement = '''bedtools slop
                    -i %(infile)s -g %(genome_fai)s -b %(ccanalyser_exclude_window)s
                    | bedtools subtract -a - -b %(ccanalyser_capture)s
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(digest_genome, build_exclusion_bed)
@transform(
    bam_to_bed,
    regex(r'ccanalysis/annotations/(.*).bam.bed.gz'),
    add_inputs(
        [
            {
                'name': 'restriction_fragment',
                'fn': 'ccanalysis/restriction_enzyme_map/genome.digest.bed.gz',
                'action': 'get',
                'overlap_fraction': 0.2,
            },
            {
                'name': 'capture',
                'fn': P.PARAMS['ccanalyser_capture'],
                'action': 'get',
                'overlap_fraction': 0.9,
            },
            {
                'name': 'exclusion',
                'fn': 'ccanalysis/annotations/exclude.bed.gz',
                'action': 'get',
                'overlap_fraction': 1e-9,
            },
            {
                'name': 'exclusion_count',
                'fn': 'ccanalysis/annotations/exclude.bed.gz',
                'action': 'count',
                'overlap_fraction': 1e-9,
            },
            {
                'name': 'capture_count',
                'fn': P.PARAMS['ccanalyser_capture'],
                'action': 'count',
                'overlap_fraction': 0.9,
            },
            {
                'name': 'blacklist',
                'fn': P.PARAMS['ccanalyser_blacklist'],
                'action': 'count',
                'overlap_fraction': 1e-9,
            },
        ]
    ),
    r"ccanalysis/annotations/\1.annotations.tsv.gz",
)
def annotate_slices(infile, outfile):

    a = infile[0]
    colnames, fnames, actions, fractions = zip(
        *[(d['name'], d['fn'], d['action'], d['overlap_fraction']) for d in infile[1]]
    )

    colnames = ' '.join(colnames)
    fnames = ' '.join([fn if fn else '-' for fn in fnames])
    actions = ' '.join(actions)
    fractions = ' '.join([str(frac) for frac in fractions])

    statement = f'''python /t1-data/user/asmith/Projects/ccanalyser/capture-c/ccanalyser/ccanalysis/generate_annotations.py
                    -a {a} -b {fnames} --actions {actions} --colnames {colnames} -f {fractions} -o {outfile}'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'], job_threads=6)


@follows(
    annotate_slices,
    mkdir('ccanalysis/captures_and_reporters'),
    mkdir('ccanalysis/stats'),
)
@transform(
    align_reads,
    regex(r'aligned/(.*).bam'),
    add_inputs(r'ccanalysis/annotations/\1.annotations.tsv.gz'),
    r'ccanalysis/stats/\1.slice.stats',
)
def ccanalyser(infiles, outfile):
    '''Processes bam files and annotations, filteres slices and outputs
       reporter slices for each capture site'''

    bam, annotations = infiles
    output_prefix = (outfile.replace('/stats/', '/captures_and_reporters/')
                            .replace('.slice.stats', '')
                    )

    stats_prefix = outfile.replace('.slice.stats', '')


    statement = '''ccanalyser ccanalysis ccanalyser
                    -i %(bam)s
                    -a %(annotations)s
                    --output_prefix %(output_prefix)s
                    --stats_out %(stats_prefix)s
                    --method %(ccanalyser_method)s'''
    P.run(
        statement,
        job_queue=P.PARAMS['run_options_queue'],
        job_memory=P.PARAMS['run_options_memory'],
    )



@active_if(P.PARAMS['ccanalyser_method'] in ['capture', 'tri'])
@follows(ccanalyser)
@collate(
    'ccanalysis/captures_and_reporters/*.tsv.gz',
    regex(r'ccanalysis/captures_and_reporters/(.*)\..*_\d+.(.*).tsv.gz'),
    r'ccanalysis/capturec_reporters_aggregated/\1.\2_raw.tsv.gz',
)
def collate_ccanalyser(infiles, outfile):
    '''Combines multiple capture site tsv files'''

    #TODO: Implement a duplication filter after aggregation

    inlist = " ".join(infiles)
    mkdir('ccanalysis/capturec_reporters_aggregated')
    
    statement = '''ccanalyser utils join_tsv
                   -f reporter_read_name
                   -i %(inlist)s
                   -o %(outfile)s
                   --method concatenate'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@transform(collate_ccanalyser,
           regex(r'ccanalysis/capturec_reporters_aggregated/(.*)_raw.tsv.gz'),
           r'ccanalysis/capturec_reporters_aggregated/\1.deduplicated.tsv.gz')
def deduplicate_aggregated_ccananalyser(infile, outfile):

    from ..ccanalysis.ccanalyser import CCSliceFilter

    df = pd.read_csv(infile, sep='\t')
    df_capture = (df.loc[:, df.columns.str.contains('capture|parent')]
                    .rename(columns=lambda col: col.split('_')[1] if 'capture' in col and len(col.split('_')) > 1 else col))
                    
    df_rep = (df.loc[:, df.columns.str.contains('reporter|parent')]
               .rename(columns=lambda col: col.split('_')[1] if 'reporter' in col else col))
    
    df_cat = pd.concat([df_rep, df_capture]).sort_values('parent_read')
    
    sf = CCSliceFilter(df_cat)
    sf.remove_duplicate_slices()
    sf.remove_duplicate_slices_pe()
    sf.merged_captures_and_reporters.to_csv(outfile, sep='\t', index=False)



@active_if(P.PARAMS['ccanalyser_method'] in ['tri', 'tiled'])
@follows(deduplicate_aggregated_ccananalyser)
@transform('ccanalysis/captures_and_reporters_aggregated/*.deduplicated.tsv.gz',
           regex(r'ccanalysis/captures_and_reporters_aggregated/(.*).deduplicated.tsv.gz'),
           r'ccanalysis/restriction_fragment_interaction_counts/\1.tsv.gz')
def count_restriction_fragment_interactions(infile, outfile):
    
    mkdir('ccanalysis/restriction_fragment_interaction_counts')
    statement = f'''python /t1-data/user/asmith/Projects/ccanalyser/capture-c/ccanalyser/ccanalysis/count_restriction_fragment_combinations.py
                    -f {infile} -o {outfile}'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'], job_threads=2)

@follows(mkdir('ccanalysis/bedgraphs'))
@transform(
    collate_ccanalyser_capturec,
    regex(r'ccanalysis/capturec_reporters_aggregated/(.*).tsv.gz'),
    add_inputs(digest_genome),
    r'ccanalysis/bedgraphs/\1.bedgraph.gz',
)
def make_bedgraph(infiles, outfile):
    '''Intersect reporters with genome restriction fragments to create bedgraph'''
    tsv_fn = infiles[0]
    re_map = infiles[1]
    statement = '''ccanalyser ccanalysis ccanalyser_to_bedgraph
                   -i %(tsv_fn)s
                   -b %(re_map)s
                   -o %(outfile)s'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir('ccanalysis/bigwigs'))
@transform(
    make_bedgraph,
    regex(r'ccanalysis/bedgraphs/(.*).bedgraph.gz'),
    r'ccanalysis/bigwigs/\1.bigWig',
)
def make_bigwig(infile, outfile):
    '''Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph'''

    tmp = infile.replace('.gz', '')
    statement = '''zcat %(infile)s > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


def write_dict_to_file(fn, dictionary):
    with open(fn, 'w') as w:
        for k, v in dictionary.items():
            w.write(f'{k} {v}\n')


def get_track_data(fn):

    track_name = ' '.join(os.path.basename(fn).replace('_', ' ').split('.')[:-1])

    track_dict = {
        'track': fn,
        'bigDataUrl': f'{P.PARAMS["hub_url"].rstrip("/")}/{(os.path.join(P.PARAMS["ASSEMBLY_DIR"], fn)).lstrip("/")}',
        'shortLabel': track_name,
        'longLabel': track_name,
        'type': f'{fn.split(".")[-1]}',
    }

    if P.PARAMS['hub_track_options']:
        try:
            options = [op.strip() for op in P.PARAMS['hub_track_options'].split()]
            options_dict = dict(zip(options[0::2], options[1::2]))
            track_dict.update(options_dict)
        except Exception as e:
            print('Invalid custom track options')

    return track_dict


@mkdir(P.PARAMS['HUB_DIR'])
@originate(os.path.join(P.PARAMS['HUB_DIR'], 'hub.txt'))
def generate_hub_metadata(outfile):

    content = {
        'hub': P.PARAMS['hub_name'],
        'shortLabel': P.PARAMS['hub_short']
        if P.PARAMS['hub_short']
        else P.PARAMS['hub_name'],
        'longLabel': P.PARAMS['hub_long']
        if P.PARAMS['hub_long']
        else P.PARAMS['hub_name'],
        'genomesFile': 'genomes.txt',
        'email': P.PARAMS['hub_email'],
        'descriptionUrl': '/'.join(
            [
                P.PARAMS["hub_url"].rstrip('/'),
                P.PARAMS['ASSEMBLY_DIR'].rstrip('/'),
                'visualise_run_statistics.html',
            ]
        ),
    }

    write_dict_to_file(outfile, content)


@follows(generate_hub_metadata)
@originate(os.path.join(P.PARAMS['HUB_DIR'], 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {
        'genome': P.PARAMS['genome_name'],
        'trackDb': os.path.join(P.PARAMS['genome_name'], 'trackDb.txt'),
    }

    write_dict_to_file(outfile, content)


@mkdir(P.PARAMS['ASSEMBLY_DIR'])
@follows(generate_hub_metadata)
@merge(make_bigwig, f'{P.PARAMS["ASSEMBLY_DIR"]}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):

    # Generate all separate tracks
    bigwig_tracks_all = [get_track_data(os.path.basename(fn)) for fn in infiles]

    # Add colours to tracks
    if not P.PARAMS['hub_colors']:
        colors = sns.color_palette('husl', len(bigwig_tracks_all))
        for track, color in zip(bigwig_tracks_all, colors):
            track['color'] = ','.join([str(c * 255) for c in color])
    else:
        for track, color in zip(
            bigwig_tracks_all, itertools.cycle(P.PARAMS['hub_colors'].split(' '))
        ):

            track['color'] = ','.join(
                [str(c * 255) for c in matplotlib.colors.to_rgb(color)]
            )

    # Write track data separated
    with open(outfile, 'w') as w:
        for track in bigwig_tracks_all:
            for label, data in track.items():
                w.write(f'{label} {data}\n')
            # Need to separate each track with a new line
            w.write('\n')

        # Group tracks by sample name and make separate combined tracks for each
        sample_key = lambda d: d['track'].split('.')[0]
        bigwig_tracks_grouped = {
            sample: list(track)
            for sample, track in itertools.groupby(
                sorted(bigwig_tracks_all, key=sample_key), key=sample_key
            )
        }

        for sample, grouped_tracks in bigwig_tracks_grouped.items():

            # Generate overlay track
            combined_track_details = {
                'track': f'{sample}_combined',
                'container': 'multiWig',
                'aggregate': 'transparentOverlay',
                'showSubtrackColorOnUi': 'on',
                'type': 'bigWig 0 250',
                'shortLabel': f'{sample}_combined',
                'longLabel': f'{sample}_all_capture_probes_combined',
                'autoScale': 'on',
                'windowingFunction': 'maximum',
            }

            # Write overlay track
            for label, data in combined_track_details.items():
                w.write(f'{label} {data}\n')
            w.write('\n')

            # Write sub-tracks
            for track in grouped_tracks:
                track['track'] = track['track'].replace('.bigWig', '_subtrack.bigWig')
                for label, data in track.items():
                    w.write(f'\t{label} {data}\n')
                w.write(f'\tparent {combined_track_details["track"]}\n')
                # Need to separate each track with a new line
                w.write('\n')


@follows(ccanalyser)
@merge(
    [
        'fastq_pre-processing/deduplicated/*.stats',
        'fastq_pre-processing/digested/*.stats',
        'ccanalysis/stats/*.slice.stats',
        'ccanalysis/stats/*.reporter.stats',
    ],
    'run_statistics/combined_stats.tsv',
)
def aggregate_stats(infiles, outfile):

    dedup = ' '.join([fn for fn in infiles if 'deduplicated/' in fn])
    digestion = ' '.join([fn for fn in infiles if 'digested/' in fn])
    slices = ' '.join([fn for fn in infiles if '.slice.stats' in fn])
    reporters = ' '.join([fn for fn in infiles if '.reporter.stats' in fn])
    outdir = os.path.dirname(outfile)

    statement = '''ccanalyser stats aggregate_stats
                    --deduplication_stats %(dedup)s
                    --digestion_stats %(digestion)s
                    --ccanalyser_stats %(slices)s
                    --reporter_stats %(reporters)s
                    --output_dir %(outdir)s
                '''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(aggregate_stats)
@merge('run_statistics/*.tsv', r'run_statistics/visualise_run_statistics.html')
def build_report(infile, outfile):
    '''Run jupyter notebook for reporting and plotting. First moves the notebook
       then converts to html'''
    statement = '''rm run_statistics/visualise_run_statistics* -f &&
                   papermill
                   %(PACKAGE_DIR)s/stats/visualise_capture-c_stats.ipynb
                   run_statistics/visualise_run_statistics.ipynb
                   -p directory $(pwd)/run_statistics/ &&
                   jupyter nbconvert
                   --no-input
                   run_statistics/visualise_run_statistics.ipynb
                   run_statistics/visualise_run_statistics.html
                   '''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(generate_trackdb_metadata)
@transform(
    [make_bigwig, build_report], regex('.*/(.*)$'), P.PARAMS['ASSEMBLY_DIR'] + r'/\1'
)
def link_files(infile, outfile):
    try:
        infile_fp = os.path.abspath(infile)
        os.symlink(infile_fp, outfile)
    except Exception as e:
        print(e)


@originate('hub_url.txt')
def write_hub_path(outfile):

    with open(outfile, 'w') as w:
        w.write(
            f'{P.PARAMS["hub_url"].rstrip("/")}/{os.path.join(P.PARAMS["HUB_DIR"], "hub.txt").lstrip("/")}\n'
        )


if __name__ == "__main__":
    main()
