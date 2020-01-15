#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Capture-C Pipeline
==================

This script processes data from the capture-c or NG capture-c sequencing 
protocols designed to identify 3d interactions in the genome. 

It takes Illumina paired-end sequencing reads in fastq format (gz compressed) 
as input and performs the following steps:

1: Quality control of raw reads (fastqc, multiqc)
2: Trimming of reads to remove adaptor sequence (trim-galore)
3: Combining overlapping read pairs (flash)
4: In silico digestion of reads (digest_fastq.py)
5: Classification of mapped read slices (ccanalyser.py)
6: Deduplication of fragment (dedup_fragments.py)

@authors: asmith, dsims
"""

# import packages
import sys
import os
import gzip
import seaborn as sns
from pysam import FastxFile
from cgatcore import pipeline as P
from ruffus import mkdir, follows, transform, merge, originate, collate, split, regex, add_inputs, suffix

# Read in parameter file
P.get_parameters('capturec_pipeline.yml')
hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])


@follows(mkdir('ccanalyser'), mkdir('ccanalyser/restriction_enzyme_map/'))
@transform(P.PARAMS['genome_fasta'], 
           regex(r'.*/(.*).fa.*'), 
           r'ccanalyser/restriction_enzyme_map/\1.digest.bed.gz')
def digest_genome(infile, outfile):
    '''Digest genome using restriction enzyme and output fragments in bed file'''
    tmp = outfile.replace('.gz', '')
    statement = '''python %(scripts_dir)s/digest_genome.py 
                   -i %(infile)s -o %(tmp)s -r %(ccanalyser_re)s
                   -l %(tmp)s.log
                   && gzip %(tmp)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'])


@follows(mkdir('fastqc'))
@transform('*.fastq.gz', 
           regex(r'(.*).fastq.gz'), 
           r'fastqc/\1_fastqc.zip')
def qc_reads(infile, outfile):
    '''Quality control of raw sequencing reads'''
    statement = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])


@follows(mkdir('report'))
@merge(qc_reads, 'report/readqc_report.html')
def multiqc_reads (infile, outfile):
    '''Collate fastqc reports into single report using multiqc'''
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -o report -n readqc_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G')


@follows(mkdir('deduplicated'))
@collate('*.fastq.gz', 
         regex(r'(.*)_[12].fastq.gz'), 
         r'deduplicated/\1_1.fastq.gz')
def deduplicate_reads(infiles, outfile):

    fq1, fq2 = infiles
    out1, out2 = outfile, outfile.replace('_1.fastq.gz', '_2.fastq.gz')
    logfile = fq1.replace('_1.fastq.gz', '.log')
    statement = '''python %(scripts_dir)s/deduplicate_fastq.py
                           -1 %(fq1)s -2 %(fq2)s
                           --out1 %(out1)s --out2 %(out2)s
                           -l %(outfile)s.log
                           -c %(compression)s
                          '''
    
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='32G')

@follows(mkdir('trim'), deduplicate_reads)
@collate(r'deduplicated/*.fastq.gz',
         regex(r'deduplicated/(.*)_[12].fastq.gz'), 
         r'trim/\1_1_val_1.fq.gz')
def trim_reads(infiles, outfile):
    '''Trim adaptor sequences using Trim-galore'''
    fastq1, fastq2 = infiles
    statement = '''trim_galore --cores %(threads)s %(trim_options)s -o trim 
                     %(fastq1)s %(fastq2)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])
    

@follows(trim_reads, mkdir('flash'))
@collate('trim/*.fq.gz', 
         regex(r'trim/(.*)_[12]_.*.fq.gz'), 
         r'flash/\1.extendedFrags.fastq.gz')
def combine_reads(infiles, outfile):
    '''Combine overlapping paired-end reads using flash'''
    fastq1, fastq2 = infiles
    output_prefix = outfile.replace('.extendedFrags.fastq.gz', '')
    statement = '''flash -z -t %(threads)s 
                    -o %(output_prefix)s 
                    %(fastq1)s %(fastq2)s'''
    P.run(statement, 
      job_queue=P.PARAMS['queue'], 
      job_threads=P.PARAMS['threads'])


@follows(mkdir('split_fastq'), combine_reads)
@transform('flash/*.fastq.gz', regex(r'flash/(.*).fastq.gz'), r'split_fastq/\1_0.fastq.gz')
def split_fastq(infile, outfile):
    '''Splits flashed fastq files into chunks for parallel processing'''
    
    #Small error in function as only processes chunksize - 1 reads
    output_prefix = outfile.replace('_0.fastq.gz', '')
    statement = '''python %(scripts_dir)s/split_fastq.py
                  -i %(infile)s -o %(output_prefix)s
                  --chunk_size %(chunksize)s -c %(compression)s '''
    P.run(statement, 
          job_queue=P.PARAMS['queue'])


@follows(mkdir('digest'), split_fastq)          
@transform('split_fastq/*.fastq.gz', 
                   regex(r'split_fastq/(.*).extendedFrags_(\d+).fastq.gz'), 
                    r'digest/\1.flashed_\2.fastq.gz')
def digest_flashed_reads(infile, outfile):
    '''In silico restriction enzyme digest of combined (flashed) readpairs'''
    statement = '''python %(scripts_dir)s/digest_fastq.py 
                   -o %(outfile)s 
                   -l %(outfile)s.log
                   -r %(ccanalyser_re)s
                   -m 18
                   -c %(compression)s
                   flashed
                   -i %(infile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'])

@follows(split_fastq)
@collate('split_fastq/*.fastq.gz', 
               regex(r'split_fastq/(.*).notCombined_[12]_(\d+).fastq.gz'), 
               r'digest/\1.pe_\2.fastq.gz')
def digest_pe_reads(infiles, outfile):
    '''In silico restriction enzyme digest of non-combined (non-flashed) readpairs'''
    
    fq1, fq2 = infiles
    
    statement = '''python %(scripts_dir)s/digest_fastq.py 
                   -l %(outfile)s.log 
                   -r %(ccanalyser_re)s
                   -o %(outfile)s
                   -c %(compression)s
                   -m 18
                   unflashed
                   -1 %(fq1)s 
                   -2 %(fq2)s  
                   '''
    
    P.run(statement, 
          job_queue=P.PARAMS['queue'])


@follows(mkdir('bam'))
@transform([digest_flashed_reads, digest_pe_reads], 
           regex(r'digest/(.*).fastq.gz'), 
           r'bam/\1.bam')
def align_reads(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    
    options = P.PARAMS['bowtie2_options'] if P.PARAMS['bowtie2_options'] else ''
        
    statement = '''bowtie2 -x %(bowtie2_index)s -U %(infile)s 
                    -p %(threads)s %(options)s 
                    | samtools view -bS > %(outfile)s 2> %(outfile)s.log
                    && samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(threads)s
                    && mv %(outfile)s.sorted.bam %(outfile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_memory='20G')

@collate(align_reads, 
         regex(r'bam/(.*)_\d+.bam'), 
         r'bam/\1.bam')
def merge_bam_files(infiles, outfile):
    '''Combines bam files (by flashed/non-flashed status and sample)'''
    fnames = ' '.join(infiles)
    
    statement = '''samtools merge %(outfile)s %(fnames)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'])
    
    
@transform(merge_bam_files, 
           regex(r'bam/(.*).bam'), 
           r'\1.picard.metrics')
def mapping_qc(infile, outfile):
    
    cmd = ['picard CollectAlignmentSummaryMetrics',
           '-R %(genome_fasta)s -I %(infile)s -O %(outfile)s',
           '&> %(outfile)s.log',
           ]
    
    statement = ' '.join(cmd)
         
    P.run(statement, 
          job_queue=P.PARAMS['queue'])


@merge(mapping_qc, 'report/mapping_report.html')
def mapping_multiqc (infile, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc bam/ -o report -n mapping_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G')


@follows(mkdir('ccanalyser/annotations'))
@transform(align_reads, 
           regex(r'bam/(.*).bam'), 
           r'ccanalyser/annotations/\1.bam.bed.gz')
def bam_to_bed(infile, outfile):
    tmp = outfile.replace('.gz', '')
    statement =  '''bedtools bamtobed 
                    -i %(infile)s | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(P.PARAMS["ccanalyser_capture"], 
           regex(P.PARAMS["ccanalyser_capture"]), 
           r'ccanalyser/annotations/exclude.bed.gz')
def build_exclusion_bed(infile, outfile):
    statement =  '''bedtools slop  
                    -i %(infile)s -g %(genome_fai)s -b %(ccanalyser_exclude_window)s
                    | bedtools subtract -a - -b %(ccanalyser_capture)s
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])

@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'), 
           r'ccanalyser/annotations/\1.annotation.capture.count.gz')
def capture_intersect_count(infile, outfile):
    '''Intersect reads with capture and exclusion files.
    report count of overlaps for each input bed using -C '''
    
    statement =  '''bedtools intersect -c -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\tcapture_count'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'), 
           r'ccanalyser/annotations/\1.annotation.exclude.count.gz')
def exclusion_intersect_count(infile, outfile):
    '''Intersect reads with capture and exclusion files.
    report count of overlaps for each input bed using -C '''
    
    statement =  '''bedtools intersect -c
                    -a %(infile)s -b ccanalyser/annotations/exclude.bed.gz
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\texclusion_count'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'), 
           r'ccanalyser/annotations/\1.annotation.capture.gz')
def capture_intersect(infile, outfile):
    statement =  '''bedtools intersect -loj -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\tcapture'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'), 
           r'ccanalyser/annotations/\1.annotation.exclude.gz')
def exclusion_intersect(infile, outfile):
    statement =  '''bedtools intersect -loj
                    -a %(infile)s -b ccanalyser/annotations/exclude.bed.gz
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\texclusion'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           add_inputs(digest_genome),
           r'ccanalyser/annotations/\1.annotation.re.gz')
def re_intersect(infiles, outfile):  
    bam, genome = infiles        
    statement =  '''bedtools intersect -loj
                    -a %(bam)s -b %(genome)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\trestriction_fragment'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])

@transform(bam_to_bed, 
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'), 
           r'ccanalyser/annotations/\1.annotation.blacklist.count.gz')
def blacklist_intersect_count(infile, outfile):
    '''Intersect reads with blacklisted regions.
    report count of overlaps for each input bed using -C '''
    statement =  '''bedtools intersect -c 
                    -a %(infile)s -b %(ccanalyser_blacklist)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\tblacklist'
                    | gzip > %(outfile)s''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@collate([capture_intersect_count, exclusion_intersect_count,
         capture_intersect, exclusion_intersect,
         re_intersect, blacklist_intersect_count], 
    regex(r'ccanalyser/annotations/(.*).annotation.*'),
    r"ccanalyser/annotations/\1.annotations.tsv.gz")
def merge_annotations(infiles, outfile):
    '''merge all intersections into a single file '''
    inlist = " ".join(infiles)
    statement = '''python %(scripts_dir)s/join_tsv.py  
                        -f read_name
                        -o %(outfile)s
                        -i %(inlist)s'''
    P.run(statement,
          job_queue=P.PARAMS['queue'])


@follows(merge_annotations, mkdir('ccanalyser/bed_files'), mkdir('ccanalyser/stats'))
@transform(align_reads, 
           regex(r'bam/(.*).bam'), 
           add_inputs(r'ccanalyser/annotations/\1.annotations.tsv.gz'),
           r'ccanalyser/bed_files/\1.reporter.bed.gz')
def ccanalyser(infiles, outfile):
    bam, annotations = infiles
    bed_out = outfile.replace('.reporter.bed.gz', '') 
    stats_out = bed_out.replace('bed_files', 'stats')
    statement =  '''python %(scripts_dir)s/ccanalyser.py 
                    -i %(bam)s
                    -a %(annotations)s
                    --bed_output %(bed_out)s 
                    --stats_out %(stats_out)s''' 
    P.run(statement,
          job_queue=P.PARAMS['queue'], 
          job_memory=P.PARAMS['ccanalyser_memory'])

@follows(ccanalyser, mkdir('ccanalyser/bed_files_combined'))
@collate('ccanalyser/bed_files/*.bed.gz',
         regex(r'ccanalyser/bed_files/(.*).digest_.*_\d+_(.*).bed.gz'),
         r'ccanalyser/bed_files_combined/\1_\2.bed.gz')
def collate_ccanalyser_output(infiles, outfile):
    infiles = ' '.join(infiles)
    statement = '''cat %(infiles)s > %(outfile)s'''
    P.run(statement,
          job_queue=P.PARAMS['queue'])

@follows(mkdir('ccanalyser/bedgraphs'))
@transform(collate_ccanalyser_output, 
           regex(r'ccanalyser/bed_files_combined/(\w+_.*).bed.gz'),
           add_inputs(digest_genome),
           r'ccanalyser/bedgraphs/\1.bedgraph.gz')

def make_bedgraph(infiles, outfile):
    '''Intersect reporters with genome restriction fragments to create bedgraph'''
    bed_fn = infiles[0]
    re_map = infiles[1]
    statement = '''bedtools coverage -counts -a %(re_map)s -b %(bed_fn)s 
                   | awk 'BEGIN{OFS="\\t"}{print $1, $2, $3, $5}'
                   | sort -k1,1 -k2,2n
                   | gzip > %(outfile)s'''
    P.run(statement,
          job_queue=P.PARAMS['queue'])


@follows(mkdir('ccanalyser/bigwigs'))
@transform(make_bedgraph,
           regex(r'ccanalyser/bedgraphs/(.*).bedgraph.gz'),
           r'ccanalyser/bigwigs/\1.bigWig')
def make_bigwig(infile, outfile):
    tmp = infile.replace('.gz', '')
    statement = '''zcat %(infile)s > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s'''
    
    P.run(statement,
          job_queue=P.PARAMS['queue'])


@follows(mkdir(hub_dir))
@originate(os.path.join(hub_dir, 'hub.txt'))
def generate_hub_metadata(outfile):

    content = {'hub': P.PARAMS['hub_name'],
               'shortLabel': P.PARAMS['hub_short'] if P.PARAMS['hub_short'] else P.PARAMS['hub_name'],
               'longLabel': P.PARAMS['hub_long'] if P.PARAMS['hub_long'] else P.PARAMS['hub_name'],
               'genomesFile': 'genomes.txt',
               'email': P.PARAMS['hub_email'],
               'descriptionUrl': f'http://userweb.molbiol.ox.ac.uk/{P.PARAMS["hub_publoc"].strip("/")}',
               }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')
            

@follows(generate_hub_metadata)
@originate(os.path.join(hub_dir, 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {'genome': P.PARAMS['hub_genome'],
                      'trackDb': os.path.join(P.PARAMS['hub_genome'], 'trackDb.txt'),
                     }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')


@follows(generate_hub_metadata, mkdir(assembly_dir))
@merge(make_bigwig, f'{assembly_dir}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):
    def get_track_data(fn):
        return {'track': fn,
                'bigDataUrl': f'http://userweb.molbiol.ox.ac.uk/{(os.path.join(assembly_dir, fn)).lstrip("/")}',
                'shortLabel': fn,
                'longLabel': fn,
                'type': f'{fn.split(".")[-1]}',
                }
    
    # Generate all separate tracks
    tracks = [get_track_data(os.path.basename(fn)) for fn in infiles]
    
    # Add colours to tracks
    colors = sns.color_palette('husl', len(tracks))
    for track, color in zip(tracks, colors):
        track['color'] = ','.join([str(c * 255) for c in color])
    
    
    # Write track data separated
    with open(outfile, 'w') as w:
        for track in tracks:
            for label, data in track.items():
                w.write(f'{label} {data}\n')
            # Need to separate each track with a new line
            w.write('\n')
        
        
        # Generate overlay track
        combined_track_details = {'track': f'{P.PARAMS["hub_name"]}_combined',
                                                   'container': 'multiWig',
                                                   'aggregate': 'transparentOverlay',
                                                   'showSubtrackColorOnUi': 'on',
                                                   'type': 'bigWig 0 1000',
                                                   'shortLabel': f'{P.PARAMS["hub_name"]}_combined',
                                                   'longLabel': f'{P.PARAMS["hub_name"]}_combined',}
        
        # Write overlay track
        for label, data in combined_track_details.items():
            w.write(f'{label} {data}\n')
        w.write('\n')
        
        # Write sub-tracks
        for track in tracks:
            track['track'] = track['track'].replace('.bigWig', '_subtrack.bigWig')
            for label, data in track.items():
                w.write(f'\t{label} {data}\n')
            w.write(f'\tparent {combined_track_details["track"]}\n')
            # Need to separate each track with a new line
            w.write('\n') 
            
            
@follows(generate_trackdb_metadata)
@transform(make_bigwig, regex('ccanalyser/bigwigs/(.*).bigWig'),  f'{assembly_dir}' + r'/\1.bigWig')
def link_bigwigs(infile, outfile):
    
    infile_fp = os.path.abspath(infile)
    os.symlink(infile_fp, outfile)



    
    

#@follows(ccanalyser)
#@transform('ccanalyser/*.stats', regex(r'ccanalyser/(.*).stats'), r'report/\1.html')
#def build_report(infile, outfile):
#    '''Run jupyter notebook for reporting and plotting'''

if __name__ == "__main__":
    sys.exit( P.main(sys.argv))
