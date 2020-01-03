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
from pysam import FastxFile
from cgatcore import pipeline as P
from ruffus import mkdir, follows, transform, merge, originate, collate, split, regex, add_inputs, suffix

# Read in parameter file
P.get_parameters('capturec_pipeline.yml')

@follows(mkdir('ccanalyser'))
@transform(P.PARAMS['genome_fasta'], 
           regex(r'.*/(.*).fa.*'), 
           r'ccanalyser/\1.digest.bed')
def digest_genome(infile, outfile):
    '''Digest genome using restriction enzyme and output fragments in bed file'''
    statement = '''python %(scripts_dir)s/digest_genome.py 
                   -i %(infile)s -o %(outfile)s -r %(ccanalyser_re)s
                   -l %(outfile)s.log'''
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

    statement = '''python %(scripts_dir)s/deduplicate_fastq.py
                           -1 %(fq1)s -2 %(fq2)s
                           --out1 %(out1)s --out2 %(out2)s
                           -l deduplicated/deduplication_logfile.txt
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
    
    chunksize = int(P.PARAMS['chunksize'])
    compression_level = int(P.PARAMS['compression'])
    fn = os.path.join('split_fastq', os.path.basename(outfile).replace('_0.fastq.gz', ''))

    split_counter = 0
    for read_counter, read in enumerate(FastxFile(infile)):

        processed_read = f'{read}\n'.encode()

        if read_counter == 0:
            out_name = f'{fn}_{split_counter}.fastq.gz'
            out_handle = gzip.open(out_name, 'wb', compresslevel=compression_level)
            out_handle.write(processed_read)
        
        elif read_counter % chunksize == 0:
            split_counter += 1
            out_handle.close()
            out_name = f'{fn}_{split_counter}.fastq.gz'
            out_handle = gzip.open(out_name, 'wb', compresslevel=compression_level)
            out_handle.write(processed_read)       
        else:
            out_handle.write(processed_read)
        
    out_handle.close()    


@follows(mkdir('digest'), split_fastq)          
@transform('split_fastq/*.fastq.gz', 
                   regex(r'split_fastq/(.*).extendedFrags_(\d+).fastq.gz'), 
                    r'digest/\1.digest_flashed_\2.fastq.gz')
def digest_flashed_reads(infile, outfile):
    '''in silico restriction enzyme digest
        Need to remove reads less than 22bp long'''
    statement = '''python %(scripts_dir)s/digest_fastq.py 
                   -o %(outfile)s 
                   -l %(outfile)s.log
                   -r %(ccanalyser_re)s
                   -m 22
                   -c %(compression)s
                   flashed
                   -i %(infile)s'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])

@follows(combine_reads, split_fastq)
@collate('split_fastq/*.fastq.gz', 
               regex(r'split_fastq/(.*).notCombined_[12]_(\d+).fastq.gz'), 
               r'digest/\1.digest_pe_\2.fastq.gz')
def digest_pe_reads(infiles, outfile):
    '''in silico restriction enzyme digest
        Need to remove reads less than 22bp long'''
    
    fq1, fq2 = infiles
    
    statement = '''python %(scripts_dir)s/digest_fastq.py 
                   -l %(outfile)s.log 
                   -r %(ccanalyser_re)s
                   -o %(outfile)s
                   -c %(compression)s
                   -m 22
                   unflashed
                   -1 %(fq1)s 
                   -2 %(fq2)s  
                   '''
    
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'])


@follows(mkdir('bam'))
@transform([digest_flashed_reads, digest_pe_reads], 
           regex(r'digest/(.*).fastq.gz'), 
           r'bam/\1.bam')
def align_reads(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    options = ''
    if P.PARAMS['bowtie2_options']:
        options = P.PARAMS['bowtie2_options']
    statement = '''bowtie2 -x %(bowtie2_index)s -U %(infile)s 
                    -p %(threads)s %(options)s 
                    | samtools view -bS > %(outfile)s 2> %(outfile)s.log'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_threads=P.PARAMS['threads'],
          job_memory='4G')


@transform(align_reads, 
           regex(r'(.*).bam'), 
           r'\1.picard.metrics')
def mapping_qc (infile, outfile):
    sorted_file = infile.replace('.bam', '.sorted.bam')
    statement = ('samtools sort %(infile)s -o %(sorted_file)s &> %(sorted_file)s.log ' 
                 '&& picard CollectAlignmentSummaryMetrics '
                 'R=%(genome_fasta)s I=%(sorted_file)s O=%(outfile)s '
                 '&> %(outfile)s.log')
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='20G')


@merge(mapping_qc, 'report/mapping_report.html')
def mapping_multiqc (infile, outfile):
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc bam/ -o report -n mapping_report.html'''
    P.run(statement, 
          job_queue=P.PARAMS['queue'], 
          job_memory='16G')


@follows(mkdir('ccanalyser'))
@transform(align_reads, 
           regex(r'bam/(.*).bam'), 
           r'ccanalyser/\1.bam.bed')
def bam_to_bed(infile, outfile):
    statement =  '''bedtools bamtobed 
                    -i %(infile)s > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(P.PARAMS["ccanalyser_capture"], 
           regex(P.PARAMS["ccanalyser_capture"]), 
           r'ccanalyser/exclude.bed')
def build_exclusion_bed(infile, outfile):
    statement =  '''bedtools slop  
                    -i %(infile)s -g %(genome_fai)s -b %(ccanalyser_exclude_window)s
                    | bedtools subtract -a - -b %(ccanalyser_capture)s
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'), 
           r'ccanalyser/\1.annotation.capture.count')
def capture_intersect_count(infile, outfile):
    '''Intersect reads with capture and exclusion files.
    report count of overlaps for each input bed using -C '''
    statement =  '''bedtools intersect -c -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\tcapture_count'
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'), 
           r'ccanalyser/\1.annotation.exclude.count')
def exclusion_intersect_count(infile, outfile):
    '''Intersect reads with capture and exclusion files.
    report count of overlaps for each input bed using -C '''
    statement =  '''bedtools intersect -c
                    -a %(infile)s -b ccanalyser/exclude.bed
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\texclusion_count'
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'), 
           r'ccanalyser/\1.annotation.capture')
def capture_intersect(infile, outfile):
    statement =  '''bedtools intersect -loj -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\tcapture'
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'), 
           r'ccanalyser/\1.annotation.exclude')
def exclusion_intersect(infile, outfile):
    statement =  '''bedtools intersect -loj
                    -a %(infile)s -b ccanalyser/exclude.bed
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\texclusion'
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])


@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'),
           add_inputs(digest_genome),
           r'ccanalyser/\1.annotation.re')
def re_intersect(infiles, outfile):   
    bam, genome = infiles        
    statement =  '''bedtools intersect -loj
                    -a %(bam)s -b %(genome)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\trestriction_fragment'
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])

@transform(bam_to_bed, 
           regex(r'ccanalyser/(.*).bam.bed'), 
           r'ccanalyser/\1.annotation.blacklist.count')
def blacklist_intersect_count(infile, outfile):
    '''Intersect reads with blacklisted regions.
    report count of overlaps for each input bed using -C '''
    statement =  '''bedtools intersect -c 
                    -a %(infile)s -b %(ccanalyser_blacklist)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\tblacklist' 
                    > %(outfile)s 2> %(outfile)s.log''' 
    P.run(statement, job_queue=P.PARAMS['queue'])

@collate([capture_intersect_count, exclusion_intersect_count,
         capture_intersect, exclusion_intersect,
         re_intersect, blacklist_intersect_count], 
    regex(r'ccanalyser/(.*).annotation.*'),
    r"ccanalyser/\1.annotations.tsv")
def merge_annotations(infiles, outfile):
    '''merge all intersections into a single file '''
    inlist = " ".join(infiles)
    statement = '''python %(scripts_dir)s/join_tsv.py  
                        -f read_name
                        -o %(outfile)s
                        -i %(inlist)s'''
    P.run(statement,
          job_queue=P.PARAMS['queue'])

@merge(merge_annotations, 'ccanalyser/merged_annotations.tsv')
def combine_annotations(infiles, outfile):
    
    fnames = ' '.join(infiles)
    statement = '''cat %(fnames)s > %(outfile)s'''
    
    P.run(statement,
          job_queue=P.PARAMS['queue'])


@follows(combine_annotations)
@transform(align_reads, 
           regex(r'bam/(.*).bam'), 
           add_inputs("ccanalyser/merged_annotations.tsv"),
           r'ccanalyser/\1.reporter.bed')
def ccanalyser(infiles, outfile):
    bam, annotations = infiles
    outfile_base = outfile.replace('.reporter.bed', '') 
    statement =  '''python %(scripts_dir)s/ccanalyser.py 
                    -i %(bam)s
                    -a %(annotations)s
                    -o %(outfile_base)s ''' 
    P.run(statement,
          job_queue=P.PARAMS['queue'], 
          job_memory=P.PARAMS['ccanalyser_memory'])


@transform(ccanalyser, suffix('.reporter.bed'), '.bedgraph')
def build_bedgraph(infile, outfile):
    '''Intersect reporters with genome restriction fragments to create bedgraph'''
    re_map = 'ccanalyser/' + os.path.basename(P.PARAMS['genome_fasta']).replace('.fasta', '.digest.bed')
    statement = '''bedtools annotate -counts -i %(re_map)s 
                    -files %(infile)s > %(outfile)s'''
    P.run(statement,
          job_queue=P.PARAMS['queue'])
    

#@follows(ccanalyser)
#@transform('ccanalyser/*.stats', regex(r'ccanalyser/(.*).stats'), r'report/\1.html')
#def build_report(infile, outfile):
#    '''Run jupyter notebook for reporting and plotting'''

if __name__ == "__main__":
    sys.exit( P.main(sys.argv))
