#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Converts the tsv output of ccanalyser.py to a bedgraph by intersecting
the reporter slices with a bed file. The pipeline uses the generated restriction
map of the genome  (ccanalysis/restriction_enzyme_map/genome.digest.bed.gz).

'''

import argparse
import os
import sys
import pandas as pd
import numpy as np
from pybedtools import BedTool
import tempfile

class CCBedgraph(object):
    def __init__(self, path=None, df=None):

        self.fn = path
        self.df = df

        if self.fn:
            self.df = pd.read_csv(self.fn,
                                  sep='\t',
                                  header=None,
                                  names=['chrom', 'start', 'end', 'score'])

    @property
    def score(self):
        return self.df.rename(columns={'score': self.fn})[self.fn]

    @property
    def coordinates(self):
        return self.df.loc[:, 'chrom':'end']

    def to_bedtool(self):
        return self.df.pipe(BedTool.from_dataframe)

    def normalise_by_region(self, region_coords, n_read_scale=1e6):

        region = BedTool('\t'.join(region_coords), from_string=True)
        intervals_in_region = (self.to_bedtool()
                                   .intersect(region, wa=True)
                                   .to_dataframe()
                              )
        scale_factor = intervals_in_region.iloc[:, -1].sum() / n_read_scale
        self.df['score'] = self.df['score'] / scale_factor
        return self

    def to_file(self, path):
        self.df.to_csv(path, sep='\t', header=None, index=None)

    def __add__(self, other):
        if isinstance(other, CCBedgraph):
            self.df['score'] = self.df['score'] + other.df['score']
            return self

        elif isinstance(other, [np.ndarray, pd.Series, int]):
            self.df['score'] = self.df['score'] + other
            return self

        else:
            return NotImplementedError()

    def __sub__(self, other):
        if isinstance(other, CCBedgraph):
            self.df['score'] = self.df['score'] - other.df['score']
            return self

        elif isinstance(other, [np.ndarray, pd.Series, int]):
            self.df['score'] = self.df['score'] - other
            return self

        else:
            return NotImplementedError()

class CCBedgraphCollection(object):
    def __init__(self, bedgraphs):
        self.bdg_fnames = bedgraphs
        self.bedgraphs = [CCBedgraph(bg) for bg in self.bdg_fnames]

    def normalise_bedgraphs(self, how='region', **kwargs):

        if how == 'region':
            self.bedgraphs = [bg.normalise_by_region(**kwargs) for bg in self.bedgraphs]
        else:
            raise NotImplementedError('No other methods')
        return self

    def get_average_bedgraph(self):
        coords = self.bedgraphs[0].coordinates
        scores = [bg.score for bg in self.bedgraphs]
        df_scores = pd.concat(scores, axis=1)
        df_scores_mean = df_scores.mean(axis=1)
        bdg = (
            pd.concat([coords, df_scores_mean], axis=1).pipe(BedTool.from_dataframe).fn
        )
        return CCBedgraph(bdg)


def make_bedgraph(reporters, bed):
    df_reporters = pd.read_csv(reporters, sep='\t')
    
    df_reporters[['reporter_start', 'reporter_end']] = df_reporters[
        ['reporter_start', 'reporter_end']
    ].astype(int)
    
    bt_reporters = BedTool.from_dataframe(
        df_reporters[
            ['reporter_chrom', 'reporter_start', 'reporter_end', 'reporter_read_name']
        ]
    )

    bt_bed = BedTool(bed)
    bedgraph = (bt_bed.intersect(bt_reporters, c=True)
                      .to_dataframe()
                      .sort_values(['chrom', 'start'])
                      .iloc[:, [0,1,2,4]]
                )

    return bedgraph


def main(
    tsv_input,
    bed,
    output='out.bedgraph',
    normalise=None,
    region=[],
    n_reads_scale=1e6,
):

    '''
    Converts reporter tsv to bedgraph

    Args:
        tsv_input - reporter tsv file name
        bed - bed file to intersect with reporters
        output - file name for output file


    '''

    # Generate bedgraph
    bedgraph = CCBedgraph(df=make_bedgraph(tsv_input, bed))

    # Run normalisation if required
    if normalise == 'region':
        bedgraph = bedgraph.normalise_by_region(region, n_read_scale=n_reads_scale)

    bedgraph.to_file(output)
