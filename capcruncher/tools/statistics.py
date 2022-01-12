import pandas as pd
import re
from collections import namedtuple
from capcruncher.utils import read_dataframes


class DeduplicationStatistics():
    def __init__(self,
                 sample: str,
                 read_type: str = 'pe',
                 reads_total: int=0,
                 reads_unique: int=0):

        self.sample = sample
        self.read_type = read_type
        self.reads_total = reads_total
        self.reads_unique = reads_unique
        self.reads_removed = reads_total - reads_unique
    
    @property
    def df(self):
        df = pd.DataFrame()
        df['stat'] = [self.reads_total, self.reads_unique, self.reads_removed]
        df['stat_type'] = ['reads_total', 'reads_unique', 'reads_removed']
        df['read_type'] = self.read_type
        df['read_number'] = 0
        df['stage'] = 'deduplication'
        df['sample'] = self.sample
        return df

DigestionStats = namedtuple('DigestionStats', ['read_type', 'read_number', 'unfiltered', 'filtered'])

def collate_histogram_data(fnames):
    return (
        pd.concat(read_dataframes(fnames))
        .groupby(["sample", "read_type", "read_number", "filtered", "n_slices"])
        .sum()
        .reset_index()
        .sort_values(["sample", "read_type", "n_slices"])
    )

def collate_read_data(fnames):

    df = (pd.concat(read_dataframes(fnames))
        .query('(read_number == 0) or (read_number == 1)')
        .groupby(["sample", "stage", "read_type", "read_number", "stat_type"])["stat"]
        .sum()
        .reset_index()
        )
    
    return df.sort_values(["sample", "stat"], ascending=[True, False])

def collate_slice_data(fnames):

    df = pd.concat(read_dataframes(fnames))
    aggregations = {col: 'sum' if not 'unique' in col else 'max' for col in df.columns
                    if not col in ['sample', 'stage', 'read_type']}

    return (df.groupby(["sample", "stage", "read_type"])
              .agg(aggregations)
              .reset_index()
              .sort_values(["sample", "unique_slices"], ascending=[True, False])
        )

def collate_cis_trans_data(fnames):
    return (pd.concat(read_dataframes(fnames))
              .groupby(['sample', 'viewpoint', 'read_type', 'cis/trans'])
              .sum()
              .reset_index()
              .sort_values(["sample", "read_type", 'count'], ascending=[True, True, False]))

def extract_trimming_stats(fn):
    stat_regexes = {
                'reads_total': re.compile(r'^Total reads processed:\s+([0-9,]+)$'),
                'adapters_removed': re.compile(r'Reads with adapters:\s+([0-9,]+).*'),
                'reads_after_filtering': re.compile(r'Reads written \(passing filters\):\s+([0-9,]+).*'),}
    
    sample_re_match = re.match(r'.*/(.*)_part\d+_(1|2).*', fn)
    
    stats = {}
    stats['sample'] = sample_re_match.group(1)
    stats['read_number'] = sample_re_match.group(2)
    stats['read_type'] = 'pe'
    
    
    with open(fn) as r:
        for line in r:
            for stat_name, pattern in stat_regexes.items():
                regex_match = pattern.match(line)
                if regex_match:
                    stats[stat_name] = int(regex_match.group(1).replace(',', ''))
    
    stats['reads_filtered'] = stats['reads_total'] - stats['reads_after_filtering']
    
    return stats
        