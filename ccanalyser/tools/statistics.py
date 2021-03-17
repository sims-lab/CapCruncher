import numpy as np
import pandas as pd
import re
from multiprocessing import Queue
from collections import defaultdict, Counter


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

class DigestionStatCollector():
    def __init__(self, statq: Queue, n_subprocess: int = 1) -> None:

        self.statq = statq
        self.n_subprocesses = n_subprocess
        self.n_subprocesses_terminated = 0
        self.stat_dict = None
    
    def _set_up_stats_dict(self, stat_sample):
        
        stat_dict = dict()
        
        for rn in stat_sample:
            stat_dict[rn] = defaultdict(list)
        
        return stat_dict

    def _append_statistics(self, frag_stats, stats_dict):
        
        for read_number, read_stats in frag_stats.items():
            for stat_type, count in read_stats.items():
                stats_dict[read_number][stat_type].append(count)
        
        return stats_dict


    def _format_stats_dict(self, stats_dict):
        for read_number, read_number_stats in stats_dict.items():
            for stat_type, counts in read_number_stats.items():
                stats_dict[read_number][stat_type] = Counter(counts)
        
        return stats_dict


    def get_collated_stats(self) -> dict:
              
        stats = self.statq.get()
        stat_dict = self._set_up_stats_dict(stats[0])
        
        while self.n_subprocesses_terminated < self.n_subprocesses:
            
           
            if stats == 'END':
                self.n_subprocesses_terminated += 1
                continue
        
            else:
                for fragment_stats in stats:
                    stat_dict = self._append_statistics(fragment_stats, stat_dict)
            
            stats = self.statq.get()

        return self._format_stats_dict(stat_dict)
  
class DigestionStatistics():
    def __init__(self,
                 sample: str,
                 read_type='pe',
                 read_number=1,  
                 slices_unfiltered: Counter = None,
                 slices_filtered: Counter = None):

        self.sample = sample
        self.read_type = read_type
        self.read_number = read_number
        self.slices_unfiltered = slices_unfiltered
        self.slices_filtered = slices_filtered


        self.n_unfiltered_slices = self._get_n_unfiltered_slices()
        self.n_unfiltered_reads = self._get_n_unfiltered_reads()
        self.n_filtered_slices = self._get_n_filtered_slices()
        self.n_filtered_reads = self._get_n_filtered_reads()
        self.filtered_histogram = self._get_filtered_histogram()
        self.unfiltered_histogram = self._get_unfiltered_histogram()
        self.slice_summary = self._get_slice_summary_df()
        self.read_summary = self._get_read_summary_df()

    def _get_n_unfiltered_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_unfiltered.items())
    
    def _get_n_unfiltered_reads(self):
        return sum(self.slices_unfiltered.values())
    
    def _get_n_filtered_slices(self):
        return sum(n_slices * count for n_slices, count in self.slices_filtered.items())

    def _get_n_filtered_reads(self):
        return sum(v for k, v in self.slices_filtered.items() if not k == 0)
    
    def _get_unfiltered_histogram(self):
        df = pd.DataFrame()
        df['number_of_slices'] = self.slices_unfiltered.keys()
        df['count'] = self.slices_unfiltered.values()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        return df[['sample', 'read_type', 'read_number', 'number_of_slices', 'count']].sort_values('number_of_slices')
    
    def _get_filtered_histogram(self):
        df = pd.DataFrame()
        df['number_of_slices'] = self.slices_filtered.keys()
        df['count'] = self.slices_filtered.values()
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        return df[['sample', 'read_type', 'read_number', 'number_of_slices', 'count']].sort_values('number_of_slices')
    
    def _get_slice_summary_df(self):
        df = pd.DataFrame()
        df['stat_type'] = ['unfiltered', 'filtered']
        df['stat'] = [self.n_unfiltered_slices, self.n_filtered_slices]
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stage'] = 'digestion'
        return df[['sample','stage', 'read_type', 'read_number', 'stat_type', 'stat']]
    
    def _get_read_summary_df(self):
        df = pd.DataFrame()
        df['stat_type'] = ['unfiltered', 'filtered']
        df['stat'] = [self.n_unfiltered_reads, self.n_filtered_reads]
        df['sample'] = self.sample
        df['read_type'] = self.read_type
        df['read_number'] = self.read_number
        df['stage'] = 'digestion'
        return df[['sample', 'stage','read_type', 'read_number', 'stat_type', 'stat']]
    
    def __add__(self, other):       
        self.unfiltered_histogram = pd.concat([self.unfiltered_histogram, other.unfiltered_histogram]).reset_index(drop=True)
        self.filtered_histogram = pd.concat([self.filtered_histogram, other.filtered_histogram]).reset_index(drop=True)
        self.slice_summary = pd.concat([self.slice_summary, other.slice_summary]).reset_index(drop=True)
        self.read_summary = pd.concat([self.read_summary, other.read_summary]).reset_index(drop=True)
        return self

def collate_histogram_data(fnames):
    return (
        pd.concat([pd.read_csv(fn) for fn in fnames])
        .groupby(["sample", "read_type", "read_number", "number_of_slices"])["count"]
        .sum()
        .reset_index()
        .sort_values(["sample", "read_type", "number_of_slices"])
    )

def collate_read_data(fnames):

    df = (pd.concat([pd.read_csv(fn) for fn in fnames])
        .query('(read_number == 0) or (read_number == 1)')
        .groupby(["sample", "stage", "read_type", "read_number", "stat_type"])["stat"]
        .sum()
        .reset_index()
        )
    
    # PE reads result in duplicate statitsics as considered independently
    # Need to halve all of the PE statistics but not for deduplication/digestion which are correct
    # df = df.assign(stat=np.where((df['read_type'] == 'pe') & 
    #                              (df['read_number'] == 0) & 
    #                              (df['stage'] != 'deduplication') &
    #                              (df['stage'] != 'digestion'), 
    #                              df['stat'] // 2, 
    #                              df['stat']))
    
    return df.sort_values(["sample", "stat"], ascending=[True, False])


def collate_slice_data(fnames):

    df = pd.concat([pd.read_csv(fn) for fn in fnames])
    aggregations = {col: 'sum' if not 'unique' in col else 'max' for col in df.columns
                    if not col in ['sample', 'stage', 'read_type']}

    return (df.groupby(["sample", "stage", "read_type"])
              .agg(aggregations)
              .reset_index()
              .sort_values(["sample", "unique_slices"], ascending=[True, False])
        )

def collate_cis_trans_data(fnames):

    return (pd.concat([pd.read_csv(fn) for fn in fnames])
              .groupby(['sample', 'capture', 'read_type', 'cis/trans'])
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
        