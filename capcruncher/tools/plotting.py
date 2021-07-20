import os
import sys
from typing import Union

import cooler
from iced import normalization
import pandas as pd
import numpy as np
import iced
import coolbox.api as cb
from coolbox.utilities import get_coverage_stack, get_feature_stack
from matplotlib.patches import Polygon
from pybedtools import BedTool
from coolbox.core.track import Track


class CCMatrix(cb.Cool):
    def __init__(self, 
                 cooler_fn: os.PathLike,
                 binsize: 5000,
                 capture_name: str,
                 remove_viewpoint=False, 
                 **kwargs):
        
        self.cooler_fn = cooler_fn
        self.binsize = binsize
        self.capture_name = capture_name
        self.remove_viewpoint = remove_viewpoint
        self.properties = dict()
        self.properties.update(kwargs)
        

        if not self._cooler_store_has_binsize:
            raise ValueError(f'Capture probe {capture_name} or resolution {binsize} not found in supplied file.')

        self.cooler_obj = cooler.Cooler(f'{cooler_fn}::{capture_name}/resolutions/{binsize}')
        self.capture_bins = self.cooler_obj.info['metadata']['capture_bins']

    def _cooler_store_has_binsize(self):
        clrs = cooler.fileops.list_coolers(self.cooler_fn)
        expected_path = f'{self.capture_name}/resolutions/{self.binsize}'
        
        if expected_path in clrs:
            return True
    
    def get_matrix(self, coordinates, field='count'):
        matrix = self.cooler_obj.matrix(field=field, balance=False).fetch(coordinates)

        offset = self.cooler_obj.offset(coordinates)
        capture_bins = [(bin - offset) for bin in self.capture_bins]

        if self.remove_viewpoint:
            matrix[capture_bins, :] = 0
            matrix[:, capture_bins] = 0
        
        return matrix 
    
    def get_matrix_normalised(self, coordinates, normalization_method=None, **normalisation_kwargs):

        methods_stored = {'n_interactions': 'count_n_interactions_norm',
                          'n_rf_n_interactions': 'count_n_rf_n_interactions_norm'
                          }
        
        if normalization_method == 'raw':
            matrix_normalised = self.get_matrix(coordinates)
        
        elif normalization_method in methods_stored:
            matrix_normalised =  self.get_matrix(coordinates, field=methods_stored[normalization_method])
        
        elif normalization_method == 'ice':
            matrix = self.get_matrix(coordinates) # Get raw matrix
            matrix_ice = iced.normalization.ICE_normalization(matrix, **normalisation_kwargs) # Get iced matrix
            matrix_normalised = matrix_ice / int(self.cooler_obj['metadata']['n_cis_interactions']) #Correct for number of interactions
        
        else:
            raise ValueError(f'Incorrect normalisation specified choose from: {" ".join(["raw", *methods_stored.keys(),"ice"])}')

        return matrix_normalised
    
    def fetch_data(self, gr: cb.GenomeRange, **kwargs) -> np.ndarray:

        norm = self.properties.get('normalization', 'raw')
        return self.get_matrix_normalised(f'{gr.chrom}:{gr.start}-{gr.end}', normalization_method=norm, **kwargs)
class CCBigWig(cb.BigWig):
    def __init__(self, file, mode='fragment', **kwargs):
        
        self.file = file
        self.mode = mode
        super(CCBigWig, self).__init__(file, **kwargs)
    
    def fetch_data(self, gr, **kwargs):
        
        if not self.mode == 'fragment':
            data = super(CCBigWig, self).fetch_data(gr, **kwargs)
        else:
            data = self.bw.fetch_intervals(gr.chrom, gr.start, gr.end)
        
        return data
    
    
    def plot_fragments(self, ax, gr, **kwargs):
        data = self.fetch_data(gr, **kwargs)
        alpha = self.properties.get('alpha', 1.0)
        threshold = self.properties.get("threshold", 0)
        offset = gr.start
        bp_proportion = 1 / (data['end'].max() - data['start'].min())

        for row in data.itertuples():
            pg = Polygon([(row.start, 0), 
                          (row.start, row.value), 
                          (row.end, row.value),
                          (row.end, 0)], 
                         color=self.properties['color'])
            ax.add_patch(pg)


        ax.set_ylim(0, data['value'].max())
        ymin, ymax = self.adjust_plot(ax, gr)
        self.plot_data_range(ax, ymin, ymax, self.properties['data_range_style'], gr)
    
    def plot(self, ax, gr, **kwargs):
        
        if not self.mode == 'fragment':
            super(CCBigWig, self).plot(ax, gr, **kwargs)
        else:
            self.plot_fragments(ax, gr, **kwargs)        
class CCBigWigCollection(cb.BigWig):
    
    DEFAULT_PROPERTIES = {
        "style": 'line',
        "fmt": "-",
        "line_width": 2.0,
        "size": 10,
        "color": "#a6cee3",
        "threshold_color": "#ff9c9c",
        "threshold": "inf",
        "cmap": "bwr",
        "alpha": 1.0,
        "orientation": None,
        "data_range_style": "y-axis",
        "min_value": "auto",
        "max_value": "auto"
    }
    
    def __init__(self, files: list, exclusions: str = None, **kwargs):
        
        self.file_names = files
        self.exclusions = exclusions
        self.bws = [cb.BigWig(fn) for fn in files]
        self.properties = {'files': self.file_names}
        self.properties['name'] = 'BigWigCollection'
        self.properties.update(kwargs)
        self.properties.update(CCBigWigCollection.DEFAULT_PROPERTIES.copy())
        
        
        self.coverages = []

        # load features from global feature stack
        features_stack = get_feature_stack()
        for features in features_stack:
            self.properties.update(features.properties)

        # load coverages from global coverages stack
        coverage_stack = get_coverage_stack()
        for coverage in coverage_stack:
            self.coverages.append(coverage)
               
    def fetch_data(self, gr, **kwargs):
        
        datasets = [bw.bw.fetch_intervals(gr.chrom, gr.start, gr.end)
                         .set_index(['chrom', 'start', 'end'])
                         .rename(columns={'value': os.path.basename(bw.properties['file'])}) 
                    for bw in self.bws]
        df = datasets[0].join(datasets[1:])
        df_summary = df.assign(mean=df.mean(axis=1), std=df.std(axis=1)).reset_index()
        
        
        intervals_to_bp = []
        for interval in df_summary.itertuples():
            interval_len = (interval.end - interval.start)

            interval_positions = np.arange(interval_len) + interval.start
            scores_mean = np.repeat(interval.mean, interval_len)
            scores_std = np.repeat(interval.std, interval_len)

            intervals_to_bp.append(np.vstack([interval_positions, scores_mean, scores_std]).T)

        
        df_intervals = pd.DataFrame(np.concatenate(intervals_to_bp), columns=['bp', 'mean', 'std'])
        
                
        if self.exclusions:
            df_intervals = pd.concat([df_intervals, self.fetch_exluded_regions(gr)]).sort_values('bp')
        
        if self.properties.get('smooth_window'):
            from scipy.signal import savgol_filter
            df_intervals['mean_smoothed'] = savgol_filter(df_intervals['mean'], 
                                                          window_length=self.properties.get('smooth_window', 1001), 
                                                          polyorder=self.properties.get('polyorder', 1))
        
        return df_intervals
    
    
    def fetch_exluded_regions(self, gr):
        
        excluded_tabix = BedTool(self.exclusions).tabix(force=True)
        df_excluded = excluded_tabix.tabix_intervals(f'{gr.chrom}:{gr.start}-{gr.end}').to_dataframe()
        
        intervals_to_bp = []
        for interval in df_excluded.itertuples():
            interval_len = (interval.end - interval.start)

            interval_positions = np.arange(interval_len) + interval.start
            scores_nan = np.repeat(np.nan, interval_len)
            intervals_to_bp.append(interval_positions)
           
        
        df_intervals = pd.Series(np.concatenate(intervals_to_bp)).to_frame('bp')
        df_intervals['mean'] = np.nan
        df_intervals['std'] = np.nan
        
        return df_intervals
    
    
    def plot(self, ax, gr, **kwargs):
        
        data = self.fetch_data(gr, **kwargs)
    
        line_width = self.properties.get('line_width', 1)
        color = self.properties.get('color')
        alpha = self.properties.get("alpha", 1.0)
        
        if self.properties.get('smooth_window'):
            scores = data['mean_smoothed']
            
        else:
            scores = data['mean']
                   
            
        ax.plot(data['bp'], scores)
        ax.fill_between(data['bp'], scores - data['std'], scores + data['std'], alpha=0.3)
        
        
        ax.set_ylim(0, round(scores.max() + data['std'].max()))
        ymin, ymax = self.adjust_plot(ax, gr)
        self.plot_data_range(ax, ymin, ymax, self.properties['data_range_style'], gr)
class ScaleBar(Track):
    def __init__(self, **kwargs):
        self.properties = dict()
        self.properties['name'] = 'Scale'
        self.properties.update(kwargs)
    
    def fetch_data(self, **kwargs):
        pass
    
    
    def get_appropriate_scale(self, length):
        
        if length <= 1e3:
            scale = 1e2
        elif 1e3 < length < 1e4:
            scale = 1e3
        elif 1e4 < length < 1e5:
            scale = 1e4
        elif 1e5 < length < 1e6:
            scale = 1e5
        elif 1e6 < length < 1e7:
            scale = 1e6
        elif 1e7 < length < 1e8:
            scale = 1e8
        
        return scale
    
    def plot(self, ax, gr, **kwargs):
        
        position  = self.properties.get('position', 'left')
        y_midpoint = 0.5
        
        if self.properties.get('scale_distance'):
            scale_distance = self.properties['scale_distance']
        
        else:
            scale_distance = self.get_appropriate_scale(gr.end - gr.start)
        
        # Determine x start and end based on position
        if position == 'left':
            x0 = gr.start
            x1 = x0 + scale_distance
        elif position == 'right':
            x0 = gr.end - scale_distance
            x1 = gr.end
        else:
            raise ValueError('Position can only be "left" or "right"')
        
        # Plot scale bar
        ax.plot([x0, x1], [y_midpoint, y_midpoint], color='black')
        ax.plot([x0, x0], [y_midpoint - 0.1, y_midpoint + 0.1], color='black', lw=1)
        ax.plot([x1, x1], [y_midpoint - 0.1, y_midpoint + 0.1], color='black', lw=1)
        
        # Add annotation
        from capcruncher.utils import get_human_readable_number_of_bp
        scale_distance_human_readable = get_human_readable_number_of_bp(scale_distance)
        
        
        ax.text((x0 + (scale_distance / 2)), y_midpoint - 0.1, scale_distance_human_readable, ha='center', va='center')
        
        
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
class SimpleBed(cb.BED):
    def __init__(self, file: str, **kwargs):
        
        self.file = file
        self.properties = dict()
        self.properties['name'] = 'BlockBed'
        self.properties.update(kwargs)
    
    def fetch_data(self, gr):
        
        bt = BedTool(self.file)
        bt_tabix = bt.tabix(force=True)
        
        return bt_tabix.tabix_intervals(f'{gr.chrom}:{gr.start}-{gr.end}').to_dataframe()
    
    
    def plot(self, ax, gr, **kwargs):
        
        data = self.fetch_data(gr)
        y_midpoint = 0.5
        
        for interval in data.itertuples():
            pg = Polygon([(interval.start, y_midpoint - 0.1), 
                     (interval.start, y_midpoint + 0.1), 
                     (interval.end, y_midpoint + 0.1), 
                     (interval.end, y_midpoint - 0.1)],
                      color=self.properties.get('color', 'black'))
            
            
            if hasattr(interval, 'name') and not self.properties.get('no_annotation'):
                interval_midpoint = interval.start + ((interval.end - interval.start) / 2)
                ax.text(interval_midpoint, y_midpoint - 0.1, interval.name, ha='center', va='center')
            
            ax.add_patch(pg)
            
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(0, 1)
            
        
    