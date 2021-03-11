import os
import sys
from typing import Union

import cooler
import pandas as pd
import numpy as np
import iced



class CCMatrix():
    def __init__(self, 
                 cooler_fn: os.PathLike,
                 binsize: 5000,
                 capture_name: str,
                 remove_capture=False):
        
        self.cooler_fn = cooler_fn
        self.binsize = binsize
        self.capture_name = capture_name
        self.remove_capture = remove_capture
        

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

        if self.remove_capture:
            matrix[capture_bins, :] = 0
            matrix[:, capture_bins] = 0
        
        return matrix 
    
    def get_matrix_normalised(self, coordinates, normalisation_method=None, **normalisation_kwargs):

        methods_stored = {'n_interactions': 'count_n_interactions_norm',
                          'n_rf_n_interactions': 'count_n_rf_n_interactions_norm'
                          }
        
        if normalisation_method in methods_stored:
            matrix_normalised =  self.get_matrix(coordinates, field=methods_stored[normalisation_method])
        
        elif normalisation_method == 'ice':
            matrix = self.get_matrix(coordinates) # Get raw matrix
            matrix_ice = iced.normalization.ICE_normalization(matrix, **normalisation_kwargs) # Get iced matrix
            matrix_normalised = matrix_ice / int(self.cooler_obj['metadata']['n_cis_interactions']) #Correct for number of interactions
        else:
            raise ValueError(f'Incorrect normalisation specified choose from: {" ".join([*methods_stored.keys(),"ice"])}')

        return matrix_normalised
            
            
