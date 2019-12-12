#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 13:47:20 2019
@author: asmith
"""

'''
TRACK HUB STRUCTURE:

hub_folder/
        hub.txt
        genomes.txt
        genome_assembly/
                        trackDB.txt
                        data.bw
hub.txt:

hub UCSCHub:
shortLabel UCSC Hub
longLabel UCSC Genome Informatics Hub for human DNase and RNAseq data
genomesFile genomes.txt
email genome@soe.ucsc.edu
descriptionUrl ucscHub.html

genomes.txt:
genome assembly_database_1
trackDb assembly_1_path/trackDb.txt
metaTab assembly_1_path/tabSeparatedFile.txt

trackDB.txt:
track track_name
bigDataUrl track_data_URL
shortLabel short_label
longLabel long_label
type track_type 
'''
import sys
import os
import seaborn as sns
from cgatcore import pipeline as P
from ruffus import *

# Read in parameter file
P.get_parameters('capturec_pipeline.yml')

hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])

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



def get_track_data(fn):
    return {'track': fn,
                'bigDataUrl': f'http://userweb.molbiol.ox.ac.uk/{(os.path.join(assembly_dir, fn)).lstrip("/")}',
                'shortLabel': fn,
                'longLabel': fn,
                'type': f'{fn.split(".")[-1]}',
                }

@follows(generate_hub_metadata, mkdir(assembly_dir))
@merge(r'bigwigs/*.bigWig', f'{assembly_dir}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):
    
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
@transform(r'bigwigs/*.bigWig', regex('bigwigs/(.*).bigWig'),  f'{assembly_dir}' + r'/\1.bigWig')
def link_bigwigs(infile, outfile):
    
    infile_fp = os.path.abspath(infile)
    os.symlink(infile_fp, outfile)








if __name__ == "__main__":
    sys.exit( P.main(sys.argv))