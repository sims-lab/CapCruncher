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
from cgatcore import pipeline as P
from ruffus import *

# Read in parameter file
P.get_parameters('capturec_pipeline.yml')

hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])

@mkdir(hub_dir)
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

@originate(os.path.join(hub_dir, 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {'genome': P.PARAMS['hub_genome'],
               'trackDb': os.path.join(P.PARAMS['hub_genome'], 'trackDb.txt'),
               }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label}: {info}\n')


@mkdir(assembly_dir)
@merge(r'*.bigWig', f'{assembly_dir}/trackDb.txt')
def generate_trackDb_metadata(infiles, outfile):

    with open(outfile, 'w') as w:
        for fn in infiles:
            fn = os.path.basename(fn)

            w.write(f'track {fn}\n')
            w.write(f'bigDataUrl http://userweb.molbiol.ox.ac.uk/{(os.path.join(assembly_dir, fn)).lstrip("/")}\n')
            w.write(f'shortLabel {fn}\n')
            w.write(f'longLabel {fn}\n')
            w.write(f'type {fn.split(".")[1]}\n')

