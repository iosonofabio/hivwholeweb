# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/03/15
content:    Main config file for the hivwholeweb web application.
'''
def find_tmp_folder():
    import os

    # Try standard webserver location first
    fn = '/home/hivwholewebu1/tmp/'
    if os.path.isdir(fn):
        return fn

    # Try OS standard
    fn = '/tmp/'
    if os.path.isdir(fn):
        return fn


DEBUG = False
CSRF_ENABLED = True
SECRET_KEY = 'a gig is a gig is a gig'

BLUEPRINTS = {}

TIMEOUT_TMP = '1d'
TMP_ROOT_FOLDER = find_tmp_folder()

DATA_SUBFOLDER = 'data'


REFERENCES = ['HXB2']

# Only patients from the paper #1 are listed
PATIENTS = ['p1', 'p2', 'p3', 'p5', 'p6', 'p8', 'p9', 'p10', 'p11']

# first 11 colors of d3 category20
COLORS_PATIENTS = {'p1': '#1f77b4',
                   'p2': '#aec7e8',
                   'p3': '#ff7f0e',
                   'p4': '#ffbb78',
                   'p5': '#2ca02c',
                   'p6': '#98df8a',
                   'p7': '#d62728',
                   'p8': '#ff9896',
                   'p9': '#9467bd',
                   'p10': '#c5b0d5',
                   'p11': '#8c564b',
                  }

# Genomic regions: haplotypes, SNP, trees, etc.
REGIONS_HAPLO = ['p17',
                 'RT1',
                 'RT2',
                 'RT3',
                 'p15',
                 'IN1',
                 'IN2',
                 'V3',
                ]

REGIONS_SNP = ['V3',
               'psi',
               'p17',
               'RRE',
              ]

REGIONS_TREE = REGIONS_SNP + [r+'_minor' for r in REGIONS_HAPLO]
