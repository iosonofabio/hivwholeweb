# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/11/14
content:    Support module that contains functions that should not be required.
            This module might disappear from one moment to the next, so do not
            rely on it heavily.
'''
# Functions
def correct_genbank_features_load(record):
    '''Prepare genbank file after loading with the correct features'''
    for feat in record.features:
        try:
            feat.id = feat.qualifiers['note'][-1]
        except KeyError, IndexError:
            pass
