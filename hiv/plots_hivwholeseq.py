# vim: fdm=indent
'''
author:     Fabio Zanini
date:       01/10/14
content:    Test D3 with our HIV data.
'''
def coverage(fragment):
    import sys
    hivwholeseq_path = '/home/fabio/university/phd/sequencing/script/mapping'
    if hivwholeseq_path not in sys.path:
        sys.path.append(hivwholeseq_path)
    general_path = '/usr/lib/python2.7/site-packages'
    if general_path not in sys.path:
        sys.path.append(general_path)
    
    from hivwholeseq.patients.patients import load_sample_sequenced
    sample = load_sample_sequenced('VL96-15555')
    cov = sample.get_coverage(fragment=fragment, PCR=1)

    import numpy as np

    x_data = np.arange(cov.shape[-1])
    y_data = cov
    data = list(y_data)

    return data

