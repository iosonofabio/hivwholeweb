# Sections of the website (for navbar et al)
sections = [{'id': 'tutorial', 'name': 'Tutorial', 'url': '/tutorial/',
             'thumbnail': '/static/images/icons/tutorial.svg'},
            {'id': 'physio', 'name': 'Viral load and CD4+ counts', 'url': '/physiological/',
             'thumbnail': '/static/images/icons/viral_load.svg'},
            {'id': 'genome', 'name': 'Genome sequences', 'url': '/genomes/',
             'thumbnail': '/static/images/icons/genome.svg'},
            {'id': 'consensi', 'name': 'Consensus sequences', 'url': '/consensi/',
             'thumbnail': '/static/images/icons/consensi.svg'},
            {'id': 'tree', 'name': 'Phylogenetic trees', 'url': '/trees/',
             'thumbnail': '/static/images/icons/tree.svg'},
            {'id': 'divdiv', 'name': 'Divergence and diversity', 'url': '/divdiv/',
             'thumbnail': '/static/images/icons/diversity.svg'},
            {'id': 'divdiv_local', 'name': 'Local divergence and diversity', 'url': '/divdiv_local/',
             'thumbnail': '/static/images/icons/divdiv_local.svg'},
            {'id': 'cov', 'name': 'Coverage', 'url': '/coverage/',
             'thumbnail': '/static/images/icons/coverage.svg'},
            {'id': 'af', 'name': 'Allele frequencies', 'url': '/allele_frequencies/',
             'thumbnail': '/static/images/icons/allele_frequencies.svg'},
            {'id': 'haplo', 'name': 'Local haplotypes', 'url': '/local_haplotypes/',
             'thumbnail': '/static/images/icons/haplotype_frequencies.svg'},
            {'id': 'sfs', 'name': 'Site frequency spectra', 'url': '/sfs/',
             'thumbnail': '/static/images/icons/SFS.svg'},
            {'id': 'prop', 'name': 'Propagators', 'url': '/propagators/',
             'thumbnail': '/static/images/icons/propagator.svg'},
            {'id': 'pfix', 'name': 'Fixation probability',
             'thumbnail': '/static/images/icons/fixation_probability.svg'},
           ]


def find_section(**kwargs):
    '''Find a section with this attributes'''
    if not kwargs:
        return ValueError('Please select at least one attribute to query')

    for section in sections:
        for (key, value) in kwargs.iteritems():
            if (key not in section) or (section[key] != value):
                break
        else:
            return section
    else:
        return ValueError('Section not found')
