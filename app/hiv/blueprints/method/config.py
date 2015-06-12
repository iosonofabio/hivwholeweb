# method sections
sections = [{'id': 'physio',
             'name': 'Viral load and CD4+ counts',
             'url': '/physiological/',
             'thumbnail': '/static/images/icons/viral_load.svg',
             'thumbnail_white': '/static/images/icons/viral_load_white.svg'},
            {'id': 'genome',
             'name': 'Genome sequences',
             'url': '/genomes/',
             'thumbnail': '/static/images/icons/genome.svg',
             'thumbnail_white': '/static/images/icons/genome_white.svg'},
            {'id': 'consensi',
             'name': 'Consensus sequences',
             'url': '/consensi/',
             'thumbnail': '/static/images/icons/consensi.svg',
             'thumbnail_white': '/static/images/icons/consensi_white.svg'},
            {'id': 'tree',
             'name': 'Phylogenetic trees',
             'url': '/trees/',
             'thumbnail': '/static/images/icons/tree.svg',
             'thumbnail_white': '/static/images/icons/tree_white.svg'},
            {'id': 'divdiv',
             'name': 'Divergence and diversity',
             'url': '/divdiv/',
             'thumbnail': '/static/images/icons/diversity.svg',
             'thumbnail_white': '/static/images/icons/diversity_white.svg'},
            {'id': 'divdiv_local',
             'name': 'Local divergence and diversity',
             'url': '/divdivLocal/',
             'thumbnail': '/static/images/icons/divdiv_local.svg',
             'thumbnail_white': '/static/images/icons/divdiv_local_white.svg'},
            {'id': 'cov',
             'name': 'Coverage',
             'url': '/coverage/',
             'thumbnail': '/static/images/icons/coverage.svg',
             'thumbnail_white': '/static/images/icons/coverage_white.svg'},
            {'id': 'af',
             'name': 'Single nucleotide polymorphisms',
             'url': '/singleNucleotidePolymorphisms/',
             'thumbnail': '/static/images/icons/allele_frequencies.svg',
             'thumbnail_white': '/static/images/icons/allele_frequencies_white.svg'},
            {'id': 'haplo',
             'name': 'Haplotypes',
             'url': '/haplotypes/',
             'thumbnail': '/static/images/icons/haplotype_frequencies.svg',
             'thumbnail_white': '/static/images/icons/haplotype_frequencies_white.svg'},
           ]
for section in sections:
    if 'url' in section:
        section['url_complete'] = '/method'+section['url']


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
