# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the download of data files.

            NOTE: this is a direct access to existing files using Flask's
            'serve_static_file' function, which calls the web server facilities.
            It is faster than the API and delivers in several formats, whereas
            the API is JSON only.

            NOTE: one COULD allow direct access to the file path from the web,
            but this blueprint creates a shim that decouples filenames on the
            hard drives from web URIs, greatly increasing flexibility in case
            file system optimizations are required (e.g. hierarchical directory
            tree). The current example of this is the agnostic behaviour of the
            download view towards underscore (_) VS slash (/).
'''
# Modules
from flask import Blueprint, render_template, abort
from ... import hiv



# Blueprint
download = Blueprint('download', __name__,
                     url_prefix='/download',
                     #static_folder='static',
                     template_folder='templates')



# Views
@download.route('/<path:path>')
def data_proxy(path):
    def send_static_file(fn):
        '''Send static files
        
        NOTE: hiv.send_static_file will guess the MIME type from the
        extension of this URL (!!)
        '''
        if fn[:8] == '/static/':
            fn = fn[8:]
        return hiv.send_static_file(fn)

    # For now accept both underscore and slash as separator,
    # and figure the data format from the extension after .
    import re
    fields = re.split('_|/', path)
    if len(fields) < 2:
        abort(404)
    dtype = fields[0]
    fields = fields[1:]
    if '.' in fields[-1]:
        fmt = fields[-1].split('.')[-1]
        fields[-1] = fields[-1][:-len(fmt) - 1]
    else:
        fmt = None

    # Switch by data type
    # TREES
    if dtype == 'tree':
        if len(fields) < 2:
            abort(404)

        pname = fields[0]
        region = fields[1]
        if fmt is None:
            fmt = 'json'

        from ...models import TreeModel
        fn = (TreeModel(pname, region)
              .get_filename(full=False, format=fmt))
        return send_static_file(fn)


    # VIRAL LOAD/CELL COUNT
    elif dtype in ['viralLoad', 'cellCount']:
        pname = fields[0]
        if fmt is None:
            fmt = 'tsv'

        from ...models import PhysioModel
        fn = (PhysioModel(pname)
              .get_physio_filename(dtype, full=False, format=fmt))
        return send_static_file(fn)


    # GENOMES
    elif dtype in ['genome', 'consensus', 'reference']:
        psname = fields[0]
        if len(fields) < 2:
            region = 'genomewide'
        else:
            region = fields[1]

        from ...models import GenomeModel
        kwargs = {}
        if fmt is not None:
            kwargs['format'] = fmt
        fn = (GenomeModel(psname, region=region)
              .get_genome_filename(full=False, **kwargs))

        return send_static_file(fn)


    # ALLELE COUNT TRAJECTORIES
    elif dtype == 'act':
        pname = fields[0]
        if fmt is None:
            fmt = 'zip'

        from ...models import AlleleFrequencyTrajectoryModel
        fn = (AlleleFrequencyTrajectoryModel(pname)
              .get_allele_counts_filename(full=False, format=fmt))
        return send_static_file(fn)


    # HAPLOTYPES
    elif dtype == 'haplotypes':
        if len(fields) < 2:
            abort(404)

        pname = fields[0]
        region = fields[1]
        if fmt is None:
            fmt = 'fasta'

        from ...models import HaplotypePrecompiledTrajectoryModel
        fn = (HaplotypePrecompiledTrajectoryModel(pname, region)
              .get_haplotype_filename(full=False, format=fmt)) 
        return send_static_file(fn)


    # READS
    elif dtype == 'reads':
        if len(fields) < 3:
            abort(404)
        pname = fields[0]
        tind = fields[1]
        fragment = fields[2]
        if fmt is None:
            fmt = 'bam'

        from ...models import ReadsTableModel
        fn = (ReadsTableModel(pname)
              .get_reads_filename(tind, fragment, full=False, format=fmt))
        return send_static_file(fn)


    # COCOUNTS
    elif dtype == 'cocounts':
        if len(fields) < 3:
            abort(404)
        pname = fields[0]
        tind = fields[1]
        fragment = fields[2]
        if fmt is None:
            fmt = 'zip'

        from ...models import CocountsTableModel
        fn = (CocountsTableModel(pname)
              .get_cocounts_filename(tind, fragment, full=False, format=fmt))
        return send_static_file(fn)


    # HLA
    elif dtype == 'HLA':
        if len(fields) < 1:
            abort(404)
        pname = fields[0]

        from ...models import HLAModel
        fn = HLAModel(pname).get_filename(full=False, format=fmt)
        return send_static_file(fn)


    # default return error code
    abort(404)
