from flask import Blueprint, render_template, abort
from ... import hiv


download = Blueprint('download', __name__,
                     url_prefix='/download',
                     #static_folder='static',
                     template_folder='templates')



# Proxy view factory for static files, for download
# NOTE: this is slightly redundant as one can directly access the
# static files, but allows us for abstraction if needed, by switching
# off static file serving and changing
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

    abort(404)
