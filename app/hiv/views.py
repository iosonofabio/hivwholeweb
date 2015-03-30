from flask import render_template, abort
from . import hiv

# Proxy view factory for static files, for download
# NOTE: this is slightly redundant as one can directly access the
# static files, but allows us for abstraction if needed, by switching
# off static file serving and changing
# TODO: fix robots.txt
@hiv.route('/download/<path:path>')
def data_proxy(path):
    import re

    sf = hiv.config['DATA_SUBFOLDER']

    fields = re.split('_|\.', path)
    dtype = fields[0]

    # Switch by data type
    # FIXME: we should be using data models to manage these paths, not
    # hardcode them!
    if dtype == 'tree':
        if len(fields) < 4:
            abort(404)

        pname = fields[1]
        fragment = fields[2]
        format = fields[3]
        fn = sf+'/trees/consensi_tree_'+pname+'_'+fragment+'.'+format

        # NOTE: send_static_file will guess the MIME type from the
        # extension of this URL (!!)
        return hiv.send_static_file(fn)

    elif dtype == 'physio':
        if len(fields) < 4:
            abort(404)

        dstypes = {'vl': 'viral_load', 'cc': 'cell_count'}
        dstype = dstypes[fields[1]]
        pname = fields[2]
        format = fields[3]
        fn = sf+'/physiological/'+dstype+'_'+pname+'.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'genome':
        if len(fields) < 3:
            abort(404)

        pname = fields[1]
        format = fields[2]
        fn = sf+'/sequences/reference_'+pname+'_genomewide.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'act':
        if len(fields) < 3:
            abort(404)
        pname = fields[1]
        format = fields[2]
        fn = sf+'/one_site/allele_counts_'+pname+'_genomewide.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'divdiv':
        if len(fields) < 5:
            abort(404)

        dstypes = {'ds': 'diversity', 'dg': 'divergence'}
        dstype = dstypes[fields[1]]
        pname = fields[2]
        fragment = fields[3]
        format = fields[4]
        fn = sf+'/divdiv/'+dstype+'_'+pname+'_'+fragment+'.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'haplotypes':
        if len(fields) < 4:
            abort(404)

        pname = fields[1]
        region = fields[2]
        format = fields[3]
        fn = sf+'/alignments/alignments_'+pname+'_'+region+'.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'consensi':
        if len(fields) < 4:
            abort(404)

        pname = fields[1]
        fragment = fields[2]
        format = fields[3]
        fn = sf+'/alignments/consensi_alignment_'+pname+'_'+fragment+'.'+format
        return hiv.send_static_file(fn)

    elif dtype == 'reads':
        if len(fields) < 5:
            abort(404)
        pname = fields[1]
        tind = fields[2]
        fragment = fields[3]
        format = fields[4]
        fn = sf+'/patients/'+pname+'/samples/'+tind+'/'+fragment+'.'+format
        print fn
        return hiv.send_static_file(fn)

    abort(404)
