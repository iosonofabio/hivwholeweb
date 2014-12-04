from re import split as resplit
from flask import (render_template, flash, redirect, request, jsonify,
                   make_response, abort)
from . import hiv
from .forms import (PatFragForm, PatForm, LocalHaplotypeForm, TreeForm,
                    PatSingleForm, PatFragSingleForm, PrecompiledHaplotypeForm,
                    PatFragGWForm)
from .models import (TreeModel, PhysioModel, DivdivModel, CoverageModel,
                     GenomeModel, AlleleFrequencyModel, SFSModel,
                     NTemplatesModel,
                     PropagatorModel, DivdivLocalModel, LocalHaplotypeModel)
from .backbone import find_section


@hiv.route('/')
@hiv.route('/index/')
def index():
    return render_template('index.html',
                           title='Home',
                           section_name='Home')


@hiv.route('/impressum/', methods=['GET'])
def impressum():
    return render_template('impressum.html',
                           title='Impressum')


@hiv.route('/contact/', methods=['GET'])
def contact():
    return render_template('contact.html',
                           title='Contact')


@hiv.route('/tutorial/', methods=['GET'])
def tutorial():
    return render_template('tutorial.html',
                           title='Tutorial')


@hiv.route('/trees/', methods=['GET', 'POST'])
def trees():
    if request.json:
        req = request.json
        tree = TreeModel(req['patient'], req['fragment']).get_newick_string()
        data = {'newick': tree}
        return jsonify(**data)

    form = TreeForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'all'
        fragment = 'F1'
    else:
        show_intro = False
        pname = form.patient.data
        fragment = form.fragment.data
        if not form.validate_on_submit():
            flash('Select one fragment and one patient!')

    data = {'pname': pname,
            'fragment': fragment,
            'name': pname+', '+fragment,
            'id': pname+'_'+fragment}

    return render_template('trees.html',
                           title='Phylogenetic trees',
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name='Phylogenetic trees',
                          )


@hiv.route(find_section(id='physio')['url'], methods=['GET', 'POST'])
def physio():
    if request.json:
        pname = request.json['patient']
        data = {'data': PhysioModel(pname).get_data()}
        return jsonify(**data)

    form = PatSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
    else:
        show_intro = False
        pname = form.patient.data
        if not form.validate_on_submit():
            flash('Select a patient!')

    data = {'pname': pname,
            'name': pname,
            'id': pname}

    return render_template('physio.html',
                           title='Viral load and CD4+ counts',
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name='Viral load and CD4+ counts',
                          )


@hiv.route('/n_templates/', methods=['POST'])
def n_templates():
    if request.json:
        pname = request.json['patient']
        data = {'data': NTemplatesModel(pname).get_data()}
        return jsonify(**data)

    else:
        abort(403)


@hiv.route('/divdiv/', methods=['GET', 'POST'])
def divdiv():
    if request.json:
        req = request.json
        pname = request.json['patient']
        fragment = request.json['fragment']
        data = {'data': DivdivModel(pname, fragment).get_data()}
        return jsonify(**data)

    form = PatFragSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
        fragment = 'F1'
    else:
        show_intro = False
        pname = form.patient.data
        fragment = form.fragment.data
        if not form.validate_on_submit():
            flash('Select a fragment and a patient!')

    data = {'pname': pname,
            'fragment': fragment,
            'name': pname+', '+fragment,
            'id': pname+'_'+fragment}

    return render_template('divdiv.html',
                           title='Divergence and diversity',
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name='Divergence and diversity',
                          )

@hiv.route(find_section(id='genome')['url'], methods=['GET', 'POST'])
def genomes():
    if request.json:
        pname = request.json['patient']
        data = {'data': GenomeModel(pname).get_data()}
        return jsonify(**data)

    section = find_section(id='genome')

    form = PatSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
    else:
        show_intro = False
        pname = form.patient.data
        if not form.validate_on_submit():
            flash('Select a patient!')

    data = {'pname': pname,
            'name': pname,
            'id': pname}

    return render_template('genome.html',
                           title=section['name'],
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name=section['name'],
                          )


@hiv.route('/coverage/', methods=['GET', 'POST'])
def coverage():
    if request.json:
        pname = request.json['patient']
        data = {'data': CoverageModel(pname).get_data()}
        return jsonify(**data)

    form = PatSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
    else:
        show_intro = False
        pname = form.patient.data
        if not form.validate_on_submit():
            flash('Select at least one patient!')

    data = {'pname': pname,
            'name': pname,
            'id': pname}

    return render_template('coverage.html',
                           title='Coverage',
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name='Coverage',
                          )


@hiv.route(find_section(id='divdiv_local')['url'], methods=['GET', 'POST'])
def divdiv_local():
    if request.json:
        req = request.json
        pname = req['patient']
        kwargs = {}
        for key in ('observables', 'itimes', 'roi'):
            if key in req:
                kwargs[key] = req[key]

        data = {'data': DivdivLocalModel(pname, **kwargs).get_data()}
        return jsonify(**data)

    section = find_section(id='divdiv_local')

    form = PatSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
    else:
        show_intro = False
        pname = form.patient.data
        if not form.validate_on_submit():
            flash('Select at least one patient!')

    data = {'pname': pname,
            'name': pname,
            'id': pname}

    return render_template('divdiv_local.html',
                           title=section['name'],
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name=section['name'],
                          )


@hiv.route(find_section(id='af')['url'], methods=['GET', 'POST'])
def allele_frequencies():
    if request.json:
        pname = request.json['patient']
        data = {'data': AlleleFrequencyModel(pname).get_data()}
        return jsonify(**data)

    section = find_section(id='af')

    form = PatSingleForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
    else:
        show_intro = False
        pname = form.patient.data
        if not form.validate_on_submit():
            flash('Select at least one patient!')

    data = {'pname': pname,
            'name': pname,
            'id': pname}

    return render_template('allele_frequencies.html',
                           title=section['name'],
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name=section['name'],
                          )


@hiv.route(find_section(id='sfs')['url'], methods=['GET', 'POST'])
def sfs():
    if request.json:
        data = {'data': SFSModel().get_data()}
        return jsonify(**data)

    section = find_section(id='sfs')

    # NOTE: we probably want to support several plots here
    show_intro = True

    return render_template('sfs.html',
                           title=section['name'],
                           show_intro=show_intro,
                           section_name=section['name'],
                          )


@hiv.route(find_section(id='prop')['url'], methods=['GET', 'POST'])
def propagators():
    if request.json:
        data = {'data': PropagatorModel().get_data()}
        return jsonify(**data)

    section = find_section(id='prop')

    # NOTE: we might want to support several plots here
    show_intro = True

    return render_template('propagator.html',
                           title=section['name'],
                           show_intro=show_intro,
                           section_name=section['name'],
                          )


@hiv.route(find_section(id='haplo')['url'], methods=['GET', 'POST'])
def haplotypes():

    form = LocalHaplotypeForm()
    formpc = PrecompiledHaplotypeForm()
    section = find_section(id='haplo')

    if request.method == 'GET':
        show_intro = True
        return render_template('haplotypes.html',
                               title=section['name'],
                               show_intro=show_intro,
                               form=form,
                               formpc=formpc,
                               section_name=section['name'],
                              )

    show_intro = False
    if (not formpc.validate_on_submit()) and (not form.validate_on_submit()):
        flash('Form incorrectly filled!')

        return render_template('haplotypes.html',
                       title=section['name'],
                       show_intro=show_intro,
                       form=form,
                       formpc=formpc,
                       section_name=section['name'],
                      )

    if formpc.validate_on_submit():
        pname = formpc.patient.data
        region = formpc.region.data
        return redirect('/download/haplotypes_'+pname+'_'+region+'.zip')

    elif form.validate_on_submit():
        pname = form.patient.data
        region = form.roi.region.data
        start = form.roi.start.data - 1 #Inclusive coordinates
        end = form.roi.end.data
        roi = (region, start, end)

        hm = LocalHaplotypeModel(pname, roi)

        # Get the data in a temporary folder/file
        tmp_folder = '/home/hivwholewebu1/tmp/'
        fn = hm.get_data(tmp_root_folder=tmp_folder)

        # Read it in...
        with open(fn, 'r') as f:
            fstr = f.read()

        # ...and delete it
        hm.clean_temporary_folders()

        # TODO: refine this to show a success page with a download link etc. (that
        # changes the policies on temporary files, storage use, etc., so watch out)
        response = make_response(fstr)
        response.headers["Content-Disposition"] = "attachment; filename=alignments.zip"
        return response


@hiv.route(find_section(id='consensi')['url'], methods=['GET', 'POST'])
def consensi():
    form = PatFragGWForm()
    section = find_section(id='consensi')

    if request.method == 'GET':
        show_intro = True
        return render_template('consensi.html',
                       title='Consensus sequences',
                       form=form,
                       show_intro=show_intro,
                       section_name=section['name'],
                      )


    if not form.validate_on_submit():
        show_intro = False
        flash('Select one fragment and one patient!')
        return render_template('consensi.html',
                       title='Consensus sequences',
                       form=form,
                       show_intro=show_intro,
                       section_name=section['name'],
                      )

    pname = form.patient.data
    fragment = form.fragment.data
    return redirect('/download/consensi_'+pname+'_'+fragment+'.fasta')


# Proxy view factory for static files, for download
# NOTE: this is slightly redundant as one can directly access the
# static files, but allows us for abstraction if needed, by switching
# off static file serving and changing
# TODO: fix robots.txt
@hiv.route('/download/<path:path>')
def data_proxy(path):
    sf = hiv.config['DATA_SUBFOLDER']

    fields = resplit('_|\.', path)
    dtype = fields[0]

    # Switch by data type
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

    abort(404)
