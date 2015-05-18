from re import split as resplit
from flask import Blueprint, render_template
from flask import (render_template, flash, redirect, request, jsonify,
                   make_response, abort)
from .forms import (LocalHaplotypeForm, TreeForm, ConsensiForm, RegionFragForm,
                    PatSingleForm, PatFragSingleForm, PrecompiledHaplotypeForm)
from ...models import (TreeModel, PhysioModel, DivdivModel, CoverageModel,
                     GenomeModel, AlleleFrequencyModel,
                     NTemplatesModel, DivdivLocalModel, LocalHaplotypeModel)
from .backbone import find_section



method = Blueprint('method', __name__,
                    url_prefix='/method',
                    static_folder='static',
                    template_folder='templates')


@method.route('/')
def index():
    return render_template('index.html',
                           title='Methods',
                           section_name='Method')


@method.route(find_section(id='af')['url'], methods=['GET', 'POST'])
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

@method.route(find_section(id='cov')['url'], methods=['GET', 'POST'])
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


@method.route(find_section(id='physio')['url'], methods=['GET', 'POST'])
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



@method.route(find_section(id='tree')['url'], methods=['GET', 'POST'])
def trees():
    if request.json:
        req = request.json
        tree = TreeModel(req['patient'], req['region'])

        # Newick returns the string, JSON returns a JSON
        data = {}
        if ('treeFormat' in req) and (req['treeFormat'].lower() == 'json'):
            data['tree'] = tree.get_json_filename()
        else:
            data['tree'] = tree.get_newick_string()

        return jsonify(**data)

    form = TreeForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'all'
        region = 'p17'
    else:
        show_intro = False
        pname = form.patient.data
        region = form.region.data
        if not form.validate_on_submit():
            flash('Select one region and one patient!')

    data = {'pname': pname,
            'region': region,
            'name': pname+', '+region,
            'id': pname+'_'+region}

    if 'minor' in region:
        plot_title = region[:-len('_minor')]+", with minor haplotypes"
    else:
        plot_title = region+", consensus sequences"

    return render_template('trees.html',
                           title='Phylogenetic trees',
                           data=data,
                           form=form,
                           plotTitle=plot_title,
                           show_intro=show_intro,
                           section_name='Phylogenetic trees',
                          )


@method.route(find_section(id='divdiv')['url'], methods=['GET', 'POST'])
def divdiv():
    if request.json:
        req = request.json
        pname = request.json['patient']
        region = request.json['region']
        data = {'data': DivdivModel(pname, region).get_data()}
        return jsonify(**data)

    form = RegionFragForm()
    if request.method == 'GET':
        show_intro = True
        pname = 'p1'
        region = 'F1'
    else:
        show_intro = False
        pname = form.patient.data
        region = form.region.data
        if not form.validate_on_submit():
            flash('Select a region and a patient!')

    data = {'pname': pname,
            'region': region,
            'name': pname+', '+region,
            'id': pname+'_'+region}

    return render_template('divdiv.html',
                           title='Divergence and diversity',
                           data=data,
                           form=form,
                           show_intro=show_intro,
                           section_name='Divergence and diversity',
                          )


@method.route(find_section(id='genome')['url'], methods=['GET', 'POST'])
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


@method.route(find_section(id='consensi')['url'], methods=['GET', 'POST'])
def consensi():
    form = ConsensiForm()
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
        flash('Select one region and one patient!')
        return render_template('consensi.html',
                       title='Consensus sequences',
                       form=form,
                       show_intro=show_intro,
                       section_name=section['name'],
                      )

    pname = form.patient.data
    region = form.region.data
    return redirect('/download/consensi_'+pname+'_'+region+'.fasta')


@method.route(find_section(id='divdiv_local')['url'], methods=['GET', 'POST'])
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


@method.route(find_section(id='haplo')['url'], methods=['GET', 'POST'])
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
        return redirect('/download/haplotypes_'+pname+'_'+region+'.fasta')

    elif form.validate_on_submit():
        # NOTE: we offer only genomewide HXB2 coordinates
        region = 'genomewide'
        pname = form.patient.data
        start = form.roi.start.data - 1 #Inclusive coordinates
        end = form.roi.end.data
        roi = (region, start, end)

        hm = LocalHaplotypeModel(pname, roi)
        try:
            hm.translate_coordinates()
        except ValueError:
            flash('No PCR fragment covers such region, please select a narrower region')

            return render_template('haplotypes.html',
                           title=section['name'],
                           show_intro=show_intro,
                           form=form,
                           formpc=formpc,
                           section_name=section['name'],
                          )


        # Get the data from a temporary folder + file
        # NOTE: the temp folders are cleaned regularly at a timeout specified
        # in the config file, so no need to to that here.
        fn = hm.get_data()

        # TODO: refine this to show a success page with a download link etc. (that
        # changes the policies on temporary files, storage use, etc., so watch out)
        with open(fn, 'r') as f:
            fstr = f.read()
        response = make_response(fstr)
        response.headers["Content-Disposition"] = ("attachment; filename="+
                                                   "haplotypes_"+pname+
                                                   "_"+str(start+1)+
                                                   "_"+str(end)+
                                                   ".zip")
        return response


@method.route('/n_templates/', methods=['POST'])
def n_templates():
    if request.json:
        pname = request.json['patient']
        data = {'data': NTemplatesModel(pname).get_data()}
        return jsonify(**data)

    else:
        abort(403)


