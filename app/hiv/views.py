from flask import render_template, flash, redirect, request, jsonify
from . import hiv
from .forms import PatFragForm, PhysioForm, CoverageForm
from .models import TreeModel, PhysioModel, DivdivModel


@hiv.route('/')
@hiv.route('/index/')
def index():
    return render_template('index.html',
                           title='Home',
                           section_name='Home',
                          )


@hiv.route('/trees/', methods=['GET', 'POST'])
def trees():
    if request.json:
        req = request.json
        tree = TreeModel(req['patient'], req['fragment']).get_newick_string()
        data = {'newick': tree}
        return jsonify(**data)

    form = PatFragForm()
    if request.method == 'GET':
        show_intro = True
        pnames = ['all']
        fragments = ['F1']
    else:
        show_intro = False
        pnames = ['p'+str(i+1) for i in xrange(11)
                  if getattr(form, 'p'+str(i+1)).data]
        fragments = ['F'+str(i+1) for i in xrange(6)
                     if getattr(form, 'F'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one fragment and patient!')

    dicts = []
    for pname in pnames:
        for fragment in fragments:
            dicts.append({'pname': pname,
                          'fragment': fragment,
                          'name': pname+', '+fragment,
                          'id': pname+'_'+fragment})

    return render_template('trees.html',
                           title='Phylogenetic trees',
                           dicts=dicts,
                           form=form,
                           show_intro=show_intro,
                           section_name='Phylogenetic trees',
                          )


@hiv.route('/physiological/', methods=['GET', 'POST'])
def physio():
    if request.json:
        pname = request.json['patient']
        data = {'data': PhysioModel(pname).get_data()}
        return jsonify(**data)

    form = PhysioForm()
    if request.method == 'GET':
        show_intro = True
        pnames = ['p1']
    else:
        show_intro = False
        pnames = ['p'+str(i+1) for i in xrange(11)
                  if getattr(form, 'p'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one patient!')

    dicts = []
    for pname in pnames:
        pm = PhysioModel(pname)
        dicts.append({'pname': pname,
                      'name': pname,
                      'id': pname})

    return render_template('physio.html',
                           title='Viral load and CD4+ counts',
                           dicts=dicts,
                           form=form,
                           show_intro=show_intro,
                           section_name='Viral load and CD4+ counts',
                          )


@hiv.route('/divdiv/', methods=['GET', 'POST'])
def divdiv():
    if request.json:
        req = request.json
        data = {'data': DivdivModel(req['patient'], req['fragment']).get_data()}
        return jsonify(**data)

    form = PatFragForm()
    if request.method == 'GET':
        show_intro = True
        pnames = ['p1']
        fragments = ['F1']
    else:
        show_intro = False
        pnames = ['p'+str(i+1) for i in xrange(11)
                  if getattr(form, 'p'+str(i+1)).data]
        fragments = ['F'+str(i+1) for i in xrange(6)
                     if getattr(form, 'F'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one fragment and patient!')

    dicts = []
    for pname in pnames:
        for fragment in fragments:
            dicts.append({'pname': pname,
                          'fragment': fragment,
                          'name': pname+', '+fragment,
                          'id': pname+'_'+fragment})

    return render_template('divdiv.html',
                           title='Divergence and diversity',
                           dicts=dicts,
                           form=form,
                           show_intro=show_intro,
                           section_name='Divergence and diversity',
                          )


@hiv.route('/coverage/', methods=['GET', 'POST'])
def coverage():
    form = CoverageForm()
    if request.method == 'GET':
        fragments = ['F1']
    else:
        fragments = ['F'+str(i+1) for i in xrange(6) if getattr(form, 'F'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one fragment!')

    from .plots_hivwholeseq import coverage
    plot_dicts = []
    for fragment in fragments:
        plot_dict = {'fragment': fragment,
                     'chart': 'chart_'+fragment}

        data = coverage(fragment=fragment)
        plot_dict['data'] = data
        plot_dicts.append(plot_dict)

    return render_template('coverage.html',
                           title='Coverage',
                           plot_dicts=plot_dicts,
                           section_name='Coverage',
                           form=form,
                          )


@hiv.route('/allele_frequencies/')
def allele_frequencies():
    data = [[3, 3, 3],
            [10, 0, 3]]

    plot_dict = {'data': data,
                 'fragment': 'F0',
                 'chart': 'chart_F0'}

    return render_template('allele_frequencies.html',
                           title='Allele frequencies',
                           plot_dict=plot_dict,
                           section_name='Allele frequencies',
                          )

    pass
