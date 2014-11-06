from flask import render_template, flash, redirect, request, jsonify
from . import hiv
from .forms import TreeForm, PhysioForm, CoverageForm
from .models import TreeModel


@hiv.route('/')
@hiv.route('/index/')
def index():
    return render_template('index.html',
                           title='Home',
                           section_name='Home',
                          )


@hiv.route('/trees/', methods=['GET', 'POST'])
def trees():
    form = TreeForm()
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
            flash('Select at least one fragment!')

    trees = []
    for pname in pnames:
        for fragment in fragments:
            trees.append({'pname': pname,
                          'fragment': fragment,
                          'name': pname+', '+fragment,
                          'id': pname+'_'+fragment})

    return render_template('trees.html',
                           title='Phylogenetic trees',
                           trees=trees,
                           form=form,
                           show_intro=show_intro,
                           section_name='Phylogenetic trees',
                          )


@hiv.route('/treedata/', methods=['POST'])
def tree_json():
    treedict = request.json
    fragment = treedict['fragment']
    pname = treedict['patient']

    tree = TreeModel(pname, fragment).get_newick_string()
    treedict = {'newick': tree}
    return jsonify(**treedict)


@hiv.route('/physiological/', methods=['GET', 'POST'])
def physio():
    form = PhysioForm()
    if request.method == 'GET':
        pnames = ['p1']
    else:
        pnames = ['p'+str(i+1) for i in xrange(5) if getattr(form, 'p'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one patient!')

    dicts = [{'url': '/static/images/physiological/'+pname+'.png'}
             for pname in pnames]

    return render_template('physio.html',
                           title='Viral load and CD4+ counts',
                           dicts=dicts,
                           form=form,
                           section_name='Viral load and CD4+ counts',
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
