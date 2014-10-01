from flask import render_template, flash, redirect, request
from hiv import hiv
from forms import TreeForm, CoverageForm
from .backbone import sections


@hiv.route('/')
@hiv.route('/index/')
def index():
    return render_template('index.html',
                           title='Home',
                           sections=sections,
                           section_name='Home',
                          )


@hiv.route('/trees/', methods=['GET', 'POST'])
def trees():
    form = TreeForm()
    if request.method == 'GET':
        fragments = ['F1']
    else:
        fragments = ['F'+str(i+1) for i in xrange(6) if getattr(form, 'F'+str(i+1)).data]
        if not form.validate_on_submit():
            flash('Select at least one fragment!')

    trees = [{'url': '/static/images/tree_consensi_'+fragment+'.png'}
             for fragment in fragments]

    return render_template('trees.html',
                           title='Phylogenetic trees',
                           trees=trees,
                           form=form,
                           sections=sections,
                           section_name='Phylogenetic trees',
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
                           sections=sections,
                           section_name='Coverage',
                           form=form,
                          )
