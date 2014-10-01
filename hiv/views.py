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


@hiv.route('/test/')
def test():
    return render_template('test.html',
                           title='BokehJS test',
                           sections=sections,
                           section_name='Test Bokeh',
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

    from test_bokeh import coverage
    bokeh_dicts = []
    for fragment in fragments:
        bokeh_dict = {'has_data': False}
        try:
            (js, div) = coverage(fragment=fragment)
            bokeh_dict['js'] = js
            bokeh_dict['div'] = div
            bokeh_dict['has_data'] = True
        except IOError:
            flash('Could not reach the source data for fragment '+fragment+'.')
        bokeh_dicts.append(bokeh_dict)

    return render_template('coverage.html',
                           title='Coverage',
                           bokeh_dicts=bokeh_dicts,
                           sections=sections,
                           section_name='Coverage',
                           form=form,
                          )


@ hiv.route('/test_csv.csv')
def test_csv():
    with open('hiv/test_csv.csv', 'r') as f:
        text = f.read()
    return text
