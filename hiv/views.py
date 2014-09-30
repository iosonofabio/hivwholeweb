from flask import render_template, flash, redirect, request
from hiv import hiv
from forms import TreeForm
from .backbone import sections


@hiv.route('/')
@hiv.route('/index/')
def index():
    return render_template('index.html',
                           title='Home',
                           sections=sections,
                          )


@hiv.route('/test/')
def test():
    return render_template('test.html',
                           title='BokehJS test',
                           sections=sections,
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
                          )


@hiv.route('/coverage/')
def coverage():
    fragment = 'F1'

    bokeh_dict = {'has_data': False}
    from test_bokeh import coverage
    try:
        (js, div) = coverage(fragment=fragment)
        bokeh_dict['js'] = js
        bokeh_dict['div'] = div
        bokeh_dict['has_data'] = True
    except IOError:
        flash('Could not reach the source data.')


    return render_template('coverage.html',
                           title='Coverage',
                           bokeh_dict=bokeh_dict,
                           sections=sections,
                          )
