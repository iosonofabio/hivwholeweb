from flask import render_template, flash, redirect
from hiv import hiv
from forms import TreeForm


@hiv.route('/')
@hiv.route('/index/')
def index():
    user = { 'nickname': 'Fabio' } # FIXME: fake user
    sections = [{'name': 'Phylogenetic trees', 'url': '/trees/'},
                {'name': 'Coverage', 'url': '/coverage/'},
                {'name': 'Allele frequencies'},
               ]
    return render_template('index.html',
                           title='Home',
                           user=user,
                           sections=sections,
                          )

@hiv.route('/trees/', methods=['GET', 'POST'])
def trees(fragment='F1'):

    form = TreeForm()
    if form.validate_on_submit():
        fragment = form.fragment.data
    else:
        # FIXME: distinguish GET from POST
        #flash('Fragment not found: '+form.fragment.data)
        pass

    trees = [{'url': '/static/images/tree_consensi_'+fragment+'.png'},
            ]
    return render_template('trees.html',
                           title='Phylogenetic trees',
                           trees=trees,
                           form=form)


@hiv.route('/coverage/')
def coverage():
    fragment = 'F1'

    from test_bokeh import coverage
    (js, div) = coverage(fragment=fragment)
    bokeh_dict = {'js': js,
                  'div': div}

    return render_template('coverage.html',
                           title='Coverage',
                           bokeh_dict=bokeh_dict)
