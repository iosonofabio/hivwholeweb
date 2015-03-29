from flask import (Blueprint, render_template, abort, flash,
                   redirect, request, jsonify)


data = Blueprint('data', __name__,
                 url_prefix='/data',
                 static_folder='static',
                 template_folder='templates')


@data.route('/', methods=['GET', 'POST'])
def index():
    # Haplotype forms
    from ..method.forms import LocalHaplotypeForm, PrecompiledHaplotypeForm
    form = LocalHaplotypeForm()
    formpc = PrecompiledHaplotypeForm()

    if request.method == 'POST':
        if ((not formpc.validate_on_submit()) and
            (not form.validate_on_submit())):
            flash('Haplotype form incorrectly filled!')

        elif formpc.validate_on_submit():
            pname = formpc.patient.data
            region = formpc.region.data
            return redirect('/download/haplotypes_'+pname+'_'+region+'.zip')

        elif form.validate_on_submit():
            flash('UNDER CONSTRUCTION!')


    # NOTE: we are using the proxy view for static downloads as a screen
    # to the real locations on the server. A nice side effect is that we need
    # not instantiate data models within this view.
    # That factory is quite poorly written though, and would need some love.
    from ...models import PatientTableModel
    data = []
    table = PatientTableModel().get_table()
    d = lambda x: '/download/'+x
    for row in table:
        pname = row['Patient']
        n_samples = row['# samples']

        data.append({'pname': pname,
                     'nsamples': n_samples,
                     'vl': d('physio_vl_'+pname+'.dat'),
                     'cd4': d('physio_cc_'+pname+'.dat'),
                     'refgb': d('genome_'+pname+'.gb'),
                     'reffa': d('genome_'+pname+'.fasta'),
                     'act': d('act_'+pname+'.zip'),
                    })

    return render_template('data.html',
                           data=data,
                           form=form,
                           formpc=formpc,
                           title='Downloads')
