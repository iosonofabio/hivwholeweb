# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the patient page.
'''
# Modules
from flask import (Blueprint, render_template, abort, request,
                   redirect, flash, make_response)
from ...models import (PatientTableModel, SampleTableModel,
                       LocalHaplotypeModel)
from .forms import RoiForm, PrecompiledHaplotypeForm, ConsensiForm



# Blueprint
patient = Blueprint('patient', __name__,
                    url_prefix='/patient',
                    static_folder='static',
                    template_folder='templates')



# Views
@patient.route('/<pname>/', methods=['GET', 'POST'])
def index(pname):
    from ... import hiv
    if pname not in hiv.config['PATIENTS']:
        abort(404)

    formco = ConsensiForm()
    formht = RoiForm()
    formpc = PrecompiledHaplotypeForm()

    if request.method == 'POST':

        # We got a form submit request, figure out which one
        subm_id = request.form['formBtn']

        if subm_id == 'submitCo':
            if formco.validate_on_submit():
                region = formco.region.data
                return redirect('/download/consensi_'+pname+'_'+region+'.fasta')
            else:
                flash('Form for consensus sequences not filled correctly!')

        elif subm_id == 'submitPc':
            if formpc.validate_on_submit():
                region = formpc.region.data
                return redirect('/download/haplotypes_'+pname+'_'+region+'.fasta')
            else:
                flash('Form for precompiled haplotype alignments not filled correctly!')

        elif subm_id == 'submitHt':
            if formht.validate_on_submit():
                # NOTE: we offer only genomewide HXB2 coordinates
                region = 'genomewide'
                start = formht.start.data - 1 #Inclusive coordinates
                end = formht.end.data
                roi = (region, start, end)

                hm = LocalHaplotypeModel(pname, roi)
                try:
                    hm.translate_coordinates()
                    cont = True
                except ValueError:
                    flash('No PCR fragment covers such region, please select a narrower region')
                    cont = False

                if cont:

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
                                                               ".fasta")
                    return response

            else:
                flash('Form for new haplotype alignments not filled correctly!')

        else:
            # It is not a known form, default to GET response
            pass


    patient_table = PatientTableModel().get_table()
    sample_table = SampleTableModel(pname).get_table()

    return render_template('patient.html',
                           pname=pname,
                           patientTable=patient_table, 
                           sampleTable=sample_table,
                           formco=formco,
                           formpc=formpc,
                           formht=formht,
                           title='Patient page')
