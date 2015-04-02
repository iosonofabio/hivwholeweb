from flask import (Blueprint, render_template, abort, request,
                   redirect, flash, make_response)
from ...models import (PatientTableModel, SampleTableModel,
                       LocalHaplotypeModel)
from .forms import RoiForm, PrecompiledHaplotypeForm, ConsensiForm

patient = Blueprint('patient', __name__,
                    url_prefix='/patient',
                    static_folder='static',
                    template_folder='templates')


@patient.route('/p<int:patient_number>/', methods=['GET', 'POST'])
def index(patient_number):
    if patient_number not in range(1, 12):
        abort(404)
    pname = 'p'+str(patient_number)

    formco = ConsensiForm()
    formht = RoiForm()
    formpc = PrecompiledHaplotypeForm()

    if request.method == 'POST':
        if ((not formpc.validate_on_submit()) and
            (not formht.validate_on_submit()) and
            (not formco.validate_on_submit())):
            flash('Form incorrectly filled!')

        elif formco.validate_on_submit():
            region = formco.region.data
            return redirect('/download/consensi_'+pname+'_'+region+'.fasta')

        elif formpc.validate_on_submit():
            region = formpc.region.data
            return redirect('/download/haplotypes_'+pname+'_'+region+'.fasta')

        elif formht.validate_on_submit():
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
                                                           ".zip")
                return response


    pnames = ['p'+str(i) for i in xrange(1, 12)]
    table = PatientTableModel().get_table()
    tree_regions = ['V3', 'V3_minor', 'p17', 'p17_minor', 'psi', 'RRE']
    divdiv_regions = ['V3', 'psi', 'p17', 'RRE']

    sample_table = SampleTableModel(pname).get_table()


    return render_template('patient.html',
                           pname=pname,
                           pnames=pnames,
                           patientTable=table, 
                           sampleTable=sample_table,
                           treeRegions=tree_regions,
                           divDivRegions=divdiv_regions,
                           formCons=formco,
                           formht=formht,
                           formpc=formpc,
                           title='Patient page')
