from flask import Blueprint, render_template, abort
from ...models import PatientTableModel, SampleTableModel

patient = Blueprint('patient', __name__,
                    url_prefix='/patient',
                    static_folder='static',
                    template_folder='templates')


@patient.route('/p<int:patient_number>/', methods=['GET'])
def index(patient_number):
    if patient_number not in range(1, 12):
        abort(404)

    pname = 'p'+str(patient_number)
    pnames = ['p'+str(i) for i in xrange(1, 12)]
    table = PatientTableModel().get_table()
    tree_regions = ['V3', 'psi', 'p17', 'RRE']
    divdiv_regions = ['V3', 'psi', 'p17', 'RRE']

    sample_table = SampleTableModel(pname).get_table()

    return render_template('patient.html',
                           pname=pname,
                           pnames=pnames,
                           patientTable=table, 
                           sampleTable=sample_table,
                           treeRegions=tree_regions,
                           divDivRegions=divdiv_regions,
                           title='Patient page')
