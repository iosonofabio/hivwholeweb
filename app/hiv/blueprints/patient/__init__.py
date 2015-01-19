from flask import Blueprint, render_template, abort


patient = Blueprint('patient', __name__,
                    url_prefix='/patient',
                    static_folder='static',
                    template_folder='templates')


@patient.route('/p<int:patient_number>/', methods=['GET'])
def index(patient_number):
    if patient_number not in range(1, 12):
        abort(404)

    pname = 'p'+str(patient_number)

    return render_template('patient.html',
                           pname=pname,
                           title='Patient page')
