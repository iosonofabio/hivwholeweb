from flask import Blueprint, render_template


patient = Blueprint('patient', __name__,
                    url_prefix='/patient',
                    static_folder='static',
                    template_folder='templates')


@patient.route('/', methods=['GET'])
def index():
    return render_template('patient.html',
                           title='Patient page')
