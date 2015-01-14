from flask import Blueprint, render_template
from ...models import PatientTableModel


welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')


@welcome.route('/', methods=['GET'])
def index():
    table = PatientTableModel().get_table()

    return render_template('welcome.html',
                           patientTable=table,
                           title='Welcome page')
