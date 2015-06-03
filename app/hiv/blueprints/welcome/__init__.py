from flask import Blueprint, render_template
from ...models import PatientTableModel


# BLUEPRINT
welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')


# VIEWS
@welcome.route('/', methods=['GET'])
@welcome.route('/index/', methods=['GET'])
@welcome.route('/welcome/', methods=['GET'])
def index():
    table = PatientTableModel().get_table()    
    return render_template('welcome.html',
                           patientTable=table,
                           region='p17',
                           title='Welcome page',
                          )
