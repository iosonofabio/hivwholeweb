from flask import Blueprint, render_template, redirect
from ...models import PatientTableModel


welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')


@welcome.route('/', methods=['GET'])
@welcome.route('/index/', methods=['GET'])
def index():
    table = PatientTableModel().get_table()

    return render_template('welcome.html',
                           patientTable=table,
                           region='V3',
                           title='Welcome page')


# FIXME for Bianca: this link should disappear, but it breaks APIs so delete it ASAP
@welcome.route('/welcome/', methods=['GET'])
def welcome_redirect():
    return redirect('/')
