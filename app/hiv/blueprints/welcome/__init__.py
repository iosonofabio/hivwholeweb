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
    
    # first 11 colors of d3 category20
    colors = { 'p1': '#1f77b4', 'p2': '#aec7e8', 'p3': '#ff7f0e', 'p4': '#ffbb78', 'p5': '#2ca02c',
               'p6': '#98df8a', 'p7': '#d62728', 'p8': '#ff9896', 'p9': '#9467bd', 'p10': '#c5b0d5', 'p11': '#8c564b' }
    
    return render_template('welcome.html',
                           patientTable=table,
                           colors=colors, 
                           region='V3',
                           title='Welcome page')


# FIXME for Bianca: this link should disappear, but it breaks APIs so delete it ASAP
@welcome.route('/welcome/', methods=['GET'])
def welcome_redirect():
    return redirect('/')
