# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the welcome page.
'''
# Modules
from flask import Blueprint, render_template
from ...models import PatientTableModel



# Blueprint
welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')



# Views
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
