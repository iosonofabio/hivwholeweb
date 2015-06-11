# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the tutorial.
'''
# Modules
from flask import Blueprint, render_template



# Blueprint
tutorial = Blueprint('tutorial', __name__,
                     url_prefix='/tutorial',
                     static_folder='static',
                     template_folder='templates')



# Views
@tutorial.route('/', methods=['GET'])
def index():
    return render_template('tutorial.html',
                           title='Tutorial')
