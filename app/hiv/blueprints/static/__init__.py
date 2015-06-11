# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for static html pages (people, impressum, etc.).
'''
# Modules
from flask import Blueprint, render_template



# Blueprint
static = Blueprint('static', __name__,
                    url_prefix='',
                    static_folder='static',
                    template_folder='templates')



# Views
@static.route('/people/', methods=['GET'])
def people():
    return render_template('people.html',
                           title='People & Funding')


@static.route('/publications/', methods=['GET'])
def publications():
    return render_template('publications.html',
                           title='Publications')


@static.route('/impressum/', methods=['GET'])
def impressum():
    return render_template('impress.html',
                           title='Impressum')


@static.route('/controls/', methods=['GET'])
def controls():
    return render_template('controls.html',
                           title='Controls')
