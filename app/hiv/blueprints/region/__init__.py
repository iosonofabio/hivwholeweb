# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the region page.
'''
# Modules
from flask import Blueprint, render_template, abort



# Blueprint
region = Blueprint('region', __name__,
                    url_prefix='/region',
                    static_folder='static',
                    template_folder='templates')



# Views
@region.route('/<regionname>/', methods=['GET'])
def index(regionname):
    from ... import hiv
    if regionname not in hiv.config['REGIONS_SNP']:
        abort(404)

    return render_template('region.html',
                           region=regionname,
                           title='Region page')
