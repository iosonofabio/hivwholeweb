from flask import Blueprint, render_template, abort

from ... import hiv



region = Blueprint('region', __name__,
                    url_prefix='/region',
                    static_folder='static',
                    template_folder='templates')


@region.route('/<regionname>/', methods=['GET'])
def index(regionname):
    if regionname not in hiv.config['REGIONS_SNP']:
        abort(404)

    pnames = ['p'+str(i) for i in xrange(1, 12)]

    return render_template('region.html',
                           region=regionname,
                           title='Region page')
