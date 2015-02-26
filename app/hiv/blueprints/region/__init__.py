from flask import Blueprint, render_template, abort
from . import backbone



region = Blueprint('region', __name__,
                    url_prefix='/region',
                    static_folder='static',
                    template_folder='templates')


@region.route('/<regionname>/', methods=['GET'])
def index(regionname):
    if regionname not in backbone.regions:
        abort(404)

    return render_template('region.html',
                           region=regionname,
                           regions=backbone.regions, 
                           title='Region page')
