from flask import Blueprint, render_template, abort

# FIXME: add all regions
regions = ('V3', 'p17', '24', 'PR', 'RT', 'vif', 'vpu', 'nef')



region = Blueprint('region', __name__,
                    url_prefix='/region',
                    static_folder='static',
                    template_folder='templates')


@region.route('/<regionname>', methods=['GET'])
def index(regionname):
    if regionname not in regions:
        abort(404)

    return render_template('region.html',
                           region=regionname,
                           title='Region page')
