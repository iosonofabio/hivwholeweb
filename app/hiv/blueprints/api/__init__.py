from flask import Blueprint

api = Blueprint('api', __name__,
                url_prefix='/api',
                static_folder='static',
                template_folder='templates')


# Views
@api.route('/test/', methods=['GET'])
def test():
    return 'API test'

