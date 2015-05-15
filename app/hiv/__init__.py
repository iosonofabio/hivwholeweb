#def apptest(environ, start_response):
#    status = '200 OK'
#    output = 'Hello World from within the hiv app, after all that!'
#
#    response_headers = [('Content-type', 'text/plain'),
#                        ('Content-Length', str(len(output)))]
#
#    start_response(status, response_headers)
#
#    return [output]
from flask import Flask


# TODO: fix static folder, this requires adapting HTML templates
hiv = Flask(__name__)
hiv.config.from_object('config')

# welcome page is the index
from .blueprints.welcome import welcome
hiv.register_blueprint(welcome, url_prefix='')

from .blueprints.tutorial import tutorial
hiv.register_blueprint(tutorial, url_prefix='/tutorial')

from .blueprints.patient import patient
hiv.register_blueprint(patient)

from .blueprints.region import region
hiv.register_blueprint(region)

from .blueprints.data import data
hiv.register_blueprint(data)

from .blueprints.method import method, backbone
hiv.register_blueprint(method)
hiv.config['BLUEPRINTS']['METHOD'] = {'SECTIONS': backbone.sections}

from .blueprints.static import static
hiv.register_blueprint(static)

# RESTful API
from .blueprints.api import api_bp
hiv.register_blueprint(api_bp)


# Import the global views
from . import views
