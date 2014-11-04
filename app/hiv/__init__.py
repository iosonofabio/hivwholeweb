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
import backbone


hiv = Flask(__name__)
hiv.config.from_object('config')
hiv.config['SECTIONS'] = backbone.sections

from . import views

