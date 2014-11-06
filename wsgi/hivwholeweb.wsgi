# Taken directly from Flask documentation
# http://flask.pocoo.org/docs/0.10/deploying/mod_wsgi/
import sys
sys.path.insert(0, '/var/www/hivwholeweb/app')
from hiv import hiv as application


## Test for now
#def application(environ, start_response):
#    status = '200 OK'
#    output = 'Hello World!'
#    output += ' '
#    output += ' '.join(sys.path)
#
#    response_headers = [('Content-type', 'text/plain'),
#                        ('Content-Length', str(len(output)))]
#
#    start_response(status, response_headers)
#
#    return [output]
