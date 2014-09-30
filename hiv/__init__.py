from flask import Flask

hiv = Flask(__name__)
from bokeh.server.app import bokeh_app
hiv.register_blueprint(bokeh_app)
hiv.config.from_object('config')

from hiv import views
