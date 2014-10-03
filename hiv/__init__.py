from flask import Flask
import backbone

hiv = Flask(__name__)
hiv.config.from_object('config')
hiv.config['SECTIONS'] = backbone.sections

from hiv import views
