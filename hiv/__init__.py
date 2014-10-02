from flask import Flask

hiv = Flask(__name__)
hiv.config.from_object('config')

from hiv import views
