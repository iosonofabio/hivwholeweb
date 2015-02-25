from flask import Blueprint, render_template

people = Blueprint('people', __name__,
                     url_prefix='/people',
                     static_folder='static',
                     template_folder='templates')

publications = Blueprint('publications', __name__,
                     url_prefix='/publications',
                     static_folder='static',
                     template_folder='templates')

impressum = Blueprint('impressum', __name__,
                     url_prefix='/impressum',
                     static_folder='static',
                     template_folder='templates')

@people.route('/', methods=['GET'])
def index():
    return render_template('people.html',
                           title='People & Funding')

@publications.route('/', methods=['GET'])
def index():
    return render_template('publications.html',
                           title='Publications')

@impressum.route('/', methods=['GET'])
def index():
    return render_template('impress.html',
                           title='Impressum')
