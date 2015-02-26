from flask import Blueprint, render_template

# People page
static = Blueprint('static', __name__,
                    url_prefix='',
                    static_folder='static',
                    template_folder='templates')

@static.route('/people/', methods=['GET'])
def people():
    return render_template('people.html',
                           title='People & Funding')


@static.route('/publications/', methods=['GET'])
def publications():
    return render_template('publications.html',
                           title='Publications')


@static.route('/impressum/', methods=['GET'])
def impressum():
    return render_template('impress.html',
                           title='Impressum')
