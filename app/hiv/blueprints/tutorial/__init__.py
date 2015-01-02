from flask import Blueprint, render_template



tutorial = Blueprint('tutorial', __name__,
                     url_prefix='/tutorial',
                     static_folder='static',
                     template_folder='templates')



@tutorial.route('/', methods=['GET'])
def index():
    return render_template('tutorial.html',
                           title='Tutorial')
