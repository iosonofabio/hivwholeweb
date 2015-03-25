from flask import Blueprint, render_template, abort


data = Blueprint('data', __name__,
                 url_prefix='/data',
                 static_folder='static',
                 template_folder='templates')


@data.route('/', methods=['GET'])
def index():
    return render_template('data.html',
                          )
