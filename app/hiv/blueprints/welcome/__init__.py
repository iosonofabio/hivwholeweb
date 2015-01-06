from flask import Blueprint, render_template


welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')


@welcome.route('/', methods=['GET'])
def index():
    return render_template('welcome.html',
                           title='Welcome page')
