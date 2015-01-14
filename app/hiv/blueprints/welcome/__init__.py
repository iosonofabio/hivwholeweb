from flask import Blueprint, render_template


welcome = Blueprint('welcome', __name__,
                    url_prefix='/welcome',
                    static_folder='static',
                    template_folder='templates')


@welcome.route('/', methods=['GET'])
def index():

    fn = __path__[0]+'/../../static/data/patients.csv'
    
    import numpy as np
    pnames = np.loadtxt(fn, usecols=[0], dtype='S3', skiprows=2)
    subtypes = np.loadtxt(fn, usecols=[1], dtype='S2', skiprows=2)
    num_samples = np.loadtxt(fn, usecols=[2], dtype='S6', skiprows=2)
    first_day = np.loadtxt(fn, usecols=[3], dtype='S5', skiprows=2)
    last_day = np.loadtxt(fn, usecols=[4], dtype='S4', skiprows=2)
    
    table = zip(pnames, subtypes, num_samples, first_day, last_day)

    return render_template('welcome.html',
                           patientTable=table,
                           title='Welcome page')
