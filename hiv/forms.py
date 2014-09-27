from flask.ext.wtf import Form
from wtforms import TextField
from wtforms.validators import Required


class TreeForm(Form):
    fragment = TextField('fragment', validators=[Required()])
