# vim: fdm=indent
'''
author:     Fabio Zanini
date:       30/03/15
content:    Forms for the patient blueprint.
'''
# Modules
from flask.ext.wtf import Form
from wtforms import BooleanField, IntegerField, SelectField, FormField
from wtforms.validators import Required

from ..method.forms import (_regions, RoiForm)


# Classes
class ConsensiForm(Form):
    region = SelectField('Region',
                         id='ConsRegField',
                         choices=([['genomewide', 'genomewide']] +
                                  [[reg] *2 for reg in _regions] +
                                  [['F'+str(i)] * 2 for i in xrange(1, 7)]))


class PrecompiledHaplotypeForm(Form):
    region = SelectField('Region',
                         choices=[[reg] *2 for reg in _regions])

