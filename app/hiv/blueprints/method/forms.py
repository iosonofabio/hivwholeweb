# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/12/14
content:    Forms for the method blueprint.
'''
# Modules
from flask.ext.wtf import Form
from wtforms import BooleanField, IntegerField, SelectField, FormField
from wtforms.validators import Required

from ... import hiv



# Globals
_regions = ('V3', 'PR', 'psi', 'vpr', 'vpu', 'p15', 'p17',
            'p6', 'p2', 'p1', 'p7', 'RRE')
patients = hiv.config['PATIENTS']



# Classes
class RoiForm(Form):
    # It only accepts one region, we only use genomewide HXB2 coordinates
    start = IntegerField('From')
    end = IntegerField('To')

    # TODO: make clean meassages for failures
    def validate(self):
        if not super(RoiForm, self).validate():
            return False

        start = self.start.data
        end = self.end.data

        if (start < 1) or (end < 1) or (end < start):
            return False

        if (end - start > 600):
            return False

        return True


class LocalHaplotypeForm(Form):
    patient = SelectField('Patient',
                          choices=[[p, p] for p in patients])
    roi = FormField(RoiForm)


class PrecompiledHaplotypeForm(Form):
    patient = SelectField('Patient',
                          choices=[[p, p] for p in patients])
    region = SelectField('Region',
                         choices=[[reg] *2 for reg in _regions])


class TreeForm(Form):
    patient = SelectField('Patient',
                          choices=[['all', 'all']] + [[p, p] for p in patients]),
    region = SelectField('Region',
                          choices=([['F'+str(i)] * 2 for i in xrange(1, 7)] +
                                   [[reg] * 2 for reg in _regions] +
                                   [[reg+'_minor', reg+' minor'] for reg in _regions]),
                        )


class PatSingleForm(Form):
    patient = SelectField('Patient',
                          choices=[[reg] *2 for reg in _regions])


class PatFragSingleForm(Form):
    patient = SelectField('Patient',
                          choices=[[p, p] for p in patients])
    fragment = SelectField('Fragment',
                           choices=[['F'+str(i)] * 2 for i in xrange(1, 7)])


class ConsensiForm(Form):
    patient = SelectField('Patient',
                          choices=[[p, p] for p in patients])
    region = SelectField('Region',
                          choices=([['genomewide', 'genomewide']] +
                                   [[reg] *2 for reg in _regions] +
                                   [['F'+str(i)] * 2 for i in xrange(1, 7)]))


class RegionFragForm(Form):
    patient = SelectField('Patient',
                          choices=[[p, p] for p in patients])
    region = SelectField('Fragment',
                          choices=([['F'+str(i)] * 2 for i in xrange(1, 7)] +
                                   [[reg] *2 for reg in _regions]))

