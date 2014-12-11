from flask.ext.wtf import Form
from wtforms import BooleanField, IntegerField, SelectField, FormField
from wtforms.validators import Required


class PatFragForm(Form):
    F1 = BooleanField('F1', default=True)
    F2 = BooleanField('F3')
    F3 = BooleanField('F3')
    F4 = BooleanField('F4')
    F5 = BooleanField('F5')
    F6 = BooleanField('F6')

    p1 = BooleanField('p1', default=True)
    p2 = BooleanField('p3')
    p3 = BooleanField('p3')
    p4 = BooleanField('p4')
    p5 = BooleanField('p5')
    p6 = BooleanField('p6')
    p7 = BooleanField('p7')
    p8 = BooleanField('p8')
    p9 = BooleanField('p9')
    p10 = BooleanField('p10')
    p11 = BooleanField('p11')

    def validate(self):
        if super(PatFragForm, self).validate():
            has_frag = any((self.F1.data, self.F2.data, self.F3.data,
                            self.F4.data, self.F5.data, self.F6.data))
            has_pat = any((self.p1.data, self.p2.data, self.p3.data,
                           self.p4.data, self.p5.data, self.p6.data,
                           self.p7.data, self.p8.data, self.p9.data,
                           self.p10.data, self.p11.data))
            return has_frag and has_pat
        else:
            return False


class PatForm(Form):
    p1 = BooleanField('p1', default=True)
    p2 = BooleanField('p3')
    p3 = BooleanField('p3')
    p4 = BooleanField('p4')
    p5 = BooleanField('p5')
    p6 = BooleanField('p6')
    p7 = BooleanField('p7')
    p8 = BooleanField('p8')
    p9 = BooleanField('p9')
    p10 = BooleanField('p10')
    p11 = BooleanField('p11')

    def validate(self):
        if super(PatForm, self).validate():
            return any((self.p1.data, self.p2.data, self.p3.data,
                        self.p4.data, self.p5.data, self.p6.data,
                        self.p7.data, self.p8.data, self.p9.data,
                        self.p10.data, self.p11.data))
        else:
            return False


class RoiForm(Form):
    # FIXME: allow all fragments, including HXB2 and genomewide
    region = SelectField('Region',
                           choices=[['F'+str(i)] * 2 for i in xrange(1, 7)])
    start = IntegerField('Start')
    end = IntegerField('End')

    # TODO: make dynamic validator that accepts different fragment
    # arguments, e.g. genomewide or HXB2
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
                          choices=[['p'+str(i)] * 2 for i in xrange(1, 12)])
    roi = FormField(RoiForm)


class PrecompiledHaplotypeForm(Form):
    patient = SelectField('Patient',
                          choices=[['p'+str(i)] * 2 for i in xrange(1, 12)])
    region = SelectField('Region',
                         choices=[('V3', 'V3'),
                                  ('PR', 'PR'),
                                  ('psi', 'psi')])


class TreeForm(Form):
    _regions = ('V3', 'PR', 'psi')
    patient = SelectField('Patient',
                          choices=[['all'] * 2] + [['p'+str(i)] * 2 for i in xrange(1, 12)])
    region = SelectField('Region',
                          choices=([['F'+str(i)] * 2 for i in xrange(1, 7)] +
                                   [[reg] *2 for reg in _regions] +
                                   [[reg+'_minor', reg+' minor'] for reg in _regions]),
                        )


class PatSingleForm(Form):
    patient = SelectField('Patient',
                          choices=[['p'+str(i)] * 2 for i in xrange(1, 12)])


class PatFragSingleForm(Form):
    patient = SelectField('Patient',
                          choices=[['p'+str(i)] * 2 for i in xrange(1, 12)])
    fragment = SelectField('Fragment',
                           choices=[['F'+str(i)] * 2 for i in xrange(1, 7)])


class ConsensiForm(Form):
    _regions = ('V3', 'PR', 'psi')
    patient = SelectField('Patient',
                          choices=[['p'+str(i)] * 2 for i in xrange(1, 12)])
    fragment = SelectField('Fragment',
                           choices=([['genomewide', 'genomewide']] +
                                    [[reg] *2 for reg in _regions] +
                                    [['F'+str(i)] * 2 for i in xrange(1, 7)]))
