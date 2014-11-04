from flask.ext.wtf import Form
from wtforms import BooleanField
from wtforms.validators import Required


class TreeForm(Form):
    F1 = BooleanField('F1', default=True)
    F2 = BooleanField('F3')
    F3 = BooleanField('F3')
    F4 = BooleanField('F4')
    F5 = BooleanField('F5')
    F6 = BooleanField('F6')


    def validate(self):
        if super(TreeForm, self).validate():
            return any((self.F1.data, self.F2.data, self.F3.data,
                        self.F4.data, self.F5.data, self.F6.data))
        else:
            return False


class PhysioForm(Form):
    p1 = BooleanField('p1', default=True)
    p2 = BooleanField('p3')
    p3 = BooleanField('p3')
    p4 = BooleanField('p4')
    p5 = BooleanField('p5')

    def validate(self):
        if super(PhysioForm, self).validate():
            return any((self.p1.data, self.p2.data, self.p3.data,
                        self.p4.data, self.p5.data))
        else:
            return False


class CoverageForm(Form):
    F1 = BooleanField('F1', default=True)
    F2 = BooleanField('F3')
    F3 = BooleanField('F3')
    F4 = BooleanField('F4')
    F5 = BooleanField('F5')
    F6 = BooleanField('F6')


    def validate(self):
        if super(CoverageForm, self).validate():
            return any((self.F1.data, self.F2.data, self.F3.data,
                        self.F4.data, self.F5.data, self.F6.data))
        else:
            return False

