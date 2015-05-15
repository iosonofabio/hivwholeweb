# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/11/14
content:    Models for the hivwholeweb site.

            We are in a particular position, because we do not collect data
            from the users and store in a database, but we need to stream data.

            This makes our 'models' a bit odd.
'''
from flask_restful import Resource, abort



class Tree(Resource):
    def get(self, pname, region):
        from ....models import TreeModel
        try:
            return TreeModel(pname, region).get_json()
        except IOError:
            abort(404,
                  message="Patient and region {}, {} doesn't exist".format(pname, region))


class Physiological(Resource):
    def get(self, pname):
        from ....models import PhysioModel
        try:
            return PhysioModel(pname).get_data()
        except IOError:
            abort(404, message="Patient {} doesn't exist".format(pname))


class NTemplates(Resource):
    def get(self, pname):
        from ....models import NTemplatesModel
        try:
            return NTemplatesModel(pname).get_data()
        except IOError:
            abort(404, message="Patient {} doesn't exist".format(pname))


class ReferenceSequence(Resource):
    def get(self, pname, region):
        from ....models import GenomeModel
        try:
            return GenomeModel(pname, region).get_data()
        except IOError:
            abort(404,
                  message="Patient and region {}, {} doesn't exist".format(pname, region))


class Coverage(Resource):
    def get(self, pname):
        from ....models import CoverageModel
        try:
            return CoverageModel(pname).get_data()
        except IOError:
            abort(404, message="Patient {} doesn't exist".format(pname))


class AlleleFrequency(Resource):
    def get(self, pname):
        from ....models import AlleleFrequencyModel
        try:
            return AlleleFrequencyModel(pname).get_data()
        except IOError:
            abort(404, message="Patient {} doesn't exist".format(pname))


class Haplotypes(Resource):
    def get(self, pname, region):
        from ....models import HaplotypePrecompiledModel
        try:
            return HaplotypePrecompiledModel(pname, region).get_data(format='json')
        except IOError:
            abort(404,
                  message="Patient and region {}, {} doesn't exist".format(pname, region))



