# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the RESTful API.

            The API is basically used to GET the data in JSON format for further
            web-friendly manipulation. No PUT, DELETE, HEAD methods are supported.
            The API gives access to a number of data including:
                - viral load
                - CD4+ cell count
                - reference sequences
                - phylogenetic trees
                - coverage
                - single nucleotide polymorphism counts
                - haplotype counts
'''
# Modules
from flask import Blueprint
from flask_restful import Api, Resource
from .resources import *



# API blueprint
api_bp = Blueprint('api', __name__,
                   url_prefix='/api',
                   template_folder='templates')

api = Api(api_bp)

api.add_resource(Tree,
                 '/data/tree/<pname>/<region>',
                )

api.add_resource(Physiological,
                 '/data/physiological/<pname>',
                )

api.add_resource(ViralLoad,
                 '/data/viralLoad/<pname>',
                )

api.add_resource(CellCount,
                 '/data/cellCount/<pname>',
                )

api.add_resource(NTemplates,
                 '/data/numberTemplates/<pname>',
                )

api.add_resource(ReferenceSequence,
                 '/data/referenceSequence/<pname>',
                 '/data/referenceSequence/<pname>/<region>',
                )

api.add_resource(DivDiv,
                 '/data/divdiv/<pname>/<region>',
                )

api.add_resource(DivDivSliding,
                 '/data/divdivSliding/<pname>',
                )

api.add_resource(Coverage,
                 '/data/coverage/<psname>',
                 '/data/coverage/<psname>/<region>',
                )

# FIXME: this should be way more flexible (only polymorphic sites are returned now)
api.add_resource(AlleleFrequencies,
                 '/data/singleNucleotidePolymorphisms/<psname>',
                 '/data/singleNucleotidePolymorphisms/<psname>/<region>',
                )

api.add_resource(Haplotypes,
                 '/data/haplotypes/<psname>/<region>',
                )
