from flask import Blueprint
from flask_restful import Api, Resource


api_bp = Blueprint('api', __name__,
                   url_prefix='/api',
                   template_folder='templates')

api = Api(api_bp)


# Resources
from .resources import *

# Data
api.add_resource(Tree, '/data/tree/<pname>/<region>')

api.add_resource(Physiological,
                 '/data/physiological/<pname>',
                )

api.add_resource(ViralLoad,
                 '/data/viralLoad/<pname>',
                )

api.add_resource(CellCount,
                 '/data/cellCount/<pname>',
                )

api.add_resource(NTemplates, '/data/numberTemplates/<pname>')

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

api.add_resource(Haplotypes, '/data/haplotypes/<psname>/<region>')

# JS for plots
# FIXME: this becomes a CDN for embedding plots into other people's websites.
# Do we want this?
