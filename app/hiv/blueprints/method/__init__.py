# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/15
content:    Blueprint for the methods page.
'''
from re import split as resplit
from flask import Blueprint, render_template
from flask import (render_template, flash, redirect, request, jsonify,
                   make_response, abort)
from .forms import (LocalHaplotypeForm, TreeForm, ConsensiForm, RegionFragForm,
                    PatSingleForm, PatFragSingleForm, PrecompiledHaplotypeForm)
from .config import sections
from ... import hiv



method = Blueprint('method', __name__,
                    url_prefix='/method',
                    static_folder='static',
                    template_folder='templates')
hiv.config['BLUEPRINTS']['METHOD'] = {'SECTIONS': sections}


@method.route('/')
def index():
    return render_template('index.html',
                           title='Methods',
                           section_name='Method')


@method.route('/singleNucleotidePolymorphisms/', methods=['GET'])
def allele_frequencies():
    return render_template('allele_frequencies.html',
                           title='Single nucleotide polymorphisms',
                           section_name='Single nucleotide polymorphisms',
                          )


@method.route('/coverage/', methods=['GET'])
def coverage():
    return render_template('coverage.html',
                           title='Coverage',
                           section_name='Coverage',
                          )


@method.route('/physiological/', methods=['GET'])
def physio():
    return render_template('physio.html',
                           title='Viral load and CD4+ counts',
                           section_name='Viral load and CD4+ counts',
                          )



@method.route('/trees/', methods=['GET'])
def trees():
    return render_template('trees.html',
                           title='Phylogenetic trees',
                           section_name='Phylogenetic trees',
                          )


@method.route('/divdiv/', methods=['GET'])
def divdiv():
    return render_template('divdiv.html',
                           title='Divergence and diversity',
                           section_name='Divergence and diversity',
                          )


@method.route('/genomes/', methods=['GET'])
def genomes():
    return render_template('genome.html',
                           title='Genome sequences',
                           section_name='Genome sequences',
                          )


@method.route('/consensi/', methods=['GET', 'POST'])
def consensi():
    form = ConsensiForm()

    if request.method == 'GET':
        show_intro = True
        return render_template('consensi.html',
                       title='Consensus sequences',
                       section_name='Consensus sequences',
                       form=form,
                       show_intro=show_intro,
                      )


    if not form.validate_on_submit():
        show_intro = False
        flash('Select one region and one patient!')
        return render_template('consensi.html',
                       title='Consensus sequences',
                       section_name='Consensus sequences',
                       form=form,
                       show_intro=show_intro,
                      )

    pname = form.patient.data
    region = form.region.data
    return redirect('/download/consensi/'+pname+'_'+region+'.fasta')


@method.route('/divdivLocal/', methods=['GET'])
def divdiv_local():
    return render_template('divdiv_local.html',
                           title='Local divergence and diversity',
                           section_name='Local divergence and diversity',
                          )


@method.route('/haplotypes/', methods=['GET', 'POST'])
def haplotypes():
    form = LocalHaplotypeForm()
    formpc = PrecompiledHaplotypeForm()

    if request.method == 'GET':
        show_intro = True
        return render_template('haplotypes.html',
                               title='Haplotypes',
                               section_name='Haplotypes',
                               show_intro=show_intro,
                               form=form,
                               formpc=formpc,
                              )

    show_intro = False
    if (not formpc.validate_on_submit()) and (not form.validate_on_submit()):
        flash('Form incorrectly filled!')

        return render_template('haplotypes.html',
                               title='Haplotypes',
                               section_name='Haplotypes',
                               show_intro=show_intro,
                               form=form,
                               formpc=formpc,
                              )

    if formpc.validate_on_submit():
        pname = formpc.patient.data
        region = formpc.region.data
        return redirect('/download/haplotypes/'+pname+'/'+region+'.fasta')

    elif form.validate_on_submit():
        from ...models import LocalHaplotypeModel

        # NOTE: we offer only genomewide HXB2 coordinates
        region = 'genomewide'
        pname = form.patient.data
        start = form.roi.start.data - 1 #Inclusive coordinates
        end = form.roi.end.data
        roi = (region, start, end)

        hm = LocalHaplotypeModel(pname, roi)
        try:
            hm.translate_coordinates()
        except ValueError:
            flash('No PCR fragment covers such region, please select a narrower region')

            return render_template('haplotypes.html',
                                   title='Haplotypes',
                                   section_name='Haplotypes',
                                   show_intro=show_intro,
                                   form=form,
                                   formpc=formpc,
                                  )


        # Get the data from a temporary folder + file
        # NOTE: the temp folders are cleaned regularly at a timeout specified
        # in the config file, so no need to to that here.
        fn = hm.get_data()

        # TODO: refine this to show a success page with a download link etc. (that
        # changes the policies on temporary files, storage use, etc., so watch out)
        with open(fn, 'r') as f:
            fstr = f.read()
        response = make_response(fstr)
        response.headers["Content-Disposition"] = ("attachment; filename="+
                                                   "haplotypes_"+pname+
                                                   "_"+str(start+1)+
                                                   "_"+str(end)+
                                                   ".fasta")
        return response
