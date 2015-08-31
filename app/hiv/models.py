# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/11/14
content:    Models for the hivwholeweb site.

            We are in a particular position, because we do not collect data
            from the users and store in a database, but we need to stream data.

            This makes our 'models' a bit odd.
'''
import os
from . import hiv
data_folder_full = hiv.static_folder+'/'+hiv.config['DATA_SUBFOLDER']+'/'
data_folder_short = '/'+os.path.basename(hiv.static_folder)+'/'+hiv.config['DATA_SUBFOLDER']+'/'
data_folder = [data_folder_short, data_folder_full]



class TreeModel(object):
    def __init__(self, pname, region):
        self.pname = pname
        if 'minor' not in region:
            self.region = region
            self.minor = False
        else:
            self.region = region[:-len('_minor')]
            self.minor = True


    def get_filename(self, full=True, format='newick'):
        if not self.minor:
            fn = 'consensi_tree_'+self.pname+'_'+self.region+'.'+format
        else:
            fn = 'haplotype_tree_'+self.pname+'_'+self.region+'.'+format

        return data_folder[full]+'trees/'+fn

    
    def get_newick_string(self):
        fn = self.get_filename(format='newick')
        with open(fn, 'r') as f:
            tree = f.read().rstrip('\n')
            # NOTE: Adding an artificial root creates a long branch.
            root_dist = tree.split(':')[-1][:-1]
            if float(root_dist) > 0.01:
                tree = tree[:tree.rfind(':')]+'0.001;'
        return tree


    def get_json_filename(self, full=False):
        return self.get_filename(full=full, format='json')


    def get_json(self):
        import json
        with open(self.get_json_filename(full=True), 'r') as f:
            return json.load(f)


class PhysioModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_viral_load_filename(self, full=True, format='tsv'):
        fn = 'viral_load_'+self.pname+'.'+format
        return data_folder[full]+'physiological/'+fn


    def get_cell_count_filename(self, full=True, format='tsv'):
        fn = 'cell_count_'+self.pname+'.'+format
        return data_folder[full]+'physiological/'+fn


    def get_physio_filename(self, kind, *args, **kwargs):
        if kind in ['viralLoad', 'vl']:
            return self.get_viral_load_filename(*args, **kwargs)
        elif kind in ['cellCount', 'cc']:
            return self.get_cell_count_filename(*args, **kwargs)
        else:
            ValueError('Data kind not understood: '+str(kind))


    def get_data(self, obs=['viral load', 'CD4+ cell count'], full_headers=False, zipped=True):
        import numpy as np
        vl = np.loadtxt(self.get_viral_load_filename(), skiprows=1)
        cc = np.loadtxt(self.get_cell_count_filename(), skiprows=1)

        if not zipped:
            vl = vl.T
            cc = cc.T
        vl = map(list, vl)
        cc = map(list, cc)
        data = {'time unit': '[days since infection]'}
        if 'viral load' in obs:
            data['viral load'] = vl
            data['viral load unit'] = '[copies/ml]'
        if 'CD4+ cell count' in obs:
            data['CD4+ cell count'] = cc
            data['CD4+ cell count unit'] = '[cells/ul]'

        return data


class NTemplatesModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_ntemplates_filename(self, full=True):
        fn = 'ntemplates_'+self.pname+'.tsv'
        return data_folder[full]+'physiological/'+fn


    def get_data(self):
        import numpy as np
        nt = np.loadtxt(self.get_ntemplates_filename(), skiprows=1)

        nt = map(list, nt)

        data = {'ntemplates': nt}
        return data


class GenomeModel(object):
    def __init__(self, psname, region='genomewide'):
        self.psname = psname
        self.region = region

        # Check whether we need a patient reference or a sample consensus
        from . import hiv
        if psname in hiv.config['PATIENTS']:
            self.type = 'patient reference'
        elif psname in hiv.config['REFERENCES']:
            self.type = 'external reference'
        else:
            self.type = 'consensus'


    def get_reference_filename(self, full=True, format='gb'):
        fn = 'reference_'+self.psname+'_'+self.region+'.'+format
        return data_folder[full]+'sequences/'+fn


    def get_consensus_filename(self, full=True, format='fasta'):
        fn = 'consensus_'+self.psname+'_'+self.region+'.'+format
        return data_folder[full]+'sequences/'+fn


    def get_genome_filename(self, *args, **kwargs):
        if self.type in ['patient reference', 'external reference']:
            return self.get_reference_filename(*args, **kwargs)
        elif self.type == 'consensus':
            return self.get_consensus_filename(*args, **kwargs)
        else:
            raise ValueError('Genome type not found')


    def get_genome(self, format='gb'):
        from Bio import SeqIO
        seq = SeqIO.read(self.get_genome_filename(format=format), format)

        if format == 'gb':
            from .cookbook import correct_genbank_features_load
            correct_genbank_features_load(seq)

        return seq


    def get_data(self):
        seq = self.get_genome()
        
        features = []
        frame_start = 0
        for fea in seq.features:
            loc = [[lp.nofuzzy_start, lp.nofuzzy_end]
                      for lp in fea.location.parts]

            f = {'name': fea.id,
                 'type': fea.type,
                 'location': loc}
            features.append(f)

            # Get all frames relative to gag (which is frame 1)
            if fea.id == 'gag':
                frame_start = loc[0][0]

        data = {'len': len(seq),
                'id': seq.id,
                'name': seq.name,
                'description': seq.description,
                'framestart': frame_start,
                'seq': ''.join(seq),
                'features': features,
               }

        return data


class DivdivModel(object):
    def __init__(self, pname, fragment):
        self.pname = pname
        self.fragment = fragment


    def get_divergence_filename(self, full=True):
        fn = 'divergence_'+self.pname+'_'+self.fragment+'.tsv'
        return data_folder[full]+'divdiv/'+fn


    def get_diversity_filename(self, full=True):
        fn = 'diversity_'+self.pname+'_'+self.fragment+'.tsv'
        return data_folder[full]+'divdiv/'+fn


    def get_data(self, zipped=True):
        import numpy as np
        dg = np.loadtxt(self.get_divergence_filename(), skiprows=1)
        ds = np.loadtxt(self.get_diversity_filename(), skiprows=1)

        if not zipped:
            dg = dg.T
            ds = ds.T

        dg = map(list, dg)
        ds = map(list, ds)

        data = {'dg': dg,
                'ds': ds}
        return data


class CoverageModel(object):
    def __init__(self, samplename, region='genomewide'):
        self.samplename = samplename
        self.region = region


    def get_coverage_filename(self, full=True):
        fn = 'coverage_'+self.samplename+'_'+self.region+'.npz'
        return data_folder[full]+'coverage/'+fn


    def get_data(self):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_coverage_filename())
        cov = npz['cov']

        # Saturate to one, the JS client would take forever to do that
        cov = list(np.maximum(0.8, cov))

        data = {'cov': cov,
                'region': self.region,
                'sample': self.samplename,
               }
        return data


class CoverageTrajectoryModel(object):
    def __init__(self, pname, region='genomewide'):
        self.pname = pname
        self.region = region


    def get_coverage_filename(self, full=True):
        fn = 'coverage_'+self.pname+'_'+self.region+'.npz'
        return data_folder[full]+'coverage/'+fn


    def get_data(self):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_coverage_filename())
        times = npz['times']
        cov = npz['cov']

        # Saturate to one, the JS client would take forever to do that
        cov = [[t, list(np.maximum(0.8, co))] for (t, co) in izip(times, cov)]

        data = {'cov': cov}
        return data


class AlleleFrequencyModel(object):
    def __init__(self, samplename, region='genomewide'):
        self.samplename = samplename
        self.region = region


    def get_allele_counts_filename(self, full=True, format='npz'):
        fn = 'allele_counts_'+self.samplename+'_'+self.region+'.'+format
        return data_folder[full]+'single_nucleotide_variants/'+fn


    def get_data(self, cov_min=100, fmt='full'):

        if fmt != 'full':
            raise ValueError('Allele frequencies are only in full format')

        from itertools import izip
        import numpy as np
        npz = np.load(self.get_allele_counts_filename())
        ac = npz['ac']
        alpha = list(npz['alpha'])

        cov = ac.sum(axis=0)
        ind = cov >= cov_min
        af = 1.0 * ac / np.maximum(cov, 1)
        af[:, -ind] = np.nan

        af = map(list, af)
        L = len(af[0])

        # Info message
        msg = ('Each list corresponds to a nucleotide, according to the alphabet. '+
               'Low coverage (<'+str(cov_min)+') positions are masked as NaNs.')

        data = {'SNPs': af,
                'length': L,
                'samplename': self.samplename,
                'region': self.region,
                'alphabet': alpha,
                'coverage_min': cov_min,
                'info': msg,
               }
        return data


class AlleleFrequencyTrajectoryModel(object):
    def __init__(self, pname, region='genomewide'):
        self.pname = pname
        self.region = region


    def get_allele_counts_filename(self, full=True, format='npz'):
        fn = 'allele_counts_'+self.pname+'_'+self.region+'.'+format
        return data_folder[full]+'single_nucleotide_variants/'+fn


    def get_data(self, cov_min=50, af_min=1e-1, fmt='full'):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_allele_counts_filename())
        times = npz['times']
        act = npz['act']
        alpha = npz['alpha']

        aft = []
        times_set = set()
        for pos in xrange(act.shape[2]):
            act_pos = act[:, :, pos]
            cov_pos = act_pos.sum(axis=1)

            # Ind are the good time points for this position. The earliest of
            # those sets the ancestral allele (not plotted)
            ind = cov_pos >= cov_min
            if ind.sum() < 2:
                continue

            act_pos = act_pos[ind]
            cov_pos = cov_pos[ind]
            aft_pos = (1.0 * act_pos.T / cov_pos).T

            if fmt == 'sparse':
                # Ignore ancestral alleles
                aft_pos[:, act_pos[0].argmax()] = 0
                cand = ((aft_pos > af_min).sum(axis=0) >= 2).nonzero()[0]

            elif fmt == 'full':
                cand = xrange(6)

            else:
                raise ValueError('Allele frequency trajectories format: sparse|full')

            for ai in cand:
                af = aft_pos[:, ai]
                # FIXME: this should be done by plots only, but it's here for
                # performance reasons
                if fmt == 'sparse':
                    af[af < 3e-3] = 1e-3

                datum = [(t, float('{:1.3f}'.format(aftmp)))
                         for (t, aftmp) in izip(times[ind], af)]

                aft.append([pos, alpha[ai], datum])

                times_set |= set(times[ind])

            
        data = {'data': aft,
                'tmax': times.max(),
                'times': sorted(times_set),
                'len': act.shape[2],
                'format': ('List of [position (bp in patient reference), '+
                           'nucleotide, '+
                           'List of pairs [time from infection (days), SNP frequency]]')
               }
        return data


class DivdivLocalModel(object):
    def __init__(self, pname,
                 observables=['dg', 'ds'],
                 itimes=None, roi=None):
        self.pname = pname
        self.observables = observables
        self.itimes = itimes
        self.roi = roi

    def get_divergence_filename(self, full=True):
        fn = 'divergence_trajectory_'+self.pname+'_genomewide.npz'
        return data_folder[full]+'divdiv/'+fn


    def get_diversity_filename(self, full=True):
        fn = 'diversity_trajectory_'+self.pname+'_genomewide.npz'
        return data_folder[full]+'divdiv/'+fn


    def trim_to_roi(self, data):
        if self.roi is None:
            return data

        # FIXME: implement this
        if self.roi[0] != 'genomewide':
            raise ValueError('Not implemented!')

        # NOTE: we use double inclusive coordinates
        return data[:, int(self.roi[1]): int(self.roi[2]) + 1]


    def get_data(self, dx=20):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_divergence_filename())
        times = npz['times']
        dg = self.trim_to_roi(npz['dg'])[:, ::dx]
        block_length = npz['block_length'][0]
        len_times = len(times)

        # NOTE: negative numbers are missing coverage, what do we do with those?
        dg = [[t, list(np.maximum(0.0001, d))] for (t, d) in izip(times, dg)]

        npz = np.load(self.get_diversity_filename())
        times = npz['times']
        ds = self.trim_to_roi(npz['ds'])[:, ::dx]
        len_times = max(len_times, len(times))

        # NOTE: we assume div and div have the same block length
        ds = [[t, list(np.maximum(0.0001, d))] for (t, d) in izip(times, ds)]

        data = {'block_length': block_length,
                'dx': dx,
                'L': int(npz['L']),
                'nTimes': len_times}

        # FIXME: this is buggy, we don't get the itime from the raw data
        if self.itimes is not None:
            dg = [dg[i] for i in self.itimes]
            ds = [ds[i] for i in self.itimes]

        if 'dg' in self.observables:
            data['dg'] = dg
        if 'ds' in self.observables:
            data['ds'] = ds

        return data


class LocalHaplotypeModel(object):
    def __init__(self, pname, roi):
        self.pname = pname
        self.roi = roi
        self.dirnames = []


    def get_bam_filename(self, i_time, fragment, full=True):
        '''Get the BAM filename of a sample and a fragment'''
        return (data_folder[full]+'reads/'+self.pname+'_sample_'+
                str(i_time+1)+'_'+fragment+'.bam')


    def get_timeline_filename(self, full=True):
        '''Get the filename of the timeline'''
        return data_folder[full]+'time_points/time_points_'+self.pname+'.tsv'


    def get_coordinate_map_filename(self, full=True, refname='HXB2', format='tsv'):
        '''Get the filename of the coordinate map'''
        fn = 'coordinate_map_'+self.pname+'_'+refname+'_genomewide.'+format
        return data_folder[full]+'coordinate_maps/'+fn


    def get_local_haplotype_filename(self, full=True):
        '''Get the filename of a temporary file with the haplotype data'''
        import random
        
        from . import hiv
        tmp_root_folder = hiv.config['TMP_ROOT_FOLDER']

        dirname = tmp_root_folder+str(random.randint(0, 10000))+'/'
        while os.path.isdir(dirname):
            dirname = tmp_root_folder+str(random.randint(0, 10000))+'/'

        self.dirnames.append(dirname)
        os.mkdir(dirname)

        fn = dirname+'haplotypes_'+self.pname+'_'+'_'.join(map(str, self.roi))+'.fasta'
        return fn


    def translate_coordinates(self):
        '''Translate HXB2 coordinate into patient coordinates'''
        import numpy as np
        # NOTE: the map is always genomewide??
        mapco = dict(np.loadtxt(self.get_coordinate_map_filename(full=True),
                                dtype=int))

        # In case the coordinate is missing, extend the region
        # Translate start position
        pos = self.roi[1]
        pmin = min(mapco.iterkeys())
        while pos not in mapco:
            pos -= 1
            if pos < pmin:
                raise ValueError('Initial position outside of sequenced area.')
        start = mapco[pos]

        # Translate end position
        pos = self.roi[2]
        pmax = max(mapco.iterkeys())
        while pos not in mapco:
            pos += 1
            if pos > pmax:
                raise ValueError('Final position outside of sequenced area.')
        end = mapco[pos]

        # Find the fragment if available
        from Bio import SeqIO
        refseq = SeqIO.read(GenomeModel(self.pname).get_reference_filename(), 'gb')
        for fea in refseq.features:
            feaname = fea.qualifiers['note'][0]
            if feaname not in ['F'+str(i) for i in xrange(1, 7)]:
                continue

            fea_start = fea.location.nofuzzy_start
            fea_end = fea.location.nofuzzy_end

            if (fea_start <= start) and (fea_end >= end):
                start -= fea_start
                end -= fea_start
                break
        else:
            raise ValueError('No fragment fully covering the haplotype')

        self.roi = (feaname, start, end)


    def clean_temporary_folders(self, timeout='default'):
        import os, shutil

        if timeout == 'default':
            from . import hiv
            timeout = hiv.config['TIMEOUT_TMP']

        if timeout is not None:
            import time
            now = time.time()

        if isinstance(timeout, basestring):
            unit = timeout[-1]
            timeout = float(timeout[:-1])
            if timeout in ['m', 'h', 'd', 'w']:
                timeout *= 60
            if timeout in ['h', 'd', 'w']:
                timeout *= 60
            if timeout in ['d', 'w']:
                timeout *= 24
            if timeout in ['w']:
                timeout *= 7

        for dirname in self.dirnames:
            if os.path.isdir(dirname):
                if timeout is not None:
                    timem = os.path.getmtime(dirname)
                    if (now - timem) < timeout:
                        continue

                shutil.rmtree(dirname)
        self.dirnames = []


    @staticmethod
    def store_alignment(ali, filename):
        '''Save alignments to file'''
        from Bio import AlignIO
        file_formats = {'stk': 'stockholm',
                        'fasta': 'fasta',
                        'phy': 'phylip-relaxed'}
    
        foldername = os.path.dirname(filename)
        if not os.path.isdir(foldername):
            raise IOError('Destination folder for file save not found')
    
        if os.path.isfile(filename):
            raise IOError('Destination file already exists on file system')
    
        file_format = filename.split('.')[-1]
        if file_format in file_formats:
            AlignIO.write(ali, filename, file_format)
    
        else:
            raise ValueError('File format not recognized')



    def get_data(self, clean=True):
        '''Get the data'''
        import numpy as np
        from .analysis.get_local_haplotypes import (get_local_haplotypes_formatted,
                                                    align_sequences)

        if clean:
            self.clean_temporary_folders()

        times = np.loadtxt(self.get_timeline_filename())

        fragment, start, end = self.roi 

        ali = []
        for i_time in xrange(len(times)):
            bamfilename = self.get_bam_filename(i_time, fragment)
            if not os.path.isfile(bamfilename):
                continue

            label_id = 'days_'+str(int(times[i_time]))+'_'
            label_description = 'days since infection: '+str(int(times[i_time]))+', '

            seqs = get_local_haplotypes_formatted(bamfilename, start, end,
                                                  label_id=label_id,
                                                  label_description=label_description,
                                                  VERBOSE=2)

            ali.extend(seqs)

        ali = align_sequences(ali)

        fn_out = self.get_local_haplotype_filename()
        self.store_alignment(ali, fn_out)

        return fn_out


class HaplotypePrecompiledModel(object):
    def __init__(self, samplename, region):
        self.samplename = samplename
        self.region = region


    def get_haplotype_filename(self, full=True, format='fasta'):
        fn = 'haplotype_alignment_'+self.samplename+'_'+self.region+'.'+format
        return data_folder[full]+'alignments/'+fn


    def get_data(self, format='json'):
        from Bio import AlignIO
        fn = self.get_haplotype_filename(full=True, format='fasta')
        ali = AlignIO.read(fn, 'fasta')

        if format == 'biopython':
            return ali

        elif format == 'json':
            alij = [{'name': s.name,
                     'description': s.description,
                     'sequence': str(s.seq),
                     'frequency [%]': float(s.name.split('_')[-1][:-1]),
                    }
                   for s in ali]

            return {'alignment': alij,
                    'samplename': self.samplename,
                    'region': self.region,
                   }

        else:
            raise ValueError('Format not understood')


class HaplotypePrecompiledTrajectoryModel(object):
    def __init__(self, pname, region):
        self.pname = pname
        self.region = region


    def get_haplotype_filename(self, full=True, format='fasta'):
        fn = 'haplotype_alignment_'+self.pname+'_'+self.region+'.'+format
        return data_folder[full]+'alignments/'+fn


    def get_data(self, format='json'):
        from Bio import AlignIO
        fn = self.get_haplotype_filename(full=True, format='fasta')
        ali = AlignIO.read(fn, 'fasta')

        if format == 'biopython':
            return ali

        elif format == 'json':
            alij = [{'name': s.name,
                     'description': s.description,
                     'sequence': str(s.seq),
                     'days since infection': int(s.name.split('_')[1]),
                     'frequency [%]': float(s.name.split('_')[-1][:-1]),
                    }
                   for s in ali]
            return alij

        else:
            raise ValueError('Format not understood')



class PatientTableModel(object):

    def get_table_filename(self, full=True, format='tsv'):
        '''Get the filename of the patient table'''
        fn = 'patients.'+format
        return data_folder[full]+'tables/'+fn


    def get_table(self):
        '''Read the table from file'''
        fn = self.get_table_filename(full=True)

        table = []
        fieldinds = []
        with open(fn, 'r') as f:
            headerfields = f.readline().rstrip('\n').split('\t')
            for field in headerfields:
                fieldinds.append(field) 

            # Second header line, empty
            f.readline()

            for line in f:
                tline = {}
                fieldsline = line.rstrip('\n').split('\t')
                if len(fieldsline) == 0:
                    continue
                tline = {fieldinds[ifi]: field for ifi, field in enumerate(fieldsline)}
                table.append(tline)

        return table


class SampleTableModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_table_filename(self, full=True, format='tsv', quantitative=False):
        '''Get the filename of the samples table'''
        fn = 'samples'
        if quantitative:
            fn = fn+'_quantitative'
        else:
            fn = fn+'_qualitative'
        fn = fn+'_'+self.pname+'.'+format
        return data_folder[full]+'tables/'+fn


    def get_table(self, fields=None):
        '''Read the table from file'''
        fn = self.get_table_filename(full=True)

        table = []
        fieldinds = []
        with open(fn, 'r') as f:
            headerfields = f.readline().rstrip('\n').split('\t')
            for field in headerfields:
                if field == 'days since infection':
                    field = 'time'
                elif field == 'RNA templates':
                    field = 'RNA'

                if (fields is None) or (field in fields):
                    fieldinds.append(field)

            for line in f:
                tline = {}
                fieldsline = line.rstrip('\n').split('\t')
                if len(fieldsline) == 0:
                    continue

                for ifi, field in enumerate(fieldsline):
                    if ifi >= len(fieldinds):
                        break

                    # Qualitative description of PCR performance
                    if fieldinds[ifi] in ['F'+str(i) for i in xrange(1, 7)]:
                        # Use reduced alphabet for website
                        if field in ['end lower', 'low diversity']:
                            datafmt = 'low'
                        else:
                            datafmt = field

                    # Provide RNA templates with 2 significant digits
                    elif fieldinds[ifi] == 'RNA':
                        from math import floor, log10
                        fmtfun = lambda x: int(round(float(x), 1 - int(floor(log10(float(x))))))
                        field = float(field)
                        if field < 1:
                            datafmt = '-'
                        else:
                            datafmt = fmtfun(field)
                            if datafmt < 1000:
                                datafmt = '{:2d}'.format(datafmt)
                            else:
                                datafmt = '{:2d} 000'.format(datafmt // 1000)                        

                    # Any other field
                    else:
                        datafmt = '{:1.0f}'.format(float(field))
                    tline[fieldinds[ifi]] = datafmt
                table.append(tline)

        return table


class ReadsTableModel(SampleTableModel):
    def get_reads_filename(self, i, fragment, full=True, format='bam'):
        '''Get the filename of the table with the presence of reads files'''
        fn = fragment+'.'+format
        fn = data_folder[full]+'reads/'+self.pname+'_sample_'+str(i)+'_'+fn
        return fn


    def get_table(self):
        import os

        table = super(ReadsTableModel, self).get_table(fields=['time'])
        for it, datum in enumerate(table, 1):
            for fragment in ['F'+str(i) for i in xrange(1, 7)]:
                datum[fragment] = os.path.isfile(self.get_reads_filename(it, fragment))

        return table

class CocountsTableModel(SampleTableModel):
    def get_cocounts_filename(self, i, fragment, full=True, format='bam'):
        '''Get the filename of the table with the presence of reads files'''
        fn = fragment+'.'+format
        fn = data_folder[full]+'pair_nucleotide_variants/cocounts_'+self.pname+'_sample_'+str(i)+'_'+fn
        return fn


    def get_table(self):
        import os

        table = super(CocountsTableModel, self).get_table(fields=['time'])
        for it, datum in enumerate(table, 1):
            for fragment in ['F'+str(i) for i in xrange(1, 7)]:
                datum[fragment] = os.path.isfile(self.get_cocounts_filename(it, fragment))

        return table


class HLAModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_filename(self, full=True, format='tsv'):
        '''Get the filename of the samples table'''
        fn = 'HLA'
        fn = fn+'_'+self.pname+'.'+format
        return data_folder[full]+'hla/'+fn

