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
data_folder_short = '/static/data/'
data_folder_full = os.path.dirname(os.path.abspath(__file__))+data_folder_short
data_folder = [data_folder_short, data_folder_full]
tmp_folder_short = '/static/tmp/'
tmp_folder_full = os.path.dirname(os.path.abspath(__file__))+tmp_folder_short
tmp_folder = [tmp_folder_short, tmp_folder_full]



class TreeModel(object):
    def __init__(self, pname, region):
        self.pname = pname
        self.region = region


    def get_filename(self, full=True, format='newick'):
        fn = 'consensi_tree_'+self.pname+'_'+self.region+'.'+format
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


    def get_json_filename(self):
        return self.get_filename(full=False, format='json')


class PhysioModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_viral_load_filename(self, full=True):
        fn = 'viral_load_'+self.pname+'.dat'
        return data_folder[full]+'physiological/'+fn


    def get_cell_count_filename(self, full=True):
        fn = 'cell_count_'+self.pname+'.dat'
        return data_folder[full]+'physiological/'+fn


    def get_data(self):
        import numpy as np
        vl = np.loadtxt(self.get_viral_load_filename(), skiprows=1)
        cc = np.loadtxt(self.get_cell_count_filename(), skiprows=1)

        vl = map(list, vl)
        cc = map(list, cc)

        data = {'vl': vl,
                'cc': cc}
        return data


class NTemplatesModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_ntemplates_filename(self, full=True):
        fn = 'ntemplates_'+self.pname+'.dat'
        return data_folder[full]+'physiological/'+fn


    def get_data(self):
        import numpy as np
        nt = np.loadtxt(self.get_ntemplates_filename(), skiprows=1)

        nt = map(list, nt)

        data = {'ntemplates': nt}
        return data


class GenomeModel(object):
    def __init__(self, pname, region='genomewide', type='patient reference'):
        self.pname = pname
        self.region = region
        self.type = type


    def get_reference_filename(self, full=True, format='gb'):
        fn = 'reference_'+self.pname+'_'+self.region+'.'+format
        return data_folder[full]+'sequences/'+fn


    def get_genome_filename(self, *args, **kwargs):
        if self.type == 'patient reference':
            return self.get_reference_filename(*args, **kwargs)
        
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
        fn = 'divergence_'+self.pname+'_'+self.fragment+'.dat'
        return data_folder[full]+'divdiv/'+fn


    def get_diversity_filename(self, full=True):
        fn = 'diversity_'+self.pname+'_'+self.fragment+'.dat'
        return data_folder[full]+'divdiv/'+fn


    def get_data(self):
        import numpy as np
        dg = np.loadtxt(self.get_divergence_filename(), skiprows=1)
        ds = np.loadtxt(self.get_diversity_filename(), skiprows=1)

        dg = map(list, dg)
        ds = map(list, ds)

        data = {'dg': dg,
                'ds': ds}
        return data


class CoverageModel(object):
    def __init__(self, pname):
        self.pname = pname


    def get_coverage_filename(self, full=True):
        fn = 'coverage_'+self.pname+'_genomewide.npz'
        return data_folder[full]+'one_site/'+fn


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
    def __init__(self, pname):
        self.pname = pname


    def get_allele_counts_filename(self, full=True, format='npz'):
        fn = 'allele_counts_'+self.pname+'_genomewide.'+format
        return data_folder[full]+'one_site/'+fn


    def get_data(self, cov_min=100, af_min=1e-1):
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
            ind = cov_pos >= cov_min
            if ind.sum() < 2:
                continue

            aft_pos = (1.0 * act_pos[ind].T / cov_pos[ind]).T

            # Ignore ancestral alleles
            aft_pos[:, act_pos[0].argmax()] = 0

            cand = ((aft_pos > af_min).sum(axis=0) >= 4).nonzero()[0]
            for ai in cand:
                af = aft_pos[:, ai]
                af[af < 2e-3] = 1e-3
                aft.append([pos, alpha[ai], zip(times[ind], af)])
                times_set |= set(times[ind])
            
        data = {'aft': aft,
                'tmax': times.max(),
                'times': sorted(times_set),
                'len': act.shape[2]}
        return data


class SFSModel(object):
    def get_sfs_filename(self, full=True):
        fn = 'sfs_allfrags_allpats.npz'
        return data_folder[full]+'one_site/'+fn


    def get_data(self):
        import numpy as np
        npz = np.load(self.get_sfs_filename())

        # We need to zip the data to make it easier in D3.js
        keys_zip = set()
        keys_nonzip = set()
        for key in npz:
            for suffix in ['_bin_centers', '_sfs']:
                if suffix in key:
                    keys_zip.add(key[:-(len(suffix))])
                    break
            else:
                keys_nonzip.add(key)

        data = {}
        for key in keys_zip:
            # We need to mask nonstring characters (nice JS)
            keynew = key.replace('.', '')
            data[keynew] = zip(npz[key+'_bin_centers'], npz[key+'_sfs'])
        for key in keys_nonzip:
            keynew = key.replace('.', '')
            data[keynew] = list(npz[key])

        return data


class PropagatorModel(object):
    def get_propagator_filename(self, dt, full=True):
        fn = 'propagator_allfrags_allpats_dt_'+'_'.join(map(str, dt))+'.npz'
        return data_folder[full]+'one_site/'+fn


    def get_data(self, dt=[500, 1000], full=True):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_propagator_filename(dt, full=full))

        # We need to zip the data to make it easier in D3.js
        keys_zip = set()
        keys_nonzip = set()
        for key in npz:
            for suffix in ['_final_frequency',
                           '_prop',
                           '_initial_frequency',
                           '_dt']:
                if suffix in key:
                    keys_zip.add(key[:-(len(suffix))])
                    break
            else:
                keys_nonzip.add(key)

        data = {}
        for key in keys_zip:
            data_key = []
            binsc = npz[key+'_final_frequency']
            dt = list(npz[key+'_dt'])
            for (xi, prop) in izip(npz[key+'_initial_frequency'],
                                   npz[key+'_prop']):

                # Zip deleting extreme bins
                prop_zipped = zip(binsc[1:-1], prop[1:-1])
                data_key.append({'dt': dt,
                                 'xi': xi,
                                 'prop': prop_zipped})
            data[key] = data_key
                
        for key in keys_nonzip:
            data[key] = list(npz[key])

        return data


class DivdivLocalModel(object):
    def __init__(self, pname, observables=['dg', 'ds'],
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
                'len': len_times}

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
        fn = self.pname+'/samples/'+str(i_time+1)+'/'+fragment+'.bam'
        return data_folder[full]+'patients/'+fn


    def get_timeline_filename(self, full=True):
        '''Get the filename of the timeline'''
	fn = self.pname+'/timeline.tsv'
	return data_folder[full]+'patients/'+fn


    def get_coordinate_map_filename(self, full=True, refname='HXB2', format='tsv'):
        '''Get the filename of the coordinate map'''
	fn = 'coordinate_map_'+self.pname+'_'+refname+'_genomewide.'+format
	return data_folder[full]+'coordinate_maps/'+fn


    def get_local_haplotype_filename(self, tmp_root_folder=None, full=True):
        '''Get the filename of a temporary file with the haplotype data'''
        import random
	if tmp_root_folder is None:
	    tmp_root_folder = tmp_folder[full]

        dirname = tmp_root_folder+str(random.randint(0, 10000))+'/'
        while os.path.isdir(dirname):
            dirname = tmp_folder[full]+str(random.randint(0, 10000))+'/'

        self.dirnames.append(dirname)
        os.mkdir(dirname)

        fn = dirname+self.pname+'_'+'_'.join(map(str, self.roi))+'.zip'
        return fn


    def translate_coordinates(self):
        import numpy as np
        # NOTE: the map is always genomewide??
        mapco = dict(np.loadtxt(self.get_coordinate_map_filename(full=True)))

        # In case the coordinate is missing, extend the region
        # Translate start position
        pos = self.roi[1]
        while pos not in mapco:
            pos -= 1
        pmin = min(mapco.iterkeys())
        pos = max(pos, pmin)
        start = mapco[pos]

        # Translate end position
        pos = self.roi[2]
        while pos not in mapco:
            pos += 1
        pmax = max(mapco.iterkeys())
        pos = min(pos, pmax)
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
    def store_alignments(alis, filename):
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
            file_format = file_formats[file_format]
	    alis = [ali['ali'] for ali in alis if ali is not None]
            AlignIO.write(alis, filename, file_format)
    
        elif file_format == 'zip':
            import StringIO
            import zipfile, zlib

            with zipfile.ZipFile(filename, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
                for i, ali in enumerate(alis):
		    if ali is None:
		        continue

                    f = StringIO.StringIO()
                    AlignIO.write(ali['ali'], f, 'fasta')
                    zf.writestr(str(ali['time'])+'_days.fasta', f.getvalue())
    
        else:
            raise ValueError('File format not recognized')



    def get_data(self, clean=True):
        '''Get the data'''
	import numpy as np
        from .analysis.get_local_haplotypes import get_local_haplotypes_aligned

        if clean:
            self.clean_temporary_folders()

	times = np.loadtxt(self.get_timeline_filename())

        # FIXME: do this better
        fragment, start, end = self.roi 

        alis = []
	for i_time in xrange(len(times)):
	    bamfilename = self.get_bam_filename(i_time, fragment)
	    if not os.path.isfile(bamfilename):
	        alis.append(None)
	        continue

            label = self.pname+'_'+str(times[i_time])+'_days_'+\
	      fragment+'_'+str(start+1)+'_'+str(end)+'_'

            ali = get_local_haplotypes_aligned(bamfilename, start, end,
                                               label=label,
                                               VERBOSE=2)

            if ali is None:
                alis.append(None)
                continue

	    alis.append({'time': times[i_time],
	                 'ali': ali})

        fn_out = self.get_local_haplotype_filename(tmp_root_folder=tmp_root_folder)
        self.store_alignments(alis, fn_out)

        return fn_out


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


    def get_table_filename(self, full=True, format='tsv'):
        '''Get the filename of the patient table'''
	fn = 'samples_'+self.pname+'.'+format
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
                    tline[fieldinds[ifi]] = '{:1.0f}'.format(float(field))
                table.append(tline)

        return table

