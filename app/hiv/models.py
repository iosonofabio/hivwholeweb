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


class TreeModel(object):
    def __init__(self, pname, fragment):
        self.pname = pname
        self.fragment = fragment


    def get_filename(self, full=True):
        fn = 'consensi_tree_'+self.pname+'_'+self.fragment+'.newick'
        return data_folder[full]+'trees/'+fn

    
    def get_newick_string(self):
        fn = self.get_filename()
        with open(fn, 'r') as f:
            tree = f.read().rstrip('\n')
            # NOTE: Adding an artificial root creates a long branch.
            root_dist = tree.split(':')[-1][:-1]
            if float(root_dist) > 0.01:
                tree = tree[:tree.rfind(':')]+'0.001;'
        return tree


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


    def get_allele_counts_filename(self, full=True):
        fn = 'allele_counts_'+self.pname+'_genomewide.npz'
        return data_folder[full]+'one_site/'+fn


    def get_data(self, cov_min=100, af_min=1e-1):
        from itertools import izip
        import numpy as np
        npz = np.load(self.get_allele_counts_filename())
        times = npz['times']
        act = npz['act']
        alpha = npz['alpha']

        aft = []
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
            
        data = {'aft': aft,
                'tmax': times.max(),
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



