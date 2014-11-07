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
data_folder_full = os.path.dirname(__file__)+data_folder_short
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

        cov = [[t, list(co)] for (t, co) in izip(times, cov)]

        data = {'cov': cov}
        return data

