# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/11/14
content:    Models for the hivwholeweb site.

            We are in a particular position, because we do not collect data
            from the users and store in a database, but we need to stream data.

            This makes our 'models' a bit odd.
'''
class TreeModel(object):

    def __init__(self, pname, fragment):
        self.pname = pname
        self.fragment = fragment


    def get_filename(self):
        import os
        folder = os.path.dirname(__file__)
        folder = folder+'/static/data/trees/'
        fn = 'consensi_tree_'+self.pname+'_'+self.fragment+'.newick'
        if self.pname == 'all':
            fn = fn.replace('_all_', '_')
        return folder+fn

    
    def get_newick_string(self):
        fn = self.get_filename()
        with open(fn, 'r') as f:
            tree = f.read().rstrip('\n')
            # NOTE: Adding an artificial root creates a long branch.
            root_dist = tree.split(':')[-1][:-1]
            if float(root_dist) > 0.01:
                tree = tree[:tree.rfind(':')]+'0.001;'
        return tree

