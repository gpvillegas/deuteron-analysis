#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 13:34:47 2022

Utility functions for root trees and histograms (up to 2d) based on uproot

not root_numpy is no longer supported.

@author: boeglinw
"""

import numpy as np
import uproot as UR
import LT.box as B

import re



#%% helper functions to prepare the root histogram titles for matplotlib 

# replace # with backslash in titles

def has_tex(s):
    return re.findall('[\#,\_,\^]', s) != []

def prep_title(t):
    ft = t.split()
    for i,ff in enumerate(ft):
        if has_tex(ff):
            ff = ff.lower()              # conver to lower case
            ff = re.sub(r'#',r'\\', ff)  # replace # with \\
            ff = r'$' + ff + '$'         # ad dollar signs
            ft[i] = ff
    return r'' + ' '.join(ft)


#%%
def save_dict(filename, dictionary, selection = None):
    """
    Save a dictionary to an npz file based on a selection 

    Parameters
    ----------
    filename : string
        filename for the output file.
    dictionary: dict
        dictionary to save
    selection : list of bools, optional
        bool arrays for the selected lines. The default is None, meaning all data.

    Returns
    -------
    None.

    """
    d = dictionary
    arg_names = list(d.keys())
    if selection is None:
        np.savez(filename, **d)
    else:
        # apply selection to all data, make temporary dictionary
        dd = dict(zip(list(d.keys()),[d[k][selection] for k in d.keys()]))
        np.savez(filename, **dd)


def save_dict_to_root(filename, dictionary, tree_name = "dict_tree", selection = None):
    """
    Save a dictionary as a root tree based on a selection. It is important
    that all key data have the same length

    Parameters
    ----------
    filename : string
        filename for the output file.
    dictionary: dict
        dictionary to save
    selection : list of bools, optional
        bool arrays for the selected lines. The default is None, meaning all data.

    Returns
    -------
    None.

    """
    d = dictionary
    arg_names = list(d.keys())
    if selection is None:
        r_file = UR.recreate(filename)
        r_file[tree_name] = d
        r_file.close()
    else:
        # apply selection to all data, make temporary dictionary
        dd = dict(zip(list(d.keys()),[d[k][selection] for k in d.keys()]))
        r_file = UR.recreate(filename)
        r_file[tree_name] = dd
        r_file.close()        
        
class read_npz:
    
    def __init__(self, filename):
        self.filename = filename
        self.data = np.load(filename)
        keys = list(self.data.keys())
        self.keys = keys
        
    def __getitem__(self, x):
        return self.data[x]
    

        

#%% Class for handlung root trees
class root_tree: 
    
    def __init__(self, file_name, trees = ['etaprpi0_Tree']):
        """
        read a root file containing trees

        Parameters
        ----------
        file_name : string 
            root file name
        trees : list of strings
            Tree names the sould be loaded The default is ['etaprpi0_Tree'].

        Returns
        -------
        None.

        """
        self.filename = file_name
        self.tree_names = trees
        self.root_file = UR.open(self.filename)
        
    def load_trees(self):
        """
        Load all requested trees and store them in a dictionary called
        tree_names. The keys are the tree names

        Returns
        -------
        None.

        """
        self.trees = {}
        for t_name in self.tree_names:
            try :
                self.trees[t_name] = self.root_file[t_name]
            except Exception as err:
                print(f'cannot get tree {t_name}: {err}')
                continue
        # select the first tree by default
        if self.trees != {}:
            k0 = list(self.trees.keys())[0]
            self.selected_tree = self.trees[k0]
            
    def select_tree (self, name):
        """
        select a tree from theas the current working tree. The default tree is the 
        first one in the list.

        Parameters
        ----------
        name : string
            tree name

        Returns
        -------
        None.

        """
        try:
            self.selected_tree = self.trees[name]
        except Exception as err:
            print(f'cannot select tree {name}: {err}')

    def list_branches(self, tree_name = ''):
        """
        list all branches in tree tree_name
        if no tree_name use the selected one        

        Parameters
        ----------
        tree_name : string, optional
            the name of the tree. The default is ''.

        Returns
        -------
        list of branch names.

        """

        if tree_name == '':
            tree = self.selected_tree
        else:
            try:
                tree = self.trees[tree_name]
            except Exception as err:
                print(f'cannot get tree {tree_name}: {err}')
                return []
        k_list=[]    
        for k in tree.keys():
            #print(k)
            k_list.append(k)
        return k_list

    def get_branches(self, branch_list = []):
        """
        
        load branch data from tree and store them either in a dictionary

        
         
         If the branch list is an empty list all branches are loaded
         
        Parameters
        ----------
        branch_list :  list of string, optional
            list of branch names to be loaded. The default is [] (all branches are loaded.
        Returns
        -------
        None.

        """
        
        tree = self.selected_tree
        if branch_list == []:
            tree_data = tree.arrays(library = "np")
        else:
            tree_data = tree.arrays(branch_list, library = "np")
        self.tree_data = tree_data
        return tree_data    

    def delete_branches(self, branch_list):
        """
        remove the list of tree branch data. This is done if you do not
        need the data anymore to free up memory        

        Parameters
        ----------
        branch_list : list of strings
            list of branches to be removed.

        Returns
        -------
        None.

        """


        for b in branch_list:
            self.tree_data.pop(b)
                


#%% read root file

class root_histos:
    """
    Class to read Root histograms and convert to LT.Box histograms:
        
    Example:
    
    Loading all histograms in a file    
    
    >>> import root_util as RU 
    >>> f = RU.root_histo('my_root_file.root')
    >>> f.load_histos()
    
    Print a list of all histograms loaded
    >>> print(f.histos)
    
    Plot a histogram with the name e.g. adc01
    
    >>> f.histos['adc01'].plot()
    
    """
    
    
    def __init__(self, file_name):
        """
        load a root file and store it in self.root_file
         self.classnames contains the list of objects and type loaded

        Parameters
        ----------
        file_name : str
            root file name.

        Returns
        -------
        None.

        """
        self.file_name = file_name
        try:
            self.root_file = UR.open(self.file_name)
            self.classnames = self.root_file.classnames()
        except Exception as e:
            print(f'Problem loading {file_name} : {e}')
            return None
        
        
    def load_histos(self):
        """
        Load all 1d and 2d histograms and convert them to LT.Box histograms. The 
        converted histograms are stored in self.histos with the same name as the original
        Root histogram

        Returns
        -------
        None.

        """
        self.histos = {}
        for k in self.classnames:
            h_name = k.split(';')[0]
            if self.classnames[k].find('TH1') >= 0:
                # most likely a 1d histo
                h = self.root_file[k]
                h_title = prep_title(h.all_members['fTitle'])
                h1 = B.histo(histogram = self.root_file[k].to_numpy(), title = h_title, xlabel = '')
                h1.bin_error = np.sqrt(self.root_file[k].variances())
                self.histos[h_name] = h1
                
            if self.classnames[k].find('TH2') >= 0:
                # most likely a 1d histo
                h = self.root_file[k]
                h_title = prep_title(h.all_members['fTitle'])
                x_label = prep_title(h.all_members['fXaxis'].all_members['fTitle'])
                y_label = prep_title(h.all_members['fYaxis'].all_members['fTitle'])
                h2 = B.histo2d(histogram = self.root_file[k].to_numpy(), title = h_title, \
                               xlabel = x_label, ylabel = y_label)
                h2.bin_error = np.sqrt(self.root_file[k].variances())
                self.histos[h_name] = h2
                
    def list_all(self, print_all = False):
        for i,k in enumerate(self.histos):
            if print_all :
                print(f"Histogram {i}: '{k}' = {self.histos[k]}")
            else:
                print(f"Histogram {i} : '{k}'")

            
               

        
        