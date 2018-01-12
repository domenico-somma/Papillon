# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:51:46 2018

@author: domenico
"""

def read_db(path, drop_comparison=[]):   
    return Papillon(path, drop_comparison)

def _FPKM(name_list):
    pass

def _vs(word1, word2=None):
    pass

def _obtain_list(genelist, path):  # To add eventually remove empty one
    pass

class PapillonBuilder:
    """Extract info from cummeRbund tables"""
    
    def __init__(self, path, drop_comparison=[]):
        try:
            self._find_galaxy()
        except:
            self._find_cummerbund()
            
    def __str__(self):
        pass

    def _find_galaxy(self):
        pass

    def _find_cummerbund(self):
        pass      
    
    def _gene_or_isoform(self, what):
        pass

    @staticmethod
    def _significant(df_detected, comparison, what):
        pass
    
    def _generate_df(self, what):
        pass
    
    def _compare(self):
        pass


class Papillon_db(PapillonBuilder):
    """Make a Papillon_db object and permit to change some values"""
    
    def __init__():
        pass
    
    def selected_exist(self, remove=False):
        pass
    
    # Modify functions
    def dropComparison(self, comparison):
        pass
 
    def change_order(self, new_order):
        pass


class Papillon(Papillon_db):
    """Select and plot genes/isoforms from a Papillon_db"""       
    # Select genes functions
    def _select(self, genelist, what, comparison, sign):
        pass
    
    def _sub_select(self, comparison, sign):
        pass
    
    def get_gene(self, genelist=None, comparison=None, sign=None, export=False):
        pass
    
    def get_isoform(self, genelist=None, comparison=None, sign=None, export=False, show_dup=False):
        pass
        
    def search(self, word, where, how="table", export=False):
        pass
    
    def _export(self, thing, export, name=None, image_extension=".png"):  # add .pdf?
        pass

    # Plot functions
    @staticmethod
    def _fusion_gene_id(df, type_selected, change_index=False):
        pass
    
    def onlyFPKM(self, return_as, **option):
        pass        
        
    @staticmethod
    def _z_score(df):
        pass
        
    def heatmap(self, z_score=True, col_cluster=False, method="complete", cmap="seismic", export=False, **options):
        pass  
    
    def plot(self, title="", legend=True, z_score=False, export=False, df=None, size=10, ci=None, **option):
        pass
    
    def __import_excel(self, filename, type_selected):
        pass  
