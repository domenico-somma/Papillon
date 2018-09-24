# -*- coding: utf-8 -*-
"""
@author: domenico.somma@glasgow.ac.uk
"""

import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import seaborn as sns
import sys
import pandas as pd
import numpy
from scipy.stats import zscore
import unittest
import imagehash
from PIL import Image

sys.path.append(os.path.abspath(os.path.join('..')))

import papillon as pp

path_to_current_file = os.path.realpath(__file__)
current_directory = os.path.dirname(path_to_current_file)
os.chdir(current_directory)

path="Test_files"
test=pp.read_folder(path)


class papillon_Test(unittest.TestCase):
    def test_different_read(self):
        with self.assertRaises(FileNotFoundError):
            pp.read_folder("Not working")
        pp.read_folder(path)
        pp.read_folder(path+"/galaxy")
        pp.read_files([path+"/gene_exp.diff",path+"/genes.fpkm_tracking",path+"/isoform_exp.diff",path+"/isoforms.fpkm_tracking"])

    def test_functions_FPKM(self):
        self.assertEqual(pp._FPKM("ciao"),"ciao_FPKM")
        self.assertEqual(pp._FPKM("ciao_FPKM"),"ciao")
        self.assertEqual(pp._FPKM(["ciao","hello"]),["ciao_FPKM","hello_FPKM"])
        self.assertEqual(pp._FPKM(["ciao_FPKM","hello_FPKM"]),["ciao","hello"])
 
    def test_functions_vs(self):
        self.assertEqual(pp._vs("ciao_vs_hello"),("ciao","hello"))
        self.assertEqual(pp._vs("ciao","hello"),"ciao_vs_hello")
        
    def test_functions__obtain_list(self):
        self.assertEqual(pp._obtain_list("ciao","fake path"),["ciao"])
        self.assertEqual(type(pp._obtain_list("test49.list",test.path)),list)
        self.assertEqual(len(pp._obtain_list("test49.list",test.path)),49)
        self.assertEqual(pp._obtain_list(["ciao","hello"],"fake path"),["ciao","hello"])
   
    def test_papillon_db(self):
        test=pp.read_folder(path)
        samples_test=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']
        self.assertTrue(test.samples==samples_test)
        comparison_test=['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 
                         'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 
                         'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']
        self.assertTrue(test.comparisons==comparison_test)
        self.assertEqual(len(test.genes_detected),5)
        self.assertEqual(len(test.isoforms_detected),28)
        a=len(test.genes_detected.columns)
        c=len(test.isoforms_detected.columns)
        self.assertTrue(a==c and c==24)
        print_test=pp.read_folder(path)
        printable="Samples: ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']\nComparison: ['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']\nGenes Detected: 5\nGenes differential expressed: 3\nIsoform Detected: 28\nIsoform differential expressed: 5\n"
#        print(print_test.__str__(),"\n",printable)
        self.assertTrue(print_test.__str__()==printable)
        del print_test
        
    def test_significant(self):
        test=pp.read_folder(path)
        self.assertEqual(len(test.Manipulate.significant(test,"gene")),3)
        self.assertEqual(len(test.Manipulate.significant(test,"isoform")),5)
        b=len(test.Manipulate.significant(test,"gene").columns)
        d=len(test.Manipulate.significant(test,"isoform").columns)
        self.assertTrue(b==d and d==24)
        
    def test_get_gene(self):
        sub=test.get_gene()
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
        
        # genelist
        sub=test.get_gene("IL17RC")
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_gene(["IL6","CCL15"])
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(sub.df.index[1],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL17RC","CCL15","CD44"])
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(sub.df.index[-1],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        # Gene-Comparison Test
        sub=test.get_gene(comparison=None)
        self.assertEqual(len(sub.df),5)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(sub.df.index[-1],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(comparison="Sample 2_vs_Sample 4")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
            
        sub=test.get_gene(comparison="Sample 1_vs_Sample 2", comparison_sign=">")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(comparison="Sample 1_vs_Sample 2", comparison_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)
          
        # Fold test
        sub=test.get_gene(fold_sign="<") # ?
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(fold_ind=1)
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[1],"IL6")
        self.assertEqual(sub.df.index[0],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(fold_ind=0.5)
        self.assertEqual(len(sub.df),3)
        self.assertEqual(sub.df.index[0],"CCL15-2")
        self.assertEqual(sub.df.index[-1],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(fold_ind=1,fold_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
               
        sub=test.get_gene(fold_ind="1",fold_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
        
        # Combinations
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign=">")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL17RC","CCL15","CD44"], comparison=None)
        self.assertEqual(len(sub.df),3)
        self.assertEqual(sub.df.index[0],"CD44")
        self.assertEqual(sub.df.index[-1],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CD44"], fold_ind="1")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)      
        
        sub=test.get_gene(["IL17RC","CCL15"], fold_ind="1", fold_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["CD44","IL6"], comparison=None, fold_ind="1")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[1],"IL6")
        self.assertEqual(sub.df.index[0],"CD44")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", fold_ind="1")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15","CD44"], comparison=None, fold_ind=3.29, fold_sign=">")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[1],"IL6")
        self.assertEqual(sub.df.index[0],"CCL15-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<", fold_ind=1, fold_sign=">")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"IL6")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15","IL17RC"], comparison="Sample 2_vs_Sample 4", comparison_sign=">", fold_ind=0.5, fold_sign=">")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["IL6","CCL15","IL17RC"], comparison="Sample 2_vs_Sample 4", comparison_sign=">", fold_ind=1, fold_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_gene(["CD44","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<", fold_ind=1, fold_sign="<")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
       
        sub=test.get_gene()
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
        
        with self.assertRaises(Exception):
            sub=test.get_gene(comparison="Sample 1_vs_Sample 2", comparison_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(comparison_sign=">")
        with self.assertRaises(Exception):
            sub=test.get_gene(comparison="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(fold_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(fold_ind=1, fold_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(fold_ind="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(fold_ind=-1)

    def test_get_isoform(self):
        
        # Final        
        sub=test.get_isoform()
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_isoform()
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)
        
        # genelist
        sub=test.get_isoform("IL6")
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform("CCL15")
        self.assertEqual(sub.df.index[0],"NM_032965.4")
        self.assertEqual(sub.df.index[1],"NM_032965.4-2")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_isoform(["IL6","CD44"])
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15"])
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
        
        # Gene-Comparison Test
        sub=test.get_isoform(comparison=None)
        self.assertEqual(len(sub.df),28)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(sub.df.index[-1],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(comparison="Sample 2_vs_Sample 4")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        # Isoform-Comparison Test        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[0],"NM_032965.4")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 4")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
        
        # Isoform-Comparison-sign Test        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign=">")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(len(sub.df.columns),24)
          
        # Fold test
        sub=test.get_isoform(fold_sign="<") # ? So far is ignored. should it return error?
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(fold_ind=1)
        self.assertEqual(len(sub.df),5)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(sub.df.index[-1],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(fold_ind=1.1)
        self.assertEqual(len(sub.df),3)
        self.assertEqual(sub.df.index[0],"NM_000610.3")
        self.assertEqual(sub.df.index[-1],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(fold_ind=1,fold_sign="<")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(sub.df.index[1],"NM_032965.4")
        self.assertEqual(len(sub.df.columns),24)
               
        sub=test.get_isoform(fold_ind="1",fold_sign="<")
        self.assertEqual(len(sub.df),2)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(sub.df.index[1],"NM_032965.4")
        self.assertEqual(len(sub.df.columns),24)
       
        # Combinations       
        sub=test.get_isoform(["IL17RC","CCL15","CD44"], comparison=None)
        self.assertEqual(len(sub.df),24)
        self.assertEqual(sub.df.index[0],"NM_000610.3")
        self.assertEqual(sub.df.index[-1],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_isoform(["CD44","IL6"], comparison=None, fold_ind="1")
        self.assertEqual(len(sub.df),8)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["CD44","IL6"], comparison=None, fold_ind="1.6")
        self.assertEqual(len(sub.df),7)
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["CD44","IL6"], comparison=None, fold_ind="2.1")
        self.assertEqual(len(sub.df),6)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<", fold_ind=1, fold_sign=">")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_000600.3")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15","IL17RC"], comparison="Sample 2_vs_Sample 4", comparison_sign="<", fold_ind=1.5, fold_sign=">")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)
        
        sub=test.get_isoform(["IL6","CCL15","IL17RC"], comparison="Sample 2_vs_Sample 4", comparison_sign="<", fold_ind=1.6, fold_sign="<")
        self.assertEqual(len(sub.df),1)
        self.assertEqual(sub.df.index[0],"NM_032965.4-2")
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_isoform(["CD44","CCL15"], comparison="Sample 1_vs_Sample 2", comparison_sign="<", fold_ind=1, fold_sign="<")
        self.assertEqual(len(sub.df),0)
        self.assertEqual(len(sub.df.columns),24)
       
        sub=test.get_isoform()
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)
        
        with self.assertRaises(Exception):
            sub=test.get_isoform(comparison="Sample 1_vs_Sample 2", comparison_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_isoform(comparison_sign=">")
        with self.assertRaises(Exception):
            sub=test.get_isoform(comparison="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_isoform(fold_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_isoform(fold_ind=1, fold_sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_isoform(fold_ind="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_isoform(fold_ind=-1)
    
    def test_onlyFPKM(self):
        test=pp.read_folder(path)
        sub=test.get_isoform()
        df=sub.onlyFPKM("df")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),5)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4-2")
        
        df=sub.onlyFPKM("gene name")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),5)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4-2")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        df=sub.onlyFPKM("array")
        self.assertTrue(type(df)==numpy.ndarray)
        self.assertEqual(len(df),5)
        self.assertEqual(list(df[1]),[0.0, 3.0, 0.0, 0.0])
        self.assertEqual(list(df[-1]),[0.0, 0.0, 0.0, 3.0])
        
        # making extra_df
        sub=test.get_isoform()
        extra_df=sub.df.iloc[:4,2:6].T.copy()
        extra_df=pd.DataFrame(data=extra_df.values, index=sub.df.index[:4], columns=sub.df.columns[2:6])
        extra_df['gene_short_name']=sub.df['gene_short_name'][:4]        

        # testing extra_df
        df=sub.plot.onlyFPKM(extra_df,sub.samples,"df")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4")
        
        df=sub.plot.onlyFPKM(extra_df,sub.samples,"gene name")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        df=sub.plot.onlyFPKM(extra_df,sub.samples,"array")
        self.assertTrue(type(df)==numpy.ndarray)
        self.assertEqual(len(df),4)
        self.assertEqual(list(df[0]),[0.0, 0.0, 4.0, 0.0])
        a=list(df[-1])
        b=[0.016800, 0.0, 0.0, 0.0]
        self.assertAlmostEqual(a[0],b[0], places=0)
        
        # testing remove_FPKM_name
        df=sub.plot.onlyFPKM(extra_df,sub.samples,"df", remove_FPKM_name=True)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(list(df.columns),sub.samples)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4")
        
        df=sub.plot.onlyFPKM(extra_df,sub.samples,"gene name", remove_FPKM_name=True)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(list(df.columns[1:]),sub.samples)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        # Final        
        sub=test.get_isoform()
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)
            
    def test_z_score(self):
        sub=test.get_isoform()
        df=sub.onlyFPKM("df")
        df1=zscore(df, axis=1, ddof=1)
        df2=sub.plot._z_score(df)
        df2=sub.plot.onlyFPKM(df2,sub.samples,"array")
        self.assertTrue(numpy.allclose(df1, df2))
    
    def test_search(self):
        search_result=test.search(word="IL6",where="genes_detected", how="list")
        self.assertTrue(type(search_result)==list)
        self.assertEqual(len(search_result),1)
        self.assertEqual(search_result[0],"IL6")
        
        search_result=test.search(word="C",where="genes_detected", how="list")
        self.assertTrue(type(search_result)==list)
        self.assertEqual(len(search_result),3)

        search_result=test.search(word="CD44",where="isoforms_detected", how="list")
        self.assertTrue(type(search_result)==list)
        self.assertEqual(len(search_result),1)
        self.assertEqual(search_result[0],"CD44")

        search_result=test.search(word="IL6",where="genes_detected", how="table")
        self.assertTrue(type(search_result)==pd.DataFrame)
        self.assertEqual(len(search_result),1)
        self.assertEqual(len(search_result.columns),24)
        self.assertEqual(search_result.index[0],"IL6")
        
        search_result=test.search(word="CD44",where="isoforms_detected", how="table")
        self.assertTrue(type(search_result)==pd.DataFrame)
        self.assertEqual(len(search_result),7)
        self.assertEqual(len(search_result.columns),24)
        self.assertEqual(search_result.index[0],"NM_000610.3")
        self.assertEqual(search_result.index[-1],"XM_006718390.1")
        
        search_result=test.search(word="IL6",where="genes_significant", how="table")
        self.assertTrue(type(search_result)==pd.DataFrame)
        self.assertEqual(len(search_result),1)
        self.assertEqual(len(search_result.columns),24)
        self.assertEqual(search_result.index[0],"IL6")
        
#        search_result=test.search(word="CCL15",where="genes_significant", how="selected")
#        self.assertEqual(len(test.df),1)
#        self.assertEqual(len(test.df.columns),24)
#        self.assertEqual(test.df.index[0],"CCL15-2")
        
#        search_result=test.search(word="CD44",where="isoforms_significant", how="selected")
#        self.assertEqual(len(test.df),2)
#        self.assertEqual(len(test.df.columns),24)
#        self.assertEqual(test.df.index[0],"NM_000610.3")
#        self.assertEqual(test.df.index[-1],"NM_001001389.1")
        
        sub=test.get_isoform()
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),24)

        sub=test.get_gene()
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),24)
    
    def test_fusion_gene_id(self):
        sub=test.get_isoform()
        m=len(sub.df.columns)
        df=sub.df.copy()
        
        df2=sub.plot._fusion_gene_id(df, sub.what, change_index=False)
        self.assertEqual(len(df2["gene/ID"]),5)
        self.assertEqual(len(df2.columns),25)
        self.assertEqual(df2["gene/ID"][0],"IL6   NM_000600.3")
        self.assertEqual(df2["gene/ID"][-1],"CCL15   NM_032965.4-2")
        self.assertEqual(df2.index[0],"NM_000600.3")
        self.assertEqual(df2.index[-1],"NM_032965.4-2")
        self.assertEqual(df2["gene_short_name"][0],"IL6")
        self.assertEqual(df2["gene_short_name"][-1],"CCL15")
        
        df2=sub.plot._fusion_gene_id(df, sub.what, change_index=True)
        self.assertEqual(len(df2),5)
        self.assertEqual(len(df2.columns),23)
        self.assertEqual(df2.index[0],"IL6   NM_000600.3")
        self.assertEqual(df2.index[-1],"CCL15   NM_032965.4-2")
        
        sub=test.get_gene()
        df=sub.df.copy()
        df2=sub.plot._fusion_gene_id(df, sub.what, change_index=False)       
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),24)
        self.assertEqual(df2.index[1],"MSTRG.10454")
        self.assertEqual(df2.index[-1],"CCL15-2")
        
        df2=sub.plot._fusion_gene_id(df, sub.what, change_index=True)
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),23)
        self.assertEqual(df2.index[0],"IL6")
        self.assertEqual(df2.index[-1],"CCL15")
        
        n=len(sub.df.columns)
        self.assertTrue(m==n)

    def test_comparison(self):
        a,b,c=test.Manipulate._compare(test)
        self.assertEqual(a,{'IL17RC'})
        self.assertEqual(b,{'CD44'})
        self.assertEqual(c,2)
    
    def test_drop_comparison(self):
        def drop(comp):
            test2=pp.read_folder(path,drop_comparison=comp)
            test3=pp.read_folder(path)
            test3.drop_comparison(comp)
            df1=test2.genes_detected.all()
            df2=test3.genes_detected.all()
            self.assertTrue(df1.all()==df2.all())
            df1=test2.isoforms_detected.all()
            df2=test3.isoforms_detected.all()
            self.assertTrue(df1.all()==df2.all())

        def multidrop(comp):
            test2=pp.read_folder(path)
            test3=pp.read_folder(path)
            test2.drop_comparison(comp)
            for c in comp:
                test3.drop_comparison(c)
            df1=test2.genes_detected.all()
            df2=test3.genes_detected.all()
            self.assertTrue(df1.all()==df2.all())
            df1=test2.isoforms_detected.all()
            df2=test3.isoforms_detected.all()
            self.assertTrue(df1.all()==df2.all())
    
        drop("Sample 1_vs_Sample 2")
        drop("Sample 1_vs_Sample 3")
        drop("Sample 1_vs_Sample 4")
        drop("Sample 2_vs_Sample 3")
        drop("Sample 2_vs_Sample 4")
        drop("Sample 3_vs_Sample 4")
        multidrop(["Sample 1_vs_Sample 2","Sample 1_vs_Sample 3"])
        multidrop(["Sample 1_vs_Sample 2","Sample 3_vs_Sample 4"])
        multidrop(["Sample 1_vs_Sample 4","Sample 2_vs_Sample 4","Sample 2_vs_Sample 3"])

        with self.assertRaises(Exception):
            pp.read_folder(path,drop_comparison="Wrong")

        test2=pp.read_folder(path)
        with self.assertRaises(Exception):
            test2.drop_comparison("Wrong")
        del test2
    
    def test_change_samples_order(self):
        test=pp.read_folder(path)
        test.change_order(["Sample 4","Sample 3","Sample 2","Sample 1"])
        samples_test=["Sample 4","Sample 3","Sample 2","Sample 1"]
        self.assertTrue(test.samples==samples_test)
        comparison_test=['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']
        self.assertTrue(test.comparisons==comparison_test)
        self.assertEqual(len(test.genes_detected),5)
        self.assertEqual(len(test.Manipulate.significant(test,"gene")),3)
        self.assertEqual(len(test.isoforms_detected),28)
        self.assertEqual(len(test.Manipulate.significant(test,"isoform")),5)
        a=len(test.genes_detected.columns)
        b=len(test.Manipulate.significant(test,"gene").columns)
        c=len(test.isoforms_detected.columns)
        d=len(test.Manipulate.significant(test,"isoform").columns)
        self.assertTrue(a==b and b==c and c==d and d==24)

        with self.assertRaises(Exception):
            test.change_order(["Sample 4","Sample 3","Sample 2"])
        with self.assertRaises(Exception):
            test.change_order(["Sample 4","Sample 3","Sample 2","Wrong"])
        
        test=pp.read_folder(path)

    def test_plots(self):
        
        def plot_maker(type_sel,z_score):
            test=pp.read_folder(path)
            if type_sel == "gene":
                sub=test.get_gene()
            elif type_sel == "isoform":
                sub=test.get_isoform()
                
            if z_score == True:
                df_ = sub.onlyFPKM(return_as="df",remove_FPKM_name=True)
                df_norm = sub.plot._z_score(df_)
                df_norm["gene_short_name"] = sub.df["gene_short_name"]
                df_ = df_norm.copy()
            elif z_score==False:        
                df_ = sub.onlyFPKM(return_as="gene name",remove_FPKM_name=True)
            
            if type_sel == "gene":
                hue = "gene_short_name"
                df_ = sub.plot._fusion_gene_id(df_, type_sel, change_index=False)
            elif type_sel == "isoform":
                hue = "gene/ID"
                df_ = sub.plot._fusion_gene_id(df_, type_sel, change_index=True)
                df_ = df_.reset_index()
            
            df = pd.melt(df_, id_vars=hue, var_name="Sample", value_name="FPKM")
            g = sns.factorplot(x="Sample", y="FPKM", hue=hue,
                           data=df, ci=None, legend=True, size=10)
            g.fig.suptitle(" Significant in AT LEAST one condition")
            g.savefig(str(test.path + "test_plot.png"))
        
        def image_check():
            im1=Image.open('Test_files/Papillon/test_plot.png')
            im2=Image.open('Test_files/Papillon/Plot.png')
            hash1 = imagehash.average_hash(im1)
            hash2 = imagehash.average_hash(im2)
#            print(hash1,hash2)
            self.assertEqual(hash1,hash2)
        
        sub=test.get_gene()
        
        plot_maker("gene",False)
        sub.lineplot(export=True)
        image_check()
        
        plot_maker("gene",True)
        sub.lineplot(export=True,z_score=True)
        image_check()
        
        sub=test.get_isoform()
        
        plot_maker("isoform",False)
        sub.lineplot(export=True)
        image_check()
        
        plot_maker("isoform",True)
        sub.lineplot(export=True,z_score=True)
        image_check()
    
    def test_heatmap(self):
        
        def heatmap_maker(z_score, type_sel):
            test=pp.read_folder(path)
            if type_sel == "gene":
                sub=test.get_gene()
            elif type_sel == "isoform":
                sub=test.get_isoform()
            df_heatmap = sub.onlyFPKM(return_as="gene name",remove_FPKM_name=True)
            df_heatmap = sub.plot._fusion_gene_id(df_heatmap, type_sel, change_index=True)
            im1 = sns.clustermap(df_heatmap, col_cluster=False, method="complete", cmap="seismic", z_score=z_score)
            im1.savefig(str(test.path + "test.png"))
              
        def image_check():
            im1=Image.open('Test_files/Papillon/test.png')
            im2=Image.open('Test_files/Papillon/small-heatmap.png')
            hash1 = imagehash.average_hash(im1)
            hash2 = imagehash.average_hash(im2)
#            print(hash1,hash2)
            self.assertEqual(hash1,hash2)
        
        sub=test.get_gene()
        
        heatmap_maker(0,"gene")        
        sub.heatmap(z_score=True,export=True)        
        image_check()
        
        heatmap_maker(None,"gene")        
        sub.heatmap(z_score=False,export=True)        
        image_check()

        sub=test.get_isoform()

        heatmap_maker(0,"isoform")        
        sub.heatmap(z_score=True,export=True)        
        image_check()
        
        heatmap_maker(None,"isoform")        
        sub.heatmap(z_score=False,export=True)        
        image_check()
    
    def test_print(self):
        test=pp.read_folder(path)
        sub=test.get_gene()
        printable2="Type of selection: gene\nNumber of gene selected: 3\nSamples: ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']\nComparison selected: ['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']\n"
        self.assertTrue(sub.__str__()==printable2)
        
    def test_add(self):
        test=pp.read_folder(path)
        # test isoform
        sub1=test.get_isoform("IL6")
        sub2=test.get_isoform("CD44")
        sub=sub1+sub2
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),len(sub1.df.columns))
        sub3=test.get_isoform("CCL15")
        sub=sum([sub1,sub2,sub3])
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),len(sub1.df.columns))
        sub4=test.get_isoform("IL6")
        sub=sub+sub4
        self.assertEqual(len(sub.df),5)
        self.assertEqual(len(sub.df.columns),len(sub1.df.columns))
        
        # test genes
        sub_g1=test.get_gene("IL6")
        sub_g2=test.get_gene("IL17RC")
        sub=sub_g1+sub_g2
        self.assertEqual(len(sub.df),2)
        self.assertEqual(len(sub.df.columns),len(sub_g1.df.columns))
        sub_g3=test.get_gene("CCL15")
        sub=sum([sub_g1,sub_g2,sub_g3])
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),len(sub_g1.df.columns))
        sub_g4=test.get_gene("IL17RC")
        sub=sub+sub_g4
        self.assertEqual(len(sub.df),3)
        self.assertEqual(len(sub.df.columns),len(sub1.df.columns))
        
        with self.assertRaises(Exception):
            sub=sub1+sub_g2
        with self.assertRaises(Exception):    
            sub_g1=test.get_gene("IL6")
            sub_g2=test.get_gene(comparison="Sample 3_vs_Sample 4")
            sub=sub_g1+sub_g2
            
    def test_list_search(self):
        test=pp.read_folder(path)
        sub=test.get_isoform()
        sub_search=sub.search("sfd")
        self.assertEqual(len(sub_search.df),0)
        sub_search=sub.search("00")
        self.assertEqual(len(sub_search.df),3)
        self.assertEqual(len(sub.df.columns),len(sub_search.df.columns))
        
        sub=test.get_gene()
        sub_search=sub.search("sfd")
        self.assertEqual(len(sub_search.df),0)
        sub_search=sub.search("il")
        self.assertEqual(len(sub_search.df),2)
        self.assertEqual(len(sub.df.columns),len(sub_search.df.columns))

    def test_sub_select(self):
        test=pp.read_folder(path)
        sub=test.get_isoform()
        a=sub.select("IL6")
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6"])
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6","wrong"])
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6","CCL15"])
        self.assertEqual(len(a.df),3)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        b=test.get_isoform(["IL6","CCL15"])
        a=sub.select(b)
        self.assertEqual(len(a.df),3)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        
        sub=test.get_gene()
        a=sub.select("IL6")
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6"])
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6","wrong"])
        self.assertEqual(len(a.df),1)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        a=sub.select(["IL6","CCL15"])
        self.assertEqual(len(a.df),2)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))
        b=test.get_isoform(["IL6","CCL15"])
        a=sub.select(b)
        self.assertEqual(len(a.df),2)
        self.assertEqual(len(a.df.columns),len(sub.df.columns))

if __name__ == '__main__':
    unittest.main()
