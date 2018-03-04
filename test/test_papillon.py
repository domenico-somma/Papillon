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
   
    def test_read_folder(self):
        test=pp.read_folder(path)
        samples_test=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']
        self.assertTrue(test.samples==samples_test)
        comparison_test=['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 
                         'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 
                         'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']
        self.assertTrue(test.comparison==comparison_test)
        self.assertEqual(len(test.genes_detect),5)
        self.assertEqual(len(test.genes_significant),3)
        self.assertEqual(len(test.isoforms_detect),28)
        self.assertEqual(len(test.isoforms_significant),5)
        a=len(test.genes_detect.columns)
        b=len(test.genes_significant.columns)
        c=len(test.isoforms_detect.columns)
        d=len(test.isoforms_significant.columns)
        self.assertTrue(a==b and b==c and c==d and d==18)
        print_test=pp.read_folder(path)
        printable="Samples: ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']\nComparison: ['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']\nGenes Detected: 5\nGenes differential expressed: 3\nIsoform Detected: 28\nIsoform differential expressed: 5\n"
        print(print_test.__str__(),"\n",printable)
        self.assertTrue(print_test.__str__()==printable)
        del print_test
        
    def test_get_gene(self):
        sub=test.get_gene()
        self.assertEqual(len(sub.selected),3)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_gene("IL17RC")
        self.assertEqual(sub.selected.index[0],"MSTRG.10454")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(len(sub.selected.columns),18)

        sub=test.get_gene(["IL6","CCL15"])
        self.assertEqual(len(sub.selected),2)
        self.assertEqual(sub.selected.index[0],"IL6")
        self.assertEqual(sub.selected.index[1],"CCL15-2")
        self.assertEqual(len(sub.selected.columns),18)
        
        # Gene-Comparison Test        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"IL6")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"CCL15-2")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(sub.selected),0)
        self.assertEqual(len(sub.selected.columns),18)
        
         # Gene-Comparison-sign Test        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign=">")
        self.assertEqual(len(sub.selected),0)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign="<")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"IL6")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_gene()
        self.assertEqual(len(sub.selected),3)
        self.assertEqual(len(sub.selected.columns),18)
        
        with self.assertRaises(Exception):
            sub=test.get_gene(comparison="Sample 1_vs_Sample 2", sign="Wrong")
        with self.assertRaises(Exception):
            sub=test.get_gene(sign=">")
        with self.assertRaises(Exception):
            sub=test.get_gene(comparison="Wrong")

    def test_get_isoform(self):
        sub=test.get_isoform()
        self.assertEqual(len(sub.selected),5)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform("IL6")
        self.assertEqual(sub.selected.index[0],"NM_000600.3")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform("CCL15")
        self.assertEqual(sub.selected.index[0],"NM_032965.4")
        self.assertEqual(sub.selected.index[1],"NM_032965.4-2")
        self.assertEqual(len(sub.selected),2)
        self.assertEqual(len(sub.selected.columns),18)

        sub=test.get_isoform(["IL6","CD44"])
        self.assertEqual(len(sub.selected),3)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform(["IL6","CCL15"])
        self.assertEqual(len(sub.selected),3)
        self.assertEqual(len(sub.selected.columns),18)
        
        # Isoform-Comparison Test        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"NM_000600.3")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(sub.selected),2)
        self.assertEqual(sub.selected.index[0],"NM_032965.4")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 4")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"NM_032965.4-2")
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(sub.selected),0)
        self.assertEqual(len(sub.selected.columns),18)
        
        # Isoform-Comparison-sign Test        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign=">")
        self.assertEqual(len(sub.selected),0)
        self.assertEqual(len(sub.selected.columns),18)
        
        sub=test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign="<")
        self.assertEqual(len(sub.selected),1)
        self.assertEqual(sub.selected.index[0],"NM_000600.3")
        self.assertEqual(len(sub.selected.columns),18)
        
        # Final        
        sub=test.get_isoform()
        self.assertEqual(len(sub.selected),5)
        self.assertEqual(len(sub.selected.columns),18)
    
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
        extra_df=sub.selected.iloc[:4,2:6].T.copy()
        extra_df=pd.DataFrame(data=extra_df.values, index=sub.selected.index[:4], columns=sub.selected.columns[2:6])
        extra_df['gene_short_name']=sub.selected['gene_short_name'][:4]        

        # testing extra_df
        df=sub.onlyFPKM("df",extra_df=extra_df)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4")
        
        df=sub.onlyFPKM("gene name",extra_df=extra_df)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        df=sub.onlyFPKM("array",extra_df=extra_df)
        self.assertTrue(type(df)==numpy.ndarray)
        self.assertEqual(len(df),4)
        self.assertEqual(list(df[0]),[0.0, 0.0, 4.0, 0.0])
        a=list(df[-1])
        b=[0.016800, 0.0, 0.0, 0.0]
        self.assertAlmostEqual(a[0],b[0], places=0)
        
        # testing remove_FPKM_name
        df=sub.onlyFPKM("df",extra_df=extra_df, remove_FPKM_name=True)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(list(df.columns),sub.samples)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4")
        
        df=sub.onlyFPKM("gene name",extra_df=extra_df, remove_FPKM_name=True)
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
        self.assertEqual(len(sub.selected),5)
        self.assertEqual(len(sub.selected.columns),18)
            
    def test_z_score(self):
        sub=test.get_isoform()
        df=sub.onlyFPKM("df")
        df1=zscore(df, axis=1, ddof=1)
        df2=sub._z_score(df)
        df2=sub.onlyFPKM("array",extra_df=df2)
        self.assertAlmostEqual(df1.all(), df2.all(), places=0)
    
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
        self.assertEqual(len(search_result.columns),18)
        self.assertEqual(search_result.index[0],"IL6")
        
        search_result=test.search(word="CD44",where="isoforms_detected", how="table")
        self.assertTrue(type(search_result)==pd.DataFrame)
        self.assertEqual(len(search_result),7)
        self.assertEqual(len(search_result.columns),18)
        self.assertEqual(search_result.index[0],"NM_000610.3")
        self.assertEqual(search_result.index[-1],"XM_006718390.1")
        
        search_result=test.search(word="IL6",where="genes_significant", how="table")
        self.assertTrue(type(search_result)==pd.DataFrame)
        self.assertEqual(len(search_result),1)
        self.assertEqual(len(search_result.columns),18)
        self.assertEqual(search_result.index[0],"IL6")
        
        search_result=test.search(word="CCL15",where="genes_significant", how="selected")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(len(test.selected.columns),18)
        self.assertEqual(test.selected.index[0],"CCL15-2")
        
        search_result=test.search(word="CD44",where="isoforms_significant", how="selected")
        self.assertEqual(len(test.selected),2)
        self.assertEqual(len(test.selected.columns),18)
        self.assertEqual(test.selected.index[0],"NM_000610.3")
        self.assertEqual(test.selected.index[-1],"NM_001001389.1")
        
        sub=test.get_isoform()
        self.assertEqual(len(sub.selected),5)
        self.assertEqual(len(sub.selected.columns),18)

        sub=test.get_gene()
        self.assertEqual(len(sub.selected),3)
        self.assertEqual(len(sub.selected.columns),18)
    
    def test_fusion_gene_id(self):
        sub=test.get_isoform()
        m=len(sub.selected.columns)
        df=sub.selected.copy()
        
        df2=sub._fusion_gene_id(df, sub.what, change_index=False)
        self.assertEqual(len(df2["gene/ID"]),5)
        self.assertEqual(len(df2.columns),19)
        self.assertEqual(df2["gene/ID"][0],"IL6   NM_000600.3")
        self.assertEqual(df2["gene/ID"][-1],"CCL15   NM_032965.4-2")
        self.assertEqual(df2.index[0],"NM_000600.3")
        self.assertEqual(df2.index[-1],"NM_032965.4-2")
        self.assertEqual(df2["gene_short_name"][0],"IL6")
        self.assertEqual(df2["gene_short_name"][-1],"CCL15")
        
        df2=sub._fusion_gene_id(df, sub.what, change_index=True)
        self.assertEqual(len(df2),5)
        self.assertEqual(len(df2.columns),17)
        self.assertEqual(df2.index[0],"IL6   NM_000600.3")
        self.assertEqual(df2.index[-1],"CCL15   NM_032965.4-2")
        
        sub=test.get_gene()
        df=sub.selected.copy()
        df2=sub._fusion_gene_id(df, sub.what, change_index=False)       
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),18)
        self.assertEqual(df2.index[1],"MSTRG.10454")
        self.assertEqual(df2.index[-1],"CCL15-2")
        
        df2=sub._fusion_gene_id(df, sub.what, change_index=True)
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),17)
        self.assertEqual(df2.index[0],"IL6")
        self.assertEqual(df2.index[-1],"CCL15")
        
        n=len(sub.selected.columns)
        self.assertTrue(m==n)

    def test_comparison(self):
        a,b,c=test._compare()
        self.assertEqual(a,{'IL17RC'})
        self.assertEqual(b,{'CD44'})
        self.assertEqual(c,2)
    
    def test_drop_comparison(self):
        def drop(comp):
            test2=pp.read_folder(path,drop_comparison=comp)
            test3=pp.read_folder(path)
            test3.dropComparison(comp)
            df1=test2.genes_significant.all()
            df2=test3.genes_significant.all()
            self.assertTrue(df1.all()==df2.all())

        def multidrop(comp):
            test2=pp.read_folder(path)
            test3=pp.read_folder(path)
            test2.dropComparison(comp)
            for c in comp:
                test3.dropComparison(c)
            df1=test2.genes_significant.all()
            df2=test3.genes_significant.all()
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
            test2.dropComparison("Wrong")
        del test2
    
    def test_change_samples_order(self):
        test=pp.read_folder(path)
        test.change_order(["Sample 4","Sample 3","Sample 2","Sample 1"])
        samples_test=["Sample 4","Sample 3","Sample 2","Sample 1"]
        self.assertTrue(test.samples==samples_test)
        comparison_test=['Sample 1_vs_Sample 2', 'Sample 1_vs_Sample 3', 'Sample 1_vs_Sample 4', 'Sample 2_vs_Sample 3', 'Sample 2_vs_Sample 4', 'Sample 3_vs_Sample 4']
        self.assertTrue(test.comparison==comparison_test)
        self.assertEqual(len(test.genes_detect),5)
        self.assertEqual(len(test.genes_significant),3)
        self.assertEqual(len(test.isoforms_detect),28)
        self.assertEqual(len(test.isoforms_significant),5)
        a=len(test.genes_detect.columns)
        b=len(test.genes_significant.columns)
        c=len(test.isoforms_detect.columns)
        d=len(test.isoforms_significant.columns)
        self.assertTrue(a==b and b==c and c==d and d==18)

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
                df_norm = sub._z_score(df_)
                df_norm["gene_short_name"] = sub.selected["gene_short_name"]
                df_ = df_norm.copy()
            elif z_score==False:        
                df_ = sub.onlyFPKM(return_as="gene name",remove_FPKM_name=True)
            
            if type_sel == "gene":
                hue = "gene_short_name"
                df_ = sub._fusion_gene_id(df_, type_sel, change_index=False)
            elif type_sel == "isoform":
                hue = "gene/ID"
                df_ = sub._fusion_gene_id(df_, type_sel, change_index=True)
                df_ = df_.reset_index()
            
#            df_ = test._fusion_gene_id(df_, type_sel, change_index=False)
            
            df = pd.melt(df_, id_vars=hue, var_name="Sample", value_name="FPKM")
            g = sns.factorplot(x="Sample", y="FPKM", hue=hue,
                           data=df, ci=None, legend=True, size=10)
            g.savefig(str(test.path + "test_plot.png"))
        
        def image_check():
            im1=Image.open('Test_files/Papillon/test_plot.png')
            im2=Image.open('Test_files/Papillon/Plot.png')
            hash1 = imagehash.average_hash(im1)
            hash2 = imagehash.average_hash(im2)
            print(hash1,hash2)
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
            df_heatmap = sub._fusion_gene_id(df_heatmap, type_sel, change_index=True)
            im1 = sns.clustermap(df_heatmap, col_cluster=False, method="complete", cmap="seismic", z_score=z_score)
            im1.savefig(str(test.path + "test.png"))
              
        def image_check():
            im1=Image.open('Test_files/Papillon/test.png')
            im2=Image.open('Test_files/Papillon/small-heatmap.png')
            hash1 = imagehash.average_hash(im1)
            hash2 = imagehash.average_hash(im2)
            print(hash1,hash2)
            self.assertEqual(hash1,hash2)
        
        sub=test.get_gene()
        
        heatmap_maker(0,"gene")        
        sub.heatmap(export=True)        
        image_check()
        
        heatmap_maker(None,"gene")        
        sub.heatmap(z_score=False,export=True)        
        image_check()

        sub=test.get_isoform()

        heatmap_maker(0,"isoform")        
        sub.heatmap(export=True)        
        image_check()
        
        heatmap_maker(None,"isoform")        
        sub.heatmap(z_score=False,export=True)        
        image_check()

if __name__ == '__main__':
    unittest.main()
