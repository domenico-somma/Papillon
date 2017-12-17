# -*- coding: utf-8 -*-
"""
@author: domenico.somma@glasgow.ac.uk
"""

# TO DO test images (plot/heatmap) (by hash?)
# queue ?
# TO DO test exception
# TO DO test export

import os
import pandas as pd
import numpy
from scipy.stats import zscore
import unittest
# import imagehash
# from PIL import Image

os.chdir('..')
import papillon as pp

os.chdir('test')
path="Test_files"
test=pp.read_db(path)

class papillon_Test(unittest.TestCase):
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
   
    def test_read_db(self):
        samples_test=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']
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

    def test_get_gene(self):
        test.get_gene()
        self.assertEqual(len(test.selected),3)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_gene("IL17RC")
        self.assertEqual(test.selected.index[0],"MSTRG.10454")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(len(test.selected.columns),18)

        test.get_gene(["IL6","CCL15"])
        self.assertEqual(len(test.selected),2)
        self.assertEqual(test.selected.index[0],"IL6")
        self.assertEqual(test.selected.index[1],"CCL15-2")
        self.assertEqual(len(test.selected.columns),18)
        
        # Gene-Comparison Test        
        test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"IL6")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_gene(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"CCL15-2")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_gene(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(test.selected),0)
        self.assertEqual(len(test.selected.columns),18)
        
         # Gene-Comparison-sign Test        
        test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign=">")
        self.assertEqual(len(test.selected),0)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_gene(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign="<")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"IL6")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_gene()
        self.assertEqual(len(test.selected),3)
        self.assertEqual(len(test.selected.columns),18)

    def test_get_isoform(self):
        test.get_isoform()
        self.assertEqual(len(test.selected),5)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform("IL6")
        self.assertEqual(test.selected.index[0],"NM_000600.3")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform("CCL15")
        self.assertEqual(test.selected.index[0],"NM_032965.4")
        self.assertEqual(test.selected.index[1],"NM_032965.4-2")
        self.assertEqual(len(test.selected),2)
        self.assertEqual(len(test.selected.columns),18)

        test.get_isoform(["IL6","CD44"])
        self.assertEqual(len(test.selected),3)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform(["IL6","CCL15"])
        self.assertEqual(len(test.selected),3)
        self.assertEqual(len(test.selected.columns),18)
        
        # Isoform-Comparison Test        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"NM_000600.3")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 3_vs_Sample 4")
        self.assertEqual(len(test.selected),2)
        self.assertEqual(test.selected.index[0],"NM_032965.4")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 4")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"NM_032965.4-2")
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 2_vs_Sample 3")
        self.assertEqual(len(test.selected),0)
        self.assertEqual(len(test.selected.columns),18)
        
        # Isoform-Comparison-sign Test        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign=">")
        self.assertEqual(len(test.selected),0)
        self.assertEqual(len(test.selected.columns),18)
        
        test.get_isoform(["IL6","CCL15"], comparison="Sample 1_vs_Sample 2", sign="<")
        self.assertEqual(len(test.selected),1)
        self.assertEqual(test.selected.index[0],"NM_000600.3")
        self.assertEqual(len(test.selected.columns),18)
        
        # Final        
        test.get_isoform()
        self.assertEqual(len(test.selected),5)
        self.assertEqual(len(test.selected.columns),18)
    
    def test_onlyFPKM(self):
        test.get_isoform()
        df=test.onlyFPKM("df")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),5)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4-2")
        
        df=test.onlyFPKM("gene name")
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),5)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4-2")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        df=test.onlyFPKM("array")
        self.assertTrue(type(df)==numpy.ndarray)
        self.assertEqual(len(df),5)
        self.assertEqual(list(df[1]),[0.0, 3.0, 0.0, 0.0])
        self.assertEqual(list(df[-1]),[0.0, 0.0, 0.0, 3.0])
        
        #making extra_df
        test.get_isoform()
        extra_df=test.selected.iloc[:4,2:6].T.copy()
        extra_df=pd.DataFrame(data=extra_df.values, index=test.selected.index[:4], columns=test.selected.columns[2:6])
        extra_df['gene_short_name']=test.selected['gene_short_name'][:4]        

        #testing extra_df
        df=test.onlyFPKM("df",extra_df=extra_df)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),4)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df.index[-1],"NM_032965.4")
        
        df=test.onlyFPKM("gene name",extra_df=extra_df)
        self.assertTrue(type(df)==pd.DataFrame)
        self.assertEqual(len(df),4)
        self.assertEqual(len(df.columns),5)
        self.assertEqual(df.index[0],"NM_000600.3")
        self.assertEqual(df["gene_short_name"][0],"IL6")
        self.assertEqual(df.index[-1],"NM_032965.4")
        self.assertEqual(df["gene_short_name"][-1],"CCL15")
        
        df=test.onlyFPKM("array",extra_df=extra_df)
        self.assertTrue(type(df)==numpy.ndarray)
        self.assertEqual(len(df),4)
        self.assertEqual(list(df[0]),[0.0, 0.0, 4.0, 0.0])
        a=list(df[-1])
        b=[0.016800, 0.0, 0.0, 0.0]
        self.failUnlessAlmostEqual(a[0],b[0], places=0)
        
        # Final        
        test.get_isoform()
        self.assertEqual(len(test.selected),5)
        self.assertEqual(len(test.selected.columns),18)
            
    def test_z_score(self):
        test.get_isoform()
        df=test.onlyFPKM("df")
        df1=zscore(df, axis=1, ddof=1)
        df2=test._z_score(df)
        df2=test.onlyFPKM("array",extra_df=df2)
        self.failUnlessAlmostEqual(df1.all(), df2.all(), places=0)
    
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
        
        test.get_isoform()
        self.assertEqual(len(test.selected),5)
        self.assertEqual(len(test.selected.columns),18)

        test.get_gene()
        self.assertEqual(len(test.selected),3)
        self.assertEqual(len(test.selected.columns),18)
    
    def test_fusion_gene_id(self):
        test.get_isoform()
        m=len(test.selected.columns)
        df=test.selected.copy()
        
        df2=test._fusion_gene_id(df, test.type_selected, change_index=False)
        self.assertEqual(len(df2["gene/ID"]),5)
        self.assertEqual(len(df2.columns),19)
        self.assertEqual(df2["gene/ID"][0],"IL6   NM_000600.3")
        self.assertEqual(df2["gene/ID"][-1],"CCL15   NM_032965.4-2")
        self.assertEqual(df2.index[0],"NM_000600.3")
        self.assertEqual(df2.index[-1],"NM_032965.4-2")
        self.assertEqual(df2["gene_short_name"][0],"IL6")
        self.assertEqual(df2["gene_short_name"][-1],"CCL15")
        
        df2=test._fusion_gene_id(df, test.type_selected, change_index=True)
        self.assertEqual(len(df2),5)
        self.assertEqual(len(df2.columns),17)
        self.assertEqual(df2.index[0],"IL6   NM_000600.3")
        self.assertEqual(df2.index[-1],"CCL15   NM_032965.4-2")
        
        test.get_gene()
        df=test.selected.copy()
        df2=test._fusion_gene_id(df, test.type_selected, change_index=False)       
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),18)
        self.assertEqual(df2.index[1],"MSTRG.10454")
        self.assertEqual(df2.index[-1],"CCL15-2")
        
        df2=test._fusion_gene_id(df, test.type_selected, change_index=True)
        self.assertEqual(len(df2),3)
        self.assertEqual(len(df2.columns),17)
        self.assertEqual(df2.index[0],"IL6")
        self.assertEqual(df2.index[-1],"CCL15")
        
        n=len(test.selected.columns)
        self.assertTrue(m==n)

    def test_comparison(self):
        a,b,c=test._compare()
        self.assertEqual(a,{'IL17RC'})
        self.assertEqual(b,{'CD44'})
        self.assertEqual(c,2)
    
    def test_drop_comparison(self):
        def drop(comp):
            test2=pp.read_db(path,drop_comparison=comp)
            test3=pp.read_db(path)
            test3.dropComparison(comp)
            df1=test2.genes_significant.all()
            df2=test3.genes_significant.all()
            self.assertTrue(df1.all()==df2.all())
        drop("Sample 1_vs_Sample 2")
        drop("Sample 1_vs_Sample 3")
        drop("Sample 1_vs_Sample 4")
        drop("Sample 2_vs_Sample 3")
        drop("Sample 2_vs_Sample 4")
        drop("Sample 3_vs_Sample 4")
    
    def test_change_samples_order(self):
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

        test.change_order(["Sample 1","Sample 2","Sample 3","Sample 4"])
        samples_test=['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4']
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

#    def test_plots(self):
#        test.get_gene()
#        test.plot(export=True)
#        hash1 = imagehash.average_hash(Image.open('Test/Papillon/Plot.png'))       
#        hash2 = imagehash.average_hash(Image.open('Test/Papillon/get_gene_plot_test1.png'))
#        self.assertEqual(hash1,hash2)
        
#        test.plot(export=True,z_score=True)
#        hash1 = imagehash.average_hash(Image.open('Test/Papillon/Plot.png'))       
#        hash2 = imagehash.average_hash(Image.open('Test/Papillon/get_gene_plot_test2.png'))
#        self.assertEqual(hash1,hash2)
    
#    def test_heatmap(self):
#        test.get_gene()
#        test.heatmap(export=True,figsize=(10,10))
#        hash1 = imagehash.phash_simple(Image.open('Test/Papillon/small-heatmap.png'))       
#        print(hash1)        
#        hash2 = imagehash.phash_simple(Image.open('Test/Papillon/test-small-heatmap.png'))
#        print(hash2)
#        self.assertEqual(hash1,hash2)
        
#        test.heatmap(z_score=False,export=True)
#        hash1 = imagehash.average_hash(Image.open('Test/Papillon/small-heatmap.png'))       
#        hash2 = imagehash.average_hash(Image.open('Test/Papillon/test2-small-heatmap.png'))
#        self.failUnlessAlmostEqual(hash1,hash2, places=0)
        
#        test.get_isoform()
#        test.heatmap(export=True)
#        hash1 = imagehash.average_hash(Image.open('Test/Papillon/small-heatmap.png'))       
#        hash2 = imagehash.average_hash(Image.open('Test/Papillon/test3-small-heatmap.png'))
#        self.assertEqual(hash1,hash2)
        
#        test.heatmap(z_score=False,export=True)
#        hash1 = imagehash.average_hash(Image.open('Test/Papillon/small-heatmap.png'))       
#        hash2 = imagehash.average_hash(Image.open('Test/Papillon/test4-small-heatmap.png'))
#        self.assertEqual(hash1,hash2)

if __name__ == '__main__':
    unittest.main()
