# -*- coding: utf-8 -*-
"""

@author: domenico
"""
import papillon as pp


class make:
    def __init__(self,gene1,l=[]):
        self.A=H293.isoform_fpkm[H293.isoform_fpkm["gene_short_name"]==gene1]
        self.B=H293.gene_fpkm[H293.gene_fpkm["gene_short_name"]==gene1]
        self.C=H293.isoform_diff[H293.isoform_diff["gene"]==gene1]
        self.D=H293.gene_diff[H293.gene_diff["gene"]==gene1]
        if l!=[]:
            for name in l:
                self.extract(name)
                
        self.export()
        print("Done")

    
    def extract(self,name):
            a=H293.isoform_fpkm[H293.isoform_fpkm["gene_short_name"]==name]
            b=H293.gene_fpkm[H293.gene_fpkm["gene_short_name"]==name]
            c=H293.isoform_diff[H293.isoform_diff["gene"]==name]
            d=H293.gene_diff[H293.gene_diff["gene"]==name]
            self.extend(a,b,c,d)
            
    def extend(self,a,b,c,d):
            self.A=self.A.append(a)
            self.B=self.B.append(b)
            self.C=self.C.append(c)
            self.D=self.D.append(d)
        
    def export(self):
            print(len(self.A))
            self.A.to_csv("isoform_fpkm.csv")
            self.B.to_csv("gene_fpkm.csv")
            self.C.to_csv("isoform_diff.csv")
            self.D.to_csv("gene_diff.csv")
            print("Exported")

H293=cP.read_db("293")    
l=["IL15","CD44","IL17RC","IL2RA","CCL14","CCL15"]
make("IL6",l)
