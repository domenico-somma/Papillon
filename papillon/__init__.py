# -*- coding: utf-8 -*-

"""A python version of CummeRbund
to read and plot Galaxy RNA-seq data"""

import os
import pandas as pd
import seaborn as sns
from distutils.version import LooseVersion
#from IPython import get_ipython
#get_ipython().run_line_magic('matplotlib', 'inline')

if LooseVersion(pd.__version__) < LooseVersion("0.17.1"):
    raise Exception("Pandas >= 0.17.1 required")
if LooseVersion(sns.__version__) < LooseVersion("0.8.1"):
    raise Exception("Seaborn >= 0.8.1 required")


def _FPKM(name):
    """Either append or remove '_FPKM' to a string/list of strings"""
    fpkm = []
    if type(name) == list:
        for n in name:
            fpkm.append(_FPKM(n))
        return fpkm
    elif type(name) == str:
        if name[-5:] == "_FPKM":
            name = name[:-5]
        else:
            name = str(name + "_FPKM")
        return name
    else:
        raise Exception("Only list or string!")


def _vs(word1, word2=None):
    """Either append or remove '_vs_' to 2 words"""
    if "_vs_" in word1 and word2 is None:
        word1, word2 = word1.split("_vs_")
        return word1, word2
    elif word2 is not None:
        fusion = str(word1 + "_vs_" + word2)
        return fusion
    else:
        if word2 is None:
            raise Exception("The string doesn't contain '_vs_'")
        else:
            raise Exception("Only strings")


def _obtain_list(genelist, path):  # To add eventually remove empty one
    """obtain a python list from a file, from a string (or from a list)"""
    gene_list = []
    if type(genelist) == str:
        if "." in genelist:
            file = open(str(path + genelist), "r")
            file = file.readlines()
            for gene in file:
                gene_list.append(gene[:-1])
        else:
            gene_list = [genelist]
    elif type(genelist) == list:
        gene_list = genelist
    elif genelist is None:
        return gene_list
    else:
        raise Exception("Only None, [list] or file_name")
    return gene_list


def read_db(path, drop_comparison=[]):
    """
    Read the cummeRbund Database.
    path - accept a str with the folder path, containing the cummeRbund files
    drop_comparison - drop comparison (str) or list of comparisons and 
                      re-calculate significant genes/isoforms
    """
    return papillon(path, drop_comparison)


class papillon:
    def __init__(self, path, drop_comparison=[]):
        """
        read cummeRbund files and return:
        self.path - files path
        self.samples - samples found
        self.comparison - comparisons found
        self.genes_detect - dataframe of genes detected
        self.genes_significant - dataframe of genes significant
        self.isoforms_detect - dataframe of isoforms detected 
        self.isoforms_significant - dataframe of isoforms significant 
        expressed
        """
        files=os.listdir(path)                                  
        galaxy=[]
        for file in files:
            if ".tabular" in file:
                galaxy.append(file)
        if len(galaxy)==4:
            for file in galaxy:
                if "transcript_FPKM_tracking" in file:
                    self.isoform_fpkm = pd.read_csv(str(path +"/"+ file), delimiter='\t', index_col=0)
                elif "gene_FPKM_tracking" in file:
                    self.gene_fpkm = pd.read_csv(str(path +"/"+ file), delimiter='\t', index_col=0)
                elif "gene_differential_expression" in file:
                    self.gene_diff = pd.read_csv(str(path +"/"+ file), delimiter='\t', index_col=0)
                elif "transcript_differential_expression" in file:
                    self.isoform_diff = pd.read_csv(str(path +"/"+ file), delimiter='\t', index_col=0)
        else:
            try:        
                self.isoform_fpkm = pd.read_csv(str(path + "/isoforms.fpkm_tracking"),
                                                delimiter='\t', index_col=0)
                self.isoform_diff = pd.read_csv(str(path + "/isoform_exp.diff"),
                                                delimiter='\t', index_col=0)
                self.gene_fpkm = pd.read_csv(str(path + "/genes.fpkm_tracking"),
                                             delimiter='\t', index_col=0)
                self.gene_diff = pd.read_csv(str(path + "/gene_exp.diff"),
                                             delimiter='\t', index_col=0)
            except:
                raise("File not found")

        self.path = str(path + "/Papillon/")
        if not os.path.exists(self.path):
            os.makedirs(self.path)  # add check if the folder is moved after creation
        print("Creating dataframe...")

        # find samples name (using isoforms, but it's the same with genes)
        self.samples = []
        print("\tsamples found: ")
        for name in list(self.isoform_fpkm.columns):
            if name[-5:] == "_FPKM":
                name = _FPKM(name)
                self.samples.append(name)
                print(name)

        # file samples in sample_1 and 2 columns
        col_sample1 = []
        col_sample2 = []
        for sample in self.samples:
            if sample in list(self.isoform_diff["sample_1"]):
                col_sample1.append(sample)
            if sample in list(self.isoform_diff["sample_2"]):
                col_sample2.append(sample)

        # generate comparisons name list
        print("\n\tcomparisons found: ")
        if type(drop_comparison)==str:
            drop_comparison=[drop_comparison]
        n=len(drop_comparison)
        self.comparison = []
        for sample1 in col_sample1:
            for sample2 in col_sample2:
                if sample1 != sample2:
                    if len(self.isoform_diff[(self.isoform_diff["sample_1"] == sample1) & (self.isoform_diff["sample_2"] == sample2)]) != 0:
                        comparison = _vs(sample1, sample2)
                        if comparison in drop_comparison:
                            n-=1
                        elif comparison not in drop_comparison:
                            self.comparison.append(comparison)
                            print(comparison)
                    else:
                        pass
                else:
                    pass
        if n != 0:
            raise Exception(drop_comparison," not found")
        self.genes_detect, self.genes_significant = self._generate_df("gene")
        self.isoforms_detect, self.isoforms_significant = self._generate_df(
            "isoform")
        self._compare()
        print("\n...Done")

        if __name__ != "__main__":
            del self.isoform_fpkm
            del self.isoform_diff
            del self.gene_fpkm
            del self.gene_diff

    def _gene_or_isoform(self, what):
        """Users should not use this function directly.
        return _fpkm and _diff tables for either gene or isoform"""
        if what == "gene":
            return self.gene_fpkm, self.gene_diff
        elif what == "isoform":
            return self.isoform_fpkm, self.isoform_diff

    @staticmethod        
    def _significant(df_detected, comparison, what):
        """Users should not use this function directly.
        Calculate significant expressed genes.
        """
        if what not in ["gene","isoform"]:
            raise Exception("what= not known")
        df_significant=df_detected[df_detected.loc[:, comparison].any(axis=1)]
        print("\n\tSignificant expressed ", what + "s: ", len(df_significant))
        return df_significant

    def _generate_df(self, what):
        """Users should not use this function directly.
        Make df for genes/isoforms detect and significant"""
        df_fpkm, df_diff = self._gene_or_isoform(what)
        columns = ["gene_short_name", "gene_id"]

        df = pd.DataFrame.copy(df_fpkm[columns])
        df[_FPKM(self.samples)] = df_fpkm[_FPKM(
            self.samples)]  # TO DO Add CI values here

        for comparison in self.comparison:
            sample1, sample2 = _vs(comparison)

            df2 = df_diff[(df_diff["sample_1"] == sample1) &
                          (df_diff["sample_2"] == sample2)]

            df[comparison] = [True if signif ==
                              "yes" else False for signif in df2["significant"]]
            df[str("q-value_" + comparison)] = df2["q_value"] 

        m = 2
        n = len(self.samples) + 2
        TrueFalseMask = df.iloc[:, m:n] > 0  # with at least 1 value>0
        df_detected = df[TrueFalseMask.any(axis=1)]
        print("\n\tDetected ", what + "s: ", len(df_detected))
        
        df_significant = self._significant(df_detected, self.comparison, what)
        return df_detected, df_significant
    
    def _compare(self):
        """ 
        Users should not use this function directly.
        Compare genes and isoforms significantly expressed
        """
        genes = list(self.genes_significant["gene_short_name"])
        isoforms = list(self.isoforms_significant["gene_short_name"])
        genes_not_found=[]
        isoforms_not_found=[]        
        for name in genes:
            if name not in isoforms:
                genes_not_found.append(name)
        genes_not_found=set(genes_not_found)
        n=len(genes_not_found)
        if n != 0:
            print("\n", n, " significant genes are not found in significant isoforms")
            if n<50 and n>0:
                if n != len(set(genes_not_found)):
                    raise Exception("Duplicate genes")
                print(set(genes_not_found))
        for name in isoforms:
            if name not in genes:
                isoforms_not_found.append(name)
        n=len(isoforms_not_found)        
        if n != 0:
            print(n, " significant isoforms are not found in significant genes")
            isoforms_not_found=set(isoforms_not_found)            
            if n<50 and n>0:
                print(set(isoforms_not_found))
        return genes_not_found, isoforms_not_found, n  # Only for tests so far.
    
    def __str__(self):
        a = "Samples: " + str(self.samples) + "\n"
        b = "Comparison: " + str(self.comparison) + "\n"
        c = "Genes Detected: " + str(len(self.genes_detect)) + "\n"
        d = "Genes differential expressed: " + \
            str(len(self.genes_significant)) + "\n"
        e = "Isoform Detected: " + str(len(self.isoforms_detect)) + "\n"
        f = "Isoform differential expressed: " + \
            str(len(self.isoforms_significant)) + "\n"
        try:
            g = str(len(self.selected)) + " " + \
                self.type_selected + " selected\n"
        except:
            g = "None of the genes is selected"
        visual = a + b + c + d + e + f + g
        return visual

    def dropComparison(self, comparison):
        """Drop Comparison (str) or list of comparisons and re-calculate 
        df_significant"""

        def dropComp(comp):
            if comp in self.comparison:
                    del self.isoforms_detect[comp]
                    del self.isoforms_detect[str("q-value_"+comp)]
                    del self.genes_detect[comp]
                    del self.genes_detect[str("q-value_"+comp)]
                    self.comparison.remove(comp)
                    self.selected_exist(remove=True)
                    print(comp, " removed")
            else:
                raise Exception(comp, " not found, please double check it")

        if type(comparison) == list:
            for comp in comparison:
                dropComp(str(comp))
        elif type(comparison) == str:
            dropComp(comparison)
        else:
            raise Exception(comparison, " not found, please double check it")
            
        self.genes_significant=self._significant(self.genes_detect,self.comparison,"gene")
        self.isoforms_significant=self._significant(self.isoforms_detect,self.comparison,"isoform")
        self._compare()
        print("...Done")
    
    def change_order(self, new_order):
        """Change the samples order"""
        self.selected_exist(remove=True)
        n_sampl=len(self.samples)
        if len(new_order) != n_sampl:
            raise Exception("Number of samples doesn't match")
        for sample in new_order:
            if sample not in self.samples:
                raise Exception(sample,"Sample not known")
        
        cols = self.genes_detect.columns.tolist()
        cols = cols[:2] + _FPKM(new_order) + cols[n_sampl+2:]
        self.samples = new_order
        self.genes_detected = self.genes_detect[cols]
        self.genes_significant = self.genes_significant[cols]
        self.isoforms_detect = self.isoforms_detect[cols]
        self.isoforms_significant = self.isoforms_significant[cols]


    def _select(self, genelist, what, comparison, sign):
        """Users should not use this function directly.
        Part of get_gene/get_isoform function
        """
        if what != "gene" and what != "isoform":
            raise Exception("Only what=gene or what=isoform admitted")
        self.type_selected = what
        gene_list = _obtain_list(genelist, path=self.path)
        if what == "gene":
            df = pd.DataFrame.copy(self.genes_significant)
        elif what == "isoform":
            df = pd.DataFrame.copy(self.isoforms_significant)

        if gene_list != []:
            df["Selected"] = [True if name in gene_list else False for name in df["gene_short_name"]]
            df = df[df["Selected"] == True].iloc[:, :-1]
            for name in gene_list:
                if name not in list(df["gene_short_name"]):
                    print("Gene name not found:\t", name)
        self.selected = df.copy()
        self._sub_select(comparison, sign)

    def _sub_select(self, comparison, sign):
        """Users should not use this function directly.
        Part of get_gene/get_isoform function
        """
        ACCEPTED_SIGN = [">", "<", None]

        if comparison is None:
            if sign is None:
                return
            elif sign is not None:
                raise Exception("Sign passed, but not comparison")
        elif comparison is not None:
            if sign not in ACCEPTED_SIGN:
                raise Exception('Only ">" "<" usable.')
            if comparison not in self.comparison:
                raise Exception("Comparison not found")
            self.selected = self.selected[self.selected[comparison] == True]
            sample1, sample2 = _vs(comparison)
            if sign == ">":
                self.selected = self.selected[self.selected[_FPKM(
                    sample1)] > self.selected[_FPKM(sample2)]]
            elif sign == "<":
                self.selected = self.selected[self.selected[_FPKM(
                    sample1)] < self.selected[_FPKM(sample2)]]
            return

    def get_gene(self, genelist=None, comparison=None, sign=None, export=False):
        """This function select genes. Create self.selected and 
        self.type_selected="gene".
        genelist - accept string (gene name), list of gene names or file
                   with a list of gene name
        comparison - accept only 1 comparison as str (already present in 
                     the data)
        sign - usable in combination with comparison, accept either ">" or 
               "<"
        export - True/False whether want or not export the dataframe of 
                 selected genes
        """
        self._select(genelist, "gene", comparison, sign)
        # self.selected.set_index(["gene_short_name"], inplace=True, verify_integrity=True) #I don't know if could be useful
        print("\nNumber of gene selected: ", len(self.selected))
        self._export(self.selected, name="selected_gene", export=export)
 
    def get_isoform(self, genelist=None, comparison=None, sign=None, export=False, show_dup=False):
        """
        This function select isoforms. Create self.selected and 
        self.type_selected="isoform" 
        genelist - accept string (gene name), list of gene names or file
                   with a list of gene name
        comparison - accept only 1 comparison as str (already present in 
                     the data)
        sign - usable in combination with comparison, accept either ">" or 
               "<"
        export - True/False whether want or not export the dataframe of 
                 selected genes
        show_dup - True/False whether want or not highlight duplicated
                   isoforms for the same gene
        """
        self._select(genelist, "isoform", comparison, sign)

        try:
            del self.selected["duplicate"]
        except:
            pass

        if show_dup == True:
            self.selected["duplicate"] = self.selected.duplicated(
                "gene_short_name", keep=False)
        else:
            pass
        # TO DO if remove_dup == True: # it'd remove the one with lower q-value. and if the q-value is the same???

        print("\nNumber of isoform selected: ", len(self.selected))
        self._export(self.selected, name="selected_isoform", export=export)
        # return proprio selected?
        # Return number genes searched, not found

    @staticmethod
    def _fusion_gene_id(df, type_selected, change_index=False):  
        """Users should not use this function directly.
        Append a "gene/ID" column to the dataframe, and use gene 
        name+id(index) as values, usable or not as index"""
        # print(df)
        if type_selected == "gene":
            if change_index == True:
                df.set_index('gene_short_name', inplace=True)
            return df
        elif type_selected == "isoform":
            df["gene/ID"] = df['gene_short_name'].map(str) + "   " + df.index
            if change_index == True:
                df.set_index("gene/ID", inplace=True)
                del df['gene_short_name']
            return df

    def onlyFPKM(self, return_as, **option):
        """Return a DataFrame with only FPKM columns, 
        return as: 
            "df" - pandas DataFrame
            "array" - numpy array
            "gene name" - pandas DataFrame containing gene names
        It uses self.selected, or an extra_df. 
        """
        self.selected_exist()
        df = self.selected.copy()
        if type(option.get("extra_df")) == pd.DataFrame:
            df = option.get("extra_df")
        if return_as == "df":
            df= df.loc[:, _FPKM(self.samples)]
        elif return_as == "array":
            df= df.loc[:, _FPKM(self.samples)].values
        elif return_as == "gene name":
            columns = ['gene_short_name'] + _FPKM(self.samples)
            df = df.loc[:, columns]
        else:
            raise Exception("Return_as not known. Only 'df','array','gene name'")
        
        if option.get("remove_FPKM_name") == True:
            mydic={}
            n=len(self.samples)
            while n!=0:
                n-=1
                mydic[_FPKM(self.samples[n])]=self.samples[n]
            df.rename(columns=mydic,inplace=True)
        return df

    def selected_exist(self,remove=False):
        """Check if self.selected exists"""
        if remove==True:
            try:
                del self.df_detected
                del self.type_selected
            except:
                pass
        elif remove==False:
            try:
                self.selected
                return True
            except AttributeError:
                raise Exception("No gene selected")
        else:
            raise Exception("Remove= value not known")

    def heatmap(self, z_score=True, col_cluster=False, method="complete", cmap="seismic", export=False, **options):
        """Generate heatmap using selected genes/isoforms
        z_score - True/False whether want or not apply z-score normalization
        col_cluster - True/False whether want or not cluster the samples
        method - clustering algorithm - default is complete-linkage
        cmap - map color
        export - True/False whether want or not export the dataframe of 
                 selected genes
        **options - all the options accepted by seaborn.clustermap
        default metric is euclidean.
        """
        if len(self.samples) > 10:
            raise Exception("High-dimensional data. Ronan et al., 2016")

        self.selected_exist()

        print("Number of genes", len(self.selected))
        if len(self.selected) == 0:
            return
        df_heatmap = self.onlyFPKM(return_as="gene name",remove_FPKM_name=True)
        
        df_heatmap = self._fusion_gene_id(
            df_heatmap, self.type_selected, change_index=True)
        
        if z_score == True: z_score = 0
        elif z_score == False: z_score = None
        small = sns.clustermap(df_heatmap, col_cluster=col_cluster, method=method, cmap=cmap, z_score=z_score, **options)
        self._export(small, name="small-heatmap", export=export)
        if len(df_heatmap) < 1000 and len(df_heatmap) > 25:
            big = sns.clustermap(df_heatmap, col_cluster=col_cluster, method=method, cmap=cmap, z_score=z_score, figsize=((len(self.samples)), int(len(df_heatmap.index) / 4)), **options)
            self._export(big, name="big-heatmap", export=export)
        elif len(df_heatmap) > 1000:
            print("Too many genes for a big heatmap")
        

    @staticmethod
    def _z_score(df):
        """
        Users should not use this function directly.
        Z-score calculation
        """
        # I could use z_score from scipy, but I don't want add scipy dependence too.
        # It would be:
        # from scipy.stats import zscore
        # zscore(a, axis=1, ddof=1)
        print("Calculating Z Score...")
        df_mean = df.mean(axis="columns")
        df_std = df.std(axis="columns")
        df = df.sub(df_mean, axis="index")
        df = df.div(df_std, axis="index")
        return df

    def plot(self, title="", legend=True, z_score=False, export=False, df=None, size=10, ci=None, **option):
        """
        LinePlot a selected dataframe of genes. Max number of genes 200
        title - accept a string as title of the plot
        legend - True/False show the legend
        z_score - True/False calculate the z-score normalization
        export - True/False whether want or not export image
        df - accept a dataframe different from self.selected
        **options - all the options accepted by seaborn.factorplot
        """
        if type(df) == pd.DataFrame:
            pass
        elif df is None:
            df = self.selected.copy()
        else:
            raise Exception("df should be a pandas df")

        if z_score == True:
            df_ = self.onlyFPKM(extra_df=df, return_as="df",remove_FPKM_name=True)
            df_norm = self._z_score(df_)
            df_norm["gene_short_name"] = df["gene_short_name"]
            df_ = df_norm.copy()
        elif z_score == False:
            df_ = self.onlyFPKM(extra_df=df, return_as="gene name",remove_FPKM_name=True)

        print("Number of genes to plot: ", len(df_))
        if len(df_) > 50 and len(df_) < 200:
            print("Too many genes. Legend not shown")
            legend = False
        elif len(df_) >= 200:
            print("Too many genes. Plot not shown")
            return

        if self.type_selected == "gene":
            hue = "gene_short_name"
            df_ = self._fusion_gene_id(
                df_, self.type_selected, change_index=False)
        elif self.type_selected == "isoform":
            hue = "gene/ID"  # Change this hue name
            df_ = self._fusion_gene_id(
                df_, self.type_selected, change_index=True)
            df_ = df_.reset_index()
        df = pd.melt(df_, id_vars=[hue], var_name="Sample", value_name="FPKM")
        g = sns.factorplot(x="Sample", y="FPKM", hue=hue,
                           data=df, ci=ci, size=size, legend=legend, **option)
        g.fig.suptitle(title)
        self._export(g, export=export, name="Plot")
        return g

    def search(self, word, where, how="table", export=False):
        """
        search among genes/isoforms names in detected and significant
        word - accept a str to search among the gene names
        where - accept:
            "genes_detected" 
            "genes_significant" 
            "isoforms_detected"
            "isoforms_significant"

        how - accept:
            "table" return the dataframe with the genes found
            "list" return a list of names, no duplicates
            "selected" put the genes found among the differential expressed 
                       genes in self.selected (to plot),
                       working only with where="significant"
        """

        def df_or_list(df_, how_):
            if how_ == "table":
                return df_
            elif how_ == "list":
                names = list(set(df_["gene_short_name"]))
                print(names)
                return names

        # Checking input
        if where not in ["genes_detected", "genes_significant", "isoforms_detected", "isoforms_significant"]:
            raise Exception("where= not known")
        elif how not in ["table", "list", "selected"]:
            raise Exception("how= not known")
        else:
            pass
        word1, word2 = where.split("_")
        if word1 == "genes":
            df = self.genes_detect[self.genes_detect["gene_short_name"].str.contains(
                word)]
            df_sig = self.genes_significant[self.genes_significant["gene_short_name"].str.contains(
                word)]
        elif word1 == "isoforms":
            df = self.isoforms_detect[self.isoforms_detect["gene_short_name"].str.contains(
                word)]
            df_sig = self.isoforms_significant[self.isoforms_significant["gene_short_name"].str.contains(
                word)]

        if len(df) == 0:
            print(word, " not found")
            return
        else:
            print(len(df), " genes/isoforms detected.")
        if len(df_sig) == 0:
            print("None of these are differentially expressed among the samples")
        else:
            print(len(df_sig), " differential expressed.")

        if how == "selected":
            if word2 != "significant":
                raise Exception('how == "selected", but only significant genes/isoforms can be selected.')
            elif word2 == "significant":
                found = df_sig.copy()
                self.selected = found.copy()
        elif word2 == "detected":
            print("Detected genes preview available")
            found = df_or_list(df, how)
        elif word2 == "significant":
            print("Differential expressed genes preview available")
            found = df_or_list(df_sig, how)

        self._export(found, export=export, name="search_result")
        return found

    def _import_excel(self, filename, type_selected):
        """Only for testing. Users should not use this function directly"""
        if type_selected not in ["gene", "isoform"]:
            raise Exception("type_selected can be only 'gene' or 'isoform'")
        self.type_selected = type_selected
        filename = str(self.path + filename)
        self.selected = pd.read_excel(filename, index_col=0)
        try:
            del self.selected["duplicate"]
        except:
            pass
        if __name__ == "__main__":
            print(
                self.selected.index,
                self.selected.columns,
                self.selected.head(),
                len(self.selected.columns)
            )

    def _export(self, thing, export, name=None, image_extension=".png"):  # add .pdf?
        """
        Manage dataframe or image export parameter.
        Users should not use this function directly"""
        if export == False:
            return
        elif export == True:
            if type(thing) == pd.DataFrame:
                thing.to_excel(str(self.path + name + '.xls'),
                               sheet_name='Sheet1')
                print("\nExported as " + name + ".xls\n")
            elif type(thing) == sns.matrix.ClusterGrid or type(thing) == sns.axisgrid.FacetGrid:
                thing.savefig(str(self.path + name + image_extension))
                print("\nExported as " + name + image_extension)
        else:
            raise Exception("export= can be only 'False' or 'True'")
