# -*- coding: utf-8 -*-

"""A python version of CummeRbund to read and plot Galaxy/cuffdiff RNA-seq data"""

import os
import warnings
from distutils.version import LooseVersion
import pandas as pd
import seaborn as sns

warnings.simplefilter('default', DeprecationWarning)

if LooseVersion(pd.__version__) < LooseVersion("0.17.1"):
    raise Exception("Pandas >= 0.17.1 required")
if LooseVersion(sns.__version__) < LooseVersion("0.8.1"):
    raise Exception("Seaborn >= 0.8.1 required")


def read_db(path, drop_comparison=None):
    """Deprecated. Use read_folder()"""
    warnings.warn(
        "read_db() deprecated. Use read_folder() instead.", DeprecationWarning)


def read_folder(path, drop_comparison=None):
    """Read the folder containing the cuffdiff/cummeRbund files, and return
    them to _papillon_builder().

    path - accept a str with the folder path, containing the cuffdiff files
    drop_comparison - drop comparison (str) or list of comparisons and
    re-calculate significant genes/isoforms"""
    if drop_comparison is None:
        drop_comparison = []
    try:
        isoform_fpkm = pd.read_csv(
            str(path + "/isoforms.fpkm_tracking"),
            delimiter='\t', index_col=0)
        isoform_diff = pd.read_csv(
            str(path + "/isoform_exp.diff"),
            delimiter='\t', index_col=0)
        gene_fpkm = pd.read_csv(
            str(path + "/genes.fpkm_tracking"),
            delimiter='\t', index_col=0)
        gene_diff = pd.read_csv(str(path + "/gene_exp.diff"),
                                delimiter='\t', index_col=0)
    except FileNotFoundError:
        files = os.listdir(path)
        galaxy = []
        for file in files:
            if ".tabular" in file:
                galaxy.append(file)
        if len(galaxy) >= 4:
            for file in galaxy:
                if "transcript_FPKM_tracking" in file:
                    isoform_fpkm = pd.read_csv(
                        str(path + "/" + file), delimiter='\t', index_col=0)
                elif "gene_FPKM_tracking" in file:
                    gene_fpkm = pd.read_csv(
                        str(path + "/" + file), delimiter='\t', index_col=0)
                elif "gene_differential_expression" in file:
                    gene_diff = pd.read_csv(
                        str(path + "/" + file), delimiter='\t', index_col=0)
                elif "transcript_differential_expression" in file:
                    isoform_diff = pd.read_csv(
                        str(path + "/" + file), delimiter='\t', index_col=0)
    try:
        return _papillon_builder(isoform_fpkm, isoform_diff, gene_fpkm, gene_diff, path, drop_comparison)
    except UnboundLocalError:
        raise Exception("File not found")


def read_files(files, path=None, drop_comparison=None):
    """Accept cuffdiff/cummeRbund files as iterable, and return
    them to _papillon_builder().

    files - accept an iterable with the cuffdiff files
    path - where export Papillon generated files
    drop_comparison - drop comparison (str) or list of comparisons and
    re-calculate significant genes/isoforms"""
    if drop_comparison is None:
        drop_comparison = []

    for file in files:
        file = pd.read_csv(file, delimiter='\t', index_col=0)
        if "significant" in file.columns:
            # It is a diff
            if file.index.all() == file["gene_id"].all():
                # It is a gene
                gene_diff = file.copy()
            else:
                # It is a isoform
                isoform_diff = file.copy()
        elif "coverage" in file.columns:
            # It is a fpkm
            if file.index.all() == file["gene_id"].all():
                gene_fpkm = file.copy()
            else:
                isoform_fpkm = file.copy()

    return _papillon_builder(isoform_fpkm, isoform_diff, gene_fpkm, gene_diff, path, drop_comparison)

def _make_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)

def _papillon_builder(isoform_fpkm, isoform_diff, gene_fpkm, gene_diff, path, drop_comparison):
    """Accept cuffdiff/cummeRbund files, check whether the files are correct,
    find samples name, comparisons, gene/isoform detected and significant
    and initiate the class papillon"""

    # Test
    if "status" not in isoform_diff.columns or "status" not in gene_diff.columns:
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")
    if isoform_diff.index.all() == isoform_diff["gene_id"].all():
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")
    if gene_diff.index.all() != gene_diff["gene_id"].all():
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")
    if "length" not in isoform_fpkm.columns or "length" not in gene_fpkm.columns:
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")
    if isoform_fpkm.index.all() == isoform_fpkm["gene_id"].all():
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")
    if gene_fpkm.index.all() != gene_fpkm["gene_id"].all():
        raise Exception(
            "Something wrong during cuffdiff/cummeRbund files reading")

    # Making folder
    if path is None:
        path = "Papillion/"
    else:
        path = str(path + "/Papillon/")
    _make_folder(path)

    print("Creating dataframe...")
    # find samples name (using isoforms, but it's the same with genes)
    samples = []
    print("\tsamples found: ")
    samples = [name[:-5]
               for name in isoform_fpkm.columns.tolist() if name[-5:] == "_FPKM"]
    [print(name) for name in samples]

    # file samples in sample_1 and 2 columns
    col_sample1 = []
    col_sample2 = []
    for sample in samples:
        if sample in list(isoform_diff["sample_1"]): # TO DO - To change with tolist?
            col_sample1.append(sample)
        if sample in list(isoform_diff["sample_2"]):
            col_sample2.append(sample)

    # generate comparisons name list
    print("\n\tcomparisons found: ")
    if isinstance(drop_comparison, str):
        drop_comparison = [drop_comparison]
    n_left = len(drop_comparison)

    comparisons = []
    for sample1 in col_sample1:
        for sample2 in col_sample2:
            if sample1 != sample2:
                if len(isoform_diff[(isoform_diff["sample_1"] == sample1) & (isoform_diff["sample_2"] == sample2)]) != 0:  # TO DO - Not so clear, can be improved?
                    comparison = _vs(sample1, sample2)
                    if comparison in drop_comparison:
                        n_left -= 1
                    elif comparison not in drop_comparison:
                        comparisons.append(comparison)
                        print(comparison)
                else:  # TO DO - Really necessary?
                    pass
            else:
                pass
    if n_left != 0:
        raise Exception(drop_comparison, " not found")
    genes_detect, genes_significant = _generate_df(
        "gene", samples, gene_fpkm, gene_diff, isoform_fpkm, isoform_diff, comparisons)
    isoforms_detect, isoforms_significant = _generate_df(
        "isoform", samples, gene_fpkm, gene_diff, isoform_fpkm, isoform_diff, comparisons)

    return Papillon_db(path, samples, comparisons, genes_detect, genes_significant, isoforms_detect, isoforms_significant)

#def    _read_folder_testing():
    #delegare e return papillon_db and cummerbund


def _generate_df(what, samples, gene_fpkm, gene_diff, isoform_fpkm, isoform_diff, comparisons):
    """Make dataframe for genes/isoforms detected and significant"""

    def gene_or_isoform(what, gene_fpkm, gene_diff, isoform_fpkm, isoform_diff):
        """return _fpkm and _diff tables for either gene or isoform"""
        if what == "gene":
            return gene_fpkm, gene_diff
        elif what == "isoform":
            return isoform_fpkm, isoform_diff

    df_fpkm, df_diff = gene_or_isoform(
        what, gene_fpkm, gene_diff, isoform_fpkm, isoform_diff)
    columns = ["gene_short_name", "gene_id"]

    df = pd.DataFrame.copy(df_fpkm[columns])
    # TO DO Add CI values here
    df[_FPKM(samples)] = df_fpkm[_FPKM(samples)]

    for comparison in comparisons:
        sample1, sample2 = _vs(comparison)

        df2 = df_diff[(df_diff["sample_1"] == sample1) &
                      (df_diff["sample_2"] == sample2)]

        df[comparison] = [True if signif ==
                          "yes" else False for signif in df2["significant"]]
        df[str("q-value_" + comparison)] = df2["q_value"]

    m = 2
    n = len(samples) + 2
    TrueFalseMask = df.iloc[:, m:n] > 0  # with at least 1 value>0
    df_detected = df[TrueFalseMask.any(axis=1)]
    print("\n\tDetected ", what + "s: ", len(df_detected))

    df_significant = _Manipulate_db._significant(df_detected, comparisons, what)
    return df_detected, df_significant


def _FPKM(name_list):
    """Either append or remove '_FPKM' to a string or an iterable"""
    try:
        if name_list.endswith("_FPKM"):
            return name_list[:-5]
        return name_list + "_FPKM"
    except AttributeError:
        if name_list[0].endswith("_FPKM"):
            return [name[:-5] for name in name_list]
        return [name + "_FPKM" for name in name_list]


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


def _obtain_list(genelist, path):  # TO DO - eventually remove empty one
    """obtain a python list from a file, from a string (or from a list)"""
    gene_list = []
    try:
        if "." in genelist:
            file = open(str(path + genelist), "r")
            gene_list = [gene[:-1] for gene in file.readlines()]
        elif isinstance(genelist, str):
            gene_list = [genelist]
        else:
            gene_list = [e for e in genelist]
    except TypeError:
        pass
    return gene_list


#class _Cummerbund:
#    def __init__(self, isoform_fpkm, isoform_diff, gene_fpkm, gene_diff):
#        self.isoform_fpkm
#        self.isoform_diff
#        self.gene_fpkm
#        self.gene_diff


class Papillon_db:
    """Make a Papillon_db object and permit to change some values

    self.path - files path
    self.samples - samples found
    self.comparisons - comparisons found
    self.genes_detect - dataframe of genes detected
    self.genes_significant - dataframe of genes significant
    self.isoforms_detect - dataframe of isoforms detected
    self.isoforms_significant - dataframe of isoforms significant
    expressed
    redefine __str__"""


    def __init__(self, path, samples, comparisons, genes_detected, genes_significant, isoforms_detected, isoform_significant):
        self.path = path
        self.samples = samples
        self.comparisons = comparisons
        self.genes_detect = genes_detected
        self.genes_significant = genes_significant
        self.isoforms_detect = isoforms_detected
        self.isoforms_significant = isoform_significant

        self.Manipulate = _Manipulate_db()
        self.Manipulate._compare(self)
        print("\n...Done")      

    def __str__(self):
        a = "Samples: " + str(self.samples) + "\n"
        b = "Comparison: " + str(self.comparisons) + "\n"
        c = "Genes Detected: " + str(len(self.genes_detect)) + "\n"
        d = "Genes differential expressed: " + \
            str(len(self.genes_significant)) + "\n"
        e = "Isoform Detected: " + str(len(self.isoforms_detect)) + "\n"
        f = "Isoform differential expressed: " + \
            str(len(self.isoforms_significant)) + "\n"
        visual = a + b + c + d + e + f
        return visual
    
    def drop_comparison(self, comparison):
        """Drop Comparison (str) or list of comparisons and re-calculate
        df_significant

        comparison: comparison (str) or list of comparisons"""
        self = self.Manipulate.dropComparison(self, comparison)
        
    def change_order(self, new_order):
        """Change the samples order

        new_order: list of samples order"""
        self = self.Manipulate.change_order(self, new_order)
        
    def get_gene(self, genelist=None, comparison=None, comparison_sign=None, fold_ind=None, fol_sign=None):
        """This function select genes. It return a Papillon_list object

        genelist - accept string (1 gene name), list of gene names or file
                   with a list of gene names
        comparison - accept only 1 comparison as str (already present in
                     the data)
        sign - usable in combination with comparison, accept either ">" or
               "<"
        fold_ind - fold induction (log2) higher then number
        """
        return self.Manipulate.get_gene(self, genelist, comparison, comparison_sign, fold_ind, fol_sign)

    def get_isoform(self, genelist=None, comparison=None, comparison_sign=None, fold_ind=None, fol_sign=None):
        """This function select isoforms. It creates a Papillon object

        genelist - accept string (gene name), list of gene names or file
                   with a list of gene names
        comparison - accept only 1 comparison as str (already present in
                     the data)
        sign - usable in combination with comparison, accept either ">" or
               "<"
        fold_ind - fold induction (log2) higher then number"""
        
        return self.Manipulate.get_isoform(self, genelist, comparison, comparison_sign, fold_ind, fol_sign)

    def search(self, word, where, how="table", export=False):
        """search among genes/isoforms names in detected and significant

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
                       working only with where="significant" """
        return self.Manipulate.search(self, word, where, how, export)

class _Manipulate_db:
    # TO DO - Add the function to export the Papillon_db (as table? as sqlite?)

    @staticmethod
    def _compare(pp):
        """Compare genes and isoforms significantly expressed"""
        genes = list(pp.genes_significant["gene_short_name"])
        isoforms = list(pp.isoforms_significant["gene_short_name"])
        genes_not_found = []
        isoforms_not_found = []
        for name in genes:
            if name not in isoforms:
                genes_not_found.append(name)
        genes_not_found = set(genes_not_found)
        n = len(genes_not_found)
        if n != 0:
            print(
                "\n", n, " significant genes are not found in significant isoforms")
            if n < 50 and n > 0:
                if n != len(set(genes_not_found)):
                    raise Exception("Duplicate genes")
                print(set(genes_not_found))
        for name in isoforms:
            if name not in genes:
                isoforms_not_found.append(name)
        n = len(isoforms_not_found)
        if n != 0:
            print(
                n, " significant isoforms are not found in significant genes")
            isoforms_not_found = set(isoforms_not_found)
            if n < 50 and n > 0:
                print(set(isoforms_not_found))
        return genes_not_found, isoforms_not_found, n  # Only for tests so far.

    def dropComparison(self, pp, comparison):
        """Drop Comparison (str) or list of comparisons and re-calculate
        df_significant

        comparison: comparison (str) or list of comparisons
        """

        def dropComp(comp):
            if comp in pp.comparisons:
                del pp.isoforms_detect[comp]
                del pp.isoforms_detect[str("q-value_" + comp)]
                del pp.genes_detect[comp]
                del pp.genes_detect[str("q-value_" + comp)]
                pp.comparisons.remove(comp)
                print(comp, " removed")
            else:
                raise Exception(comp, " not found, please double check it")

        if isinstance(comparison, str):
            dropComp(comparison)
        else:
            for comp in comparison:
                dropComp(comp)

        pp.genes_significant = self._significant(
            pp.genes_detect, pp.comparisons, "gene")
        pp.isoforms_significant = self._significant(
            pp.isoforms_detect, pp.comparisons, "isoform")
        self._compare(pp)
        print("...Done")
        return pp

    @staticmethod
    def _significant(df_detected, comparison, what):
        """Calculate significant expressed genes."""
        if what not in ["gene", "isoform"]:
            raise Exception("what= not known")
        df_significant = df_detected[
            df_detected.loc[:, comparison].any(axis=1)]
        print("\n\tSignificant expressed ", what + "s: ", len(df_significant))
        return df_significant

    @staticmethod
    def change_order(pp, new_order):
        """Change the samples order

        new_order: list of samples order"""
        n_sampl = len(pp.samples)
        if len(new_order) != n_sampl:
            raise Exception("Number of samples doesn't match")
        for sample in new_order:
            if sample not in pp.samples:
                raise Exception(sample, "Sample not known")

        cols = pp.genes_detect.columns.tolist()
        cols = cols[:2] + _FPKM(new_order) + cols[n_sampl + 2:]
        pp.samples = new_order
        pp.genes_detected = pp.genes_detect[cols]
        pp.genes_significant = pp.genes_significant[cols]
        pp.isoforms_detect = pp.isoforms_detect[cols]
        pp.isoforms_significant = pp.isoforms_significant[cols]
        return pp
     
    @staticmethod
    def _select(pp, genelist, what):#, comparison, sign):
        """Part of get_gene/get_isoform function"""
        
        if what != "gene" and what != "isoform":
            raise Exception("Only what=gene or what=isoform admitted")
        gene_list = _obtain_list(genelist, path=pp.path)
        if what == "gene":
            df = pd.DataFrame.copy(pp.genes_significant)
        elif what == "isoform":
            df = pd.DataFrame.copy(pp.isoforms_significant)      
        n=0
        if gene_list != []:
            df["Selected"] = [True if name in gene_list else False for name in df["gene_short_name"]]
            df = df[df["Selected"] == True].iloc[:, :-1]
            for name in gene_list:
                if name not in list(df["gene_short_name"]):
                    print("Gene name not found:\t", name) # return not found list 
                    n+=1
            print("Number of gene not found: ",n)
        return df
        
    def search(self, pp ,word, where, how="table", export=False):
        """search among genes/isoforms names in detected and significant

        word - accept a str to search among the gene names
        where - accept:
            "genes_detected"
            "genes_significant"
            "isoforms_detected"
            "isoforms_significant"

        how - accept:
            "table" return the dataframe with the genes found
            "list" return a list of names, no duplicates"""
#            "selected" put the genes found among the differential expressed
#                       genes in self.selected (to plot),
#                       working only with where="significant" """

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
            df = pp.genes_detect[pp.genes_detect["gene_short_name"].str.contains(
                word)]
            df_sig = pp.genes_significant[pp.genes_significant["gene_short_name"].str.contains(
                word)]
        elif word1 == "isoforms":
            df = pp.isoforms_detect[pp.isoforms_detect["gene_short_name"].str.contains(
                word)]
            df_sig = pp.isoforms_significant[pp.isoforms_significant["gene_short_name"].str.contains(
                word)]

        if len(df) == 0:
            print(word, " not found")
            return
        else:
            print(len(df), " genes/isoforms detected.")
        if len(df_sig) == 0:
            print(
                "None of these are differentially expressed among the samples")
        else:
            print(len(df_sig), " differential expressed.")

        if how == "selected":
            if word2 != "significant":
                raise Exception(
                    'how == "selected", but only significant genes/isoforms can be selected.')
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
        print(found)
        return found
    
    def get_gene(self, pp, genelist=None, comparison=None, comparison_sign=None, fold_ind=None, fold_sign=None):
        """This function select genes. It creates a Papillon object

        genelist - accept string (1 gene name), list of gene names or file
                   with a list of gene names
        comparison - accept only 1 comparison as str (already present in
                     the data)
        sign - usable in combination with comparison, accept either ">" or
               "<"
        """
        df = self._select(pp, genelist, "gene")#, comparison, comparison_sign) # Why I have done this pre-selection?
        return Papillon_list(df, "gene", pp.comparisons, pp.path, pp.samples ,comparison, comparison_sign, fold_sign, fold_ind)
        # To do - Return number genes not found

    def get_isoform(self, pp, genelist=None, comparison=None, comparison_sign=None, fold_ind=None, fold_sign=None):
        """This function select isoforms. It creates a Papillon object

        genelist - accept string (gene name), list of gene names or file
                   with a list of gene names
        comparison - accept only 1 comparison as str (already present in
                     the data)
        sign - usable in combination with comparison, accept either ">" or
               "<"
        export - True/False whether want or not export the dataframe of
                 selected genes"""
        df = self._select(pp, genelist, "isoform")#, comparison, comparison_sign) # Why I have done this pre-selection?
        return Papillon_list(df, "isoform", pp.comparisons, pp.path, pp.samples, comparison, comparison_sign, fold_sign, fold_ind)
    
    def _export(self, thing, export, name=None):
        """Manage dataframe or image export parameter."""
        if export is False: 
            return
        elif export is True:
            _make_folder(self.path)
            thing.to_excel(str(self.path + name + '.xls'), sheet_name='Sheet1')
            print("\nExported as " + name + ".xls\n")
        else:
            raise Exception("export= can be only 'False' or 'True'")

    
class Papillon_list:
    def __init__(self, df, what, comparisons, path, samples, comparison=None, comparison_sign=None, fold_ind=None, fold_sign=">"):
        self.df=df
        self.what=what
        self.path=path
        self.samples=samples
#        self.genelist=genelist
        if comparison is None: 
            if comparison_sign is not None:
                raise Exception("Sign passed but not comparison")
            else:
                self.comparison=comparisons
        elif comparison is not None:
            if comparison not in comparisons:
                raise Exception("Comparison not found")
            else:
                self.comparison=comparison
        self.comparison_sign=comparison_sign
        self.fold_ind=fold_ind
        self.fold_sign=fold_sign
        
        self.plot=_Plot()
        self.Manipulate = _Manipulate_list()
        self.Manipulate._sub_select(self, comparison, comparison_sign, fold_ind, fold_sign)
        
    def __str__(self):
        a = "Number of "+ self.what + " selected: "+ str(len(self.df)) + "\n"
        b = "Samples: " + str(self.samples) + "\n"
        visual = a + b
        if self.comparison_sign is not None:
            w1,w2=_vs(self.comparison)
            n = "Comparison selected: " + w1 + self.comparison_sign + w2 + "\n"
        else:
            n = "Comparison selected: " + str(self.comparison) + "\n"
        visual = visual + n
        if self.fold_ind is not None:
            n = "Fold induction log2" + self.comparison_sign + str(self.fold_ind) + "\n"
            visual = visual + n
        self.show()
        return visual
        
    def __add__(self, other):
        if self.what != other.what:
            raise Exception("Impossible, one is gene, the other isoform")
        elif self.samples != other.samples or self.comparison != other.comparison or self.path != other.path:
            raise Exception("The two elements seems to have different origins")
        df= pd.merge(self.df, other.df, how='outer')
        return Papillon_list(df, what=self.what, comparisons=self.comparison, comparison=None, path=self.path, samples=self.samples)
            
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other) # sum([T1, T2, T3])
        
    def show(self):
        self.Manipulate.show(self)

#    def __getattr__(self, arg): # TO DO - Add signle functions with descriptions
#        _plot=_Plot()
#        setattr(_plot, "pp", self)
#        return getattr(_plot, arg)
    
    def onlyFPKM(self, return_as, remove_FPKM_name=False):
        """Take a Papillon dataframe and a list of samples, return only FPKM columns.

        return as:
            "df" - pandas DataFrame
            "array" - numpy array
            "gene name" - pandas DataFrame containing gene names"""
        return self.plot.onlyFPKM(self.df, self.samples, return_as, remove_FPKM_name)
    
    def plot(**parameter):
        """Deprecated. Use self.lineplot()"""
        warnings.warn('self.plot() is deprecated. Use self.lineplot() instead.', DeprecationWarning)

    def heatmap(self, z_score=False, col_cluster=False, method="complete", cmap="seismic", export=False, **options):
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
        self.plot.heatmap(self, z_score, col_cluster, method, cmap, export, **options)

    def lineplot(self, title="", legend=True, z_score=False, export=False, df=None, size=10, ci=None, **option):
        """
        LinePlot selected genes expression levels. Max number of genes 200

        title - accept a str as title of the plot
        legend - True/False show the legend
        z_score - True/False calculate the z-score normalization
        export - True/False whether or not export the image
        df - accept an exernal dataframe, different from self.selected
        **options - all the options accepted by seaborn.factorplot"""
        self.plot.lineplot(self, title, legend, z_score, export, df, size, ci, **option)

    
class _Manipulate_list:
    def _sub_select(self, pp, comparison, comparison_sign, fold_ind, fold_sign):
        """ Part of get_gene/get_isoform function"""
        
        ACCEPTED_SIGN = [">", "<", None]

        if comparison_sign is not None:
            if comparison_sign not in ACCEPTED_SIGN:
                raise Exception('Only ">" "<" usable.')
            if comparison is None:
                raise Exception("Comparison_sign passed, but not comparison")
        
        selected = pp.df
        if comparison is not None:
            selected = selected[selected[comparison] == True]
            sample1, sample2 = _vs(comparison)
            if comparison_sign == ">":
                selected = selected[selected[_FPKM(sample1)] > selected[_FPKM(sample2)]]
            elif comparison_sign == "<":
                selected = selected[selected[_FPKM(sample1)] < selected[_FPKM(sample2)]]
        if fold_ind is not None:
                fi=[str("fi_log2_")+comp for comp in self.comparison]
                if fold_sign is ">":
                    TrueFalseMask1=selected[fi]>fold_ind# or [df[fi]<-fold_ind]
                    TrueFalseMask2=selected[fi]<-fold_ind
                if fold_sign is "<":
                    TrueFalseMask1=selected[fi]<fold_ind# or [df[fi]<-fold_ind]
                    TrueFalseMask2=selected[fi]>-fold_ind

                TrueFalseMask1=TrueFalseMask1.any(axis=1)
                TrueFalseMask2=TrueFalseMask2.any(axis=1)
#                TrueFalseMask=pd.merge(TrueFalseMask1,TrueFalseMask2,how="left")
                TrueFalseMask=selected.copy()
                TrueFalseMask["A"]=TrueFalseMask1
                TrueFalseMask["B"]=TrueFalseMask2
                TrueFalseMask=TrueFalseMask.loc[:,["A","B"]]
                TrueFalseMask=TrueFalseMask.any(axis=1)
                selected=selected[TrueFalseMask]
        pp.df=selected.copy()
                    
        print("\nNumber of ", pp.what," selected: ", len(pp.df))
    
    def show(self, pp):
        print(pp.df)
    
    def _export(self, thing, export, name=None): # TO DO - Not working. To activate
        """Manage dataframe or image export parameter."""
        if export is False: 
            return
        elif export is True:
            _make_folder(self.path)
            thing.to_excel(str(self.path + name + '.xls'), sheet_name='Sheet1')
            print("\nExported as " + name + ".xls\n")
        else:
            raise Exception("export= can be only 'False' or 'True'")
        
    def swap_gene_isoform(self):
        pass
        #TO DO

class _Plot:
    """  """

    @staticmethod
    def _fusion_gene_id(df, what, change_index=False):
        """Append a "gene/ID" column to the dataframe, and use gene
        name+id(index) as values, usable or not as index"""
        # print(df)
        if what == "gene":
            if change_index is True:
                df.set_index('gene_short_name', inplace=True)
            return df
        elif what == "isoform":
            df["gene/ID"] = df['gene_short_name'].map(str) + "   " + df.index
            if change_index is True:
                df.set_index("gene/ID", inplace=True)
                del df['gene_short_name']
            return df

    def onlyFPKM(self, df, samples, return_as, remove_FPKM_name=False):
        """Take a Papillon dataframe and a list of samples, return only FPKM columns.

        return as:
            "df" - pandas DataFrame
            "array" - numpy array
            "gene name" - pandas DataFrame containing gene names"""
        if return_as == "df":
            df = df.loc[:, _FPKM(samples)]
        elif return_as == "array":
            df = df.loc[:, _FPKM(samples)].values
        elif return_as == "gene name":
            columns = ['gene_short_name'] + _FPKM(samples)
            df = df.loc[:, columns]
        else:
            raise Exception(
                "Return_as not known. Only 'df','array','gene name'")

        if remove_FPKM_name is True:
            mydic = {}
            n = len(samples)
            while n != 0:
                n -= 1
                mydic[_FPKM(samples[n])] = samples[n]
                df.rename(columns=mydic, inplace=True)
        return df

    def heatmap(self, pp, z_score=True, col_cluster=False, method="complete", cmap="seismic", export=False, **options):
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
        if len(pp.samples) > 10:
            raise Exception("High-dimensional data. Ronan et al., 2016")
        print("Number of genes", len(pp.df))
        if len(pp.df) == 0:
            return
        df_heatmap = self.onlyFPKM(pp.df, pp.samples, return_as="gene name", remove_FPKM_name=True)

        df_heatmap = self._fusion_gene_id(
            df_heatmap, pp.what, change_index=True)

        if z_score is True: 
            z_score = 0
        elif z_score is False: 
            z_score = None
        small = sns.clustermap(
            df_heatmap, col_cluster=col_cluster, method=method, cmap=cmap, z_score=z_score, **options)
        self._export(small, pp.path, name="small-heatmap", export=export)
        if len(df_heatmap) < 1000 and len(df_heatmap) > 25:
            big = sns.clustermap(
                df_heatmap, col_cluster=col_cluster, method=method, cmap=cmap,
                z_score=z_score, figsize=((len(pp.samples)), int(len(df_heatmap.index) / 4)), **options)
            self._export(big, pp.path, name="big-heatmap", export=export)
        elif len(df_heatmap) > 1000:
            print("Too many genes for a big heatmap")

    @staticmethod
    def _z_score(df):
        """Z-score calculation"""
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

    def lineplot(self, pp, title="", legend=True, z_score=False, export=False, df=None, size=10, ci=None, **option):
        """
        LinePlot selected genes expression levels. Max number of genes 200

        title - accept a str as title of the plot
        legend - True/False show the legend
        z_score - True/False calculate the z-score normalization
        export - True/False whether or not export the image
        df - accept an exernal dataframe, different from self.selected
        **options - all the options accepted by seaborn.factorplot"""

        if df is None:
            df = pp.df.copy()
            samples=pp.samples

        if z_score is True:
            df_ = self.onlyFPKM(df, samples, return_as="df", remove_FPKM_name=True)
            df_norm = self._z_score(df_)
            df_norm["gene_short_name"] = df["gene_short_name"]
            df_ = df_norm.copy()
        elif z_score is False:
            df_ = self.onlyFPKM(df, samples, return_as="gene name", remove_FPKM_name=True)

        print("Number of genes to plot: ", len(df_))
        if len(df_) > 50 and len(df_) < 200:
            print("Too many genes. Legend not shown")
            legend = False
        elif len(df_) >= 200:
            print("Too many genes. Plot not shown")
            return

        if pp.what == "gene":
            hue = "gene_short_name"
            df_ = self._fusion_gene_id(
                df_, pp.what, change_index=False)
        elif pp.what == "isoform":
            hue = "gene/ID"  # Change this hue name
            df_ = self._fusion_gene_id(
                df_, pp.what, change_index=True)
            df_ = df_.reset_index()
        df = pd.melt(df_, id_vars=[hue], var_name="Sample", value_name="FPKM")
        g = sns.factorplot(x="Sample", y="FPKM", hue=hue,
                           data=df, ci=ci, size=size, legend=legend, **option)
        g.fig.suptitle(title)
        self._export(g, pp.path, export=export, name="Plot")
        return g
    
    def _export(self, thing, path, export, name=None, image_extension=".png"):  # add .pdf?
        """Manage dataframe or image export parameter."""
        if export is False: 
            return
        elif export is True:
            _make_folder(path)
            thing.savefig(str(path + name + image_extension))
            print("\nExported as " + name + image_extension)
        else:
            raise Exception("export= can be only 'False' or 'True'")
