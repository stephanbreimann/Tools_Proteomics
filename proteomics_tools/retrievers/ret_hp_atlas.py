"""
This is a script for comparing different genes using Human Protein Atlas
    or to get summary statistics of target genes
"""
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

import proteomics_tools._utils as ut


# I Helper Functions
def _add_hek_cells(df_cell, list_genes=None, folder_in=None, filter_col="NX"):
    """Add HEK cells to cell line df"""
    file = "rna_celline.tsv"
    df = pd.read_csv(folder_in + file, sep="\t")
    df = df[(df["Gene name"].isin(list_genes)) & (df["Cell line"] == "HEK 293")]
    df = df.set_index("Gene name").rename({filter_col: "HEK 293"}, axis=1)["HEK 293"]
    df_hek = df_cell.transpose().join(df)
    df_hek = df_hek[["HEK 293"] + list(df_hek)[0:-1]].sort_index()
    return df_hek.transpose()


# II Main Functions
class HumanProteinAtlas:
    """Class to retrieve data from Human Protein Atlas https://www.proteinatlas.org/"""

    def __init__(self, folder_in=None):
        if folder_in is None:
            folder_in = ut.FOLDER_DATA + "HumanProteinAtlas/"
        self._folder_in = folder_in
        # File overview
        self.list_files = [x for x in os.listdir(self._folder_in) if "#" not in x]
        self.list_tissue_cell_type = ["rna_single_cell_type_tissue.tsv", "normal_tissue.tsv", "proteinatlas.tsv"]
        self.all_ran_expression_files = [x for x in self.list_files if x not in self.list_tissue_cell_type]
        # Data overview
        self.list_target_col = ["Tissue", "Cell type", "Blood cell", "Brain region", "Cell line"]
        self.levels = ["Not detected", "Low", "Medium", "High"]
        self._dict_level_str_int = {level: i for i, level in enumerate(self.levels)}
        self._dict_level_int_str = {self._dict_level_str_int[i]: i for i in self._dict_level_str_int}
        self._dict_level_int_str[0] = "ND"

    # Helper methods
    def _check_input_file(self, input_file):
        """Check if input file exists and if has .tsv"""
        if "tsv" not in input_file:
            input_file += ".tsv"
        if input_file not in self.list_files:
            raise TypeError("'input_file' does not exists. Chose one of following: {}".format(self.list_files))
        return input_file

    @staticmethod
    def _check_if_col_in_df(df, col):
        """Check if target column in df"""
        if col not in list(df):
            raise TypeError("{} not in df. Chose one of following: {}".format(col, list(df)))

    def _check_target_col(self, df, target_col, verbose=False):
        """Check if target column in df or select matching column"""
        if target_col not in list(df):
            for tc in self.list_target_col:
                if tc in list(df):
                    target_col = tc
                    break
            if verbose:
                print("'{}' was selected as 'target_col'".format(target_col))
        return target_col
    
    @staticmethod
    def _check_genes(list_genes=None, co_expressed_genes=None):
        """Check if gene input data"""
        if list_genes is None:
            raise ValueError("'list_genes' must be given")
        if co_expressed_genes is None:
            co_expressed_genes = list_genes
        list_genes = list(dict.fromkeys(list_genes))  # Remove duplicates while preserving order
        co_expressed_genes = list(dict.fromkeys(co_expressed_genes))
        list_genes = co_expressed_genes + [g for g in list_genes if g not in co_expressed_genes]
        return list_genes, co_expressed_genes
    
    # Filtering methods
    def _filter_rna_expression(self, input_file="rna_celline.tsv", list_genes=None, co_expressed_genes=None,
                               target_col="Cell line", filter_col="NX", min_rna=0.5, verbose=False):
        """Filter RNA expression based on 'min_rna' threshold for co expressed genes"""
        # Check input
        list_genes, co_expressed_genes = self._check_genes(list_genes=list_genes,
                                                           co_expressed_genes=co_expressed_genes)
        df = self.load_file(input_file)
        self._check_if_col_in_df(df, filter_col)
        target_col = self._check_target_col(df, target_col, verbose=verbose)
        list_filter_col = ["TPM", "pTPM", "NX"]
        if filter_col not in list_filter_col:
            raise TypeError("'filter_col' must be one of {}".format(list_filter_col))
        # Filtering
        list_cl = []
        for gene in co_expressed_genes:
            df_filter = df[(df["Gene name"] == gene) & (df[filter_col] >= min_rna)]
            list_cl.append(set(df_filter[target_col]))
        list_matches = list(set.intersection(*list_cl))
        if len(list_matches) == 0:
            raise ValueError("No co-expression in {} for {}".format(input_file, co_expressed_genes))
        list_df = [df[(df["Gene name"] == gene)] for gene in list_genes]
        list_df_filtered = []
        for df in list_df:
            df_filtered = df[df[target_col].isin(list_matches)][[target_col, filter_col]]
            df_filtered.set_index(target_col, inplace=True)
            list_df_filtered.append(df_filtered)
        df_res = pd.concat(list_df_filtered, axis=1)
        df_res.columns = list_genes
        return df_res

    def _filter_normal_tissue(self, input_file=None, list_genes=None, co_expressed_genes=None, min_level="Low", int_out=True):
        """Filter based on Protein level for normal_tissue.tsv"""
        list_genes, co_expressed_genes = self._check_genes(list_genes=list_genes,
                                                           co_expressed_genes=co_expressed_genes)
        df = self.load_file(input_file)
        self._check_if_col_in_df(df, "Level")
        levels = self._dict_level_str_int.keys()
        if min_level not in levels:
            raise TypeError("'min_level' must be one of {}".format(levels))
        df["Level"] = [self._dict_level_str_int.get(l, 0) for l in df["Level"]]
        # Filtering based on co-expressed genes
        list_cl = []
        for gene in co_expressed_genes:
            df_filter = df[(df["Gene name"] == gene) & (df["Level"] >= self._dict_level_str_int[min_level])]
            list_cl.append(set(zip(df_filter["Tissue"], df_filter["Cell type"])))
        # list of df for all genes
        list_df = [df[(df["Gene name"] == gene)] for gene in list_genes]
        list_matches = sorted(list(set.intersection(*list_cl)))
        # Join filtered data
        list_data = []
        for df in list_df:
            list_levels = []
            for i in list_matches:
                level_int = df[(df["Tissue"] == i[0]) & (df["Cell type"] == i[1])]["Level"]
                if int_out:
                    list_levels.extend(level_int)
                else:
                    list_levels.extend([self._dict_level_int_str[i] for i in level_int])
            list_data.append(list_levels)
        df_res = pd.DataFrame(data=list_data, index=list_genes).transpose()
        df_res.insert(0, "Cell type", [i[1] for i in list_matches])
        df_res.insert(0, "Tissue", [i[0] for i in list_matches])
        df_res.fillna(0, inplace=True)
        return df_res

    def _filter_rna_single_cell_type_tissue(self, input_file=None, list_genes=None, co_expressed_genes=None, min_ptpm=5):
        """Filter ran_single_cell_type_tissue file using pTPM threshold"""
        list_genes, co_expressed_genes = self._check_genes(list_genes=list_genes,
                                                           co_expressed_genes=co_expressed_genes)
        df = self.load_file(input_file)
        self._check_if_col_in_df(df, "pTPM")
        # Filtering
        list_cl = []
        list_df = []
        for gene in co_expressed_genes:
            df_filter = df[(df["Gene name"] == gene) & (df["pTPM"] >= min_ptpm)]
            list_df.append(df_filter)
            list_cl.append(set(zip(df_filter["Tissue"], df_filter["Cell type"])))
        list_matches = sorted(list(set.intersection(*list_cl)))
        list_df = [df[(df["Gene name"] == gene)] for gene in list_genes]
        # Join filtered data
        list_data = []
        for df in list_df:
            list_exp = []
            for i in list_matches:
                # Calculate ptpm mean to account for multiple clusters
                ptpm_mean = df[(df["Tissue"] == i[0]) & (df["Cell type"] == i[1])]["pTPM"].mean()
                list_exp.append(ptpm_mean)
            list_data.append(list_exp)
        df_res = pd.DataFrame(data=list_data, index=list_genes).transpose()
        df_res.insert(0, "Cell type", [i[1] for i in list_matches])
        df_res.insert(0, "Tissue", [i[0] for i in list_matches])
        return df_res
    
    # Plotting Heatmap
    @staticmethod
    def _get_title(file, filter_col="NX"):
        """Set title for RNA expression"""
        if ".tsv" in file:
            file = file.replace(".tsv", "")
        title = "Expression for " + file.replace("rna", "RNA").replace("_", " ") + " [{}]".format(filter_col)
        # FANTOM5: Functional annotation of the mammalian genome https://fantom.gsc.riken.jp/5/
        # Via Cap Analysis of Gene Expression (CAGE)
        # GTEx: Genotype-Tissue Expression https://gtexportal.org/home/
        # collected from 54 non-diseased tissue sites across nearly 1000 individual (WGS, WES, RNA-Seq)
        # HPA: Human Protein Atlas via RNA-seq
        dict_sources = {"fantom": "FANTOM5",
                        "gtex": "GTEx",
                        "hpa": "HPA"}
        for key in dict_sources:
            if key in title:
                title = title.replace(key, dict_sources[key])
        return title

    @staticmethod
    def _heatmap_expression(df, title=None, annot=True, max_annot=20):
        """Plot heatmap for df_nx"""
        if len(df) == 0:
            raise TypeError("Empty df")
        fmt = ".3g"
        if type(annot) is bool:
            if annot and len(df) >= max_annot:
                annot = False
        else:
            fmt = ""
        plt.tight_layout()
        sns.set_style("whitegrid")
        # colors: https://chrisalbon.com/python/data_visualization/seaborn_color_palettes/
        cmap = sns.diverging_palette(220, 20, as_cmap=True)
        ax = sns.heatmap(df, annot=annot, cmap=cmap, xticklabels=True, yticklabels=True, fmt=fmt)
        plt.title(title, size=11, weight="bold")
        y_label = df.index.name
        plt.ylabel(y_label, size=12, weight="bold")
        return ax
    
    # Main methods
    def file_info(self):
        """Get overview of all available files"""
        list_keys = ["tissue", "single_cell",
                     "brain", "blood",
                     "consensu", "fantom", "gtex", "hpa",
                     "celline", "proteinatlas"]
        dict_files = {key: sorted([x for x in self.list_files if key in x]) for key in list_keys}
        for key in dict_files["tissue"]:
            if "single_cell" in key:
                dict_files["tissue"].remove(key)
        return dict_files

    def load_file(self, input_file):
        """Load input file"""
        input_file = self._check_input_file(input_file)
        df = pd.read_csv(self._folder_in + input_file, sep="\t")
        return df
    
    def rna_expression(self, file="rna_consensus", list_genes=None, co_expressed_genes=None, filter_col="NX",
                       min_rna=0.1, min_rna_step=0.1, max_n=20, heatmap=False, annot=True, add_hek=True):
        """Get RNA expression for list of genes in given input file
        In: a) file: input file downloaded from Human Protein Atlas
            b) list_genes: list of genes for which expression is shown
            c) co_expressed_genes: list of genes that are used for filtering
            d) filter_col: filtering column {'NX', 'pTPM'}
            e1) min_rna: minimum threshold for filtering based on filter_col
            e2) max_n: maximum of selected features
            f1) heatmap: boolean whether heatmap should be created
            f2) annot: boolean whether heatmap should be labeled with values
        Out:a) df_res: result df with comparison
        Notes
        -----
            'NX': normalized RNA expression
            'pTPM': RNA level (protein coding transcript per million)'
        """
        kwargs = {"input_file": file, "list_genes": list_genes, "filter_col": filter_col,
                  "co_expressed_genes": co_expressed_genes}

        df_res = self._filter_rna_expression(**kwargs, min_rna=min_rna)
        while len(df_res) > max_n:
            min_rna += min_rna_step
            df_res = self._filter_rna_expression(**kwargs, min_rna=min_rna)
        if len(df_res) != 0:
            if "rna_celline" in file and add_hek:
                df_res = _add_hek_cells(df_res,
                                        folder_in=self._folder_in,
                                        list_genes=list_genes,
                                        filter_col=filter_col)
            file_out = file.replace(".tsv", "")
            file_out += "_" + "_".join(co_expressed_genes)
            if heatmap:
                title = self._get_title(file=file, filter_col=filter_col)
                ax = self._heatmap_expression(df_res, title=title, annot=annot, max_annot=30)
                return df_res, ax
            return df_res

    def all_rna_expression(self, list_genes=None, co_expressed_genes=None, files_excluded=None,
                           min_rna=0.1, max_n=25, annot=True):
        """Plot all RNA expression files that have not tissue and cell type (c.f. rna_expression)"""
        MIN_RNA = min_rna   # Constant input min_rna
        if files_excluded is None:
            files_excluded = self.list_tissue_cell_type
        elif type(files_excluded) is str:
            files_excluded = [files_excluded] + self.list_tissue_cell_type
        else:
            files_excluded += self.list_tissue_cell_type
        list_files = [x for x in self.list_files if x not in files_excluded]
        list_data = []
        for file in list_files:
            data = self.rna_expression(file=file,
                                       list_genes=list_genes,
                                       co_expressed_genes=co_expressed_genes,
                                       min_rna=MIN_RNA,
                                       max_n=max_n,
                                       annot=annot,
                                       heatmap=True)

    def tissue_type_expression(self, list_genes=None, co_expressed_genes=None, level="Protein", min_th=None,
                               heatmap=True, annot=True):
        """Heatmap for normal tissue or single sequence RNA data for tissue and cell type
        In: a) list_genes: list of genes for which expression is shown
            b) co_expressed_genes: list of genes that are used for filtering
            d) level: string to decide for protein or single cell RNA level {'Protein', 'ssRNA'}
            e1) min_th: minimum threshold for filtering for level resp. pTPM int
            f1) heatmap: boolean whether heatmap should be created
            f2) annot: boolean whether heatmap should be labeled with values
        Out:a) df_res: result df with comparison"""
        list_level = ["ssRNA", "Protein"]
        if level not in list_level:
            raise TypeError("'level' should be one of following: {}".format(list_level))
        kwargs = {"list_genes": list_genes,
                  "co_expressed_genes": co_expressed_genes}
        if level == "Protein":
            if min_th is None:
                min_th = "Low"
            kwargs.update({"min_level": min_th,
                           "input_file": "normal_tissue.tsv"})
            df_res = self._filter_normal_tissue(**kwargs)
            df_res.set_index(["Tissue", "Cell type"], inplace=True)
            file_out = "normal_tissue"
            title = "Level of protein expression in " + file_out.replace("_", " ")
            annot = df_res.applymap(lambda x: self._dict_level_int_str.get(x, 0)).values
        else:
            if min_th is not None:
                kwargs.update({"min_ptpm": min_th})
            kwargs.update({"input_file": "rna_single_cell_type_tissue.tsv"})
            df_res = self._filter_rna_single_cell_type_tissue(**kwargs)
            df_res.set_index(["Tissue", "Cell type"], inplace=True)
            file_out = "rna_single_cell_type_tissue"
            title = self._get_title(file=file_out, filter_col="pTPM")
        if heatmap:
            plt.figure(figsize=(7, 4))
            ax = self._heatmap_expression(df_res, title=title, annot=annot, max_annot=30)
            if level == "Protein":
                colorbar = ax.collections[0].colorbar
                start, stop = int(df_res.min().min()), int(df_res.max().max())
                levels = self.levels[start:stop+1]
                ticks = list(range(start, stop + 1))
                colorbar.set_ticks(ticks)
                colorbar.set_ticklabels(levels)
            plt.ylabel("Tissue|Cell type")
            return df_res, ax
        return df_res


