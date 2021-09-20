"""
This is a script for retrieving data from Proteomics DB
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

import proteomics_tools._utils as ut


# I Helper Functions
def check_input_file(input_file=None, list_files=None):
    """Check if input file exists and if has .tsv"""
    list_files = [x for x in list_files]
    if input_file not in list_files:
        raise TypeError("'{}' does not exists. Chose from following: {}. "
                        "\n or retrieve expression data from ProteomicsDB: "
                        "https://www.proteomicsdb.org/proteomicsdb/#overview "
                        "".format(input_file, list_files))
    return input_file


def check_path(path):
    """"""
    if not os.path.exists(path):
        raise ValueError(f"Following path does not exist: {path}")


def check_genes(genes_target=None, genes_filter=None):
    """"""
    if genes_target is None:
        raise TypeError("'list_genes' should not be None")
    if type(genes_target) is str:
        genes_target = [genes_target]
    if type(genes_filter) is str:
        genes_filter = [genes_filter]
    if genes_filter is None:
        genes_filter = genes_target
    else:
        genes_target = genes_target + [gene for gene in genes_filter if gene not in genes_target]
    return genes_target, genes_filter


# II Main Functions
class ProteomicsDB:
    """Class to retrieve downloaded data from Human Protein Atlas https://www.proteinatlas.org/
    Expression can not be downloaded programmatically!"""

    def __init__(self, folder_in=None):
        if folder_in is None:
            folder_in = ut.FOLDER_DATA #.folder #os.path.dirname(os.path.realpath(__file__)) + "/data/ProteomicsDB/"
        self._folder_in = folder_in + "ProteomicsDB" + ut.SEP
        check_path(path=self._folder_in)

    # Helper methods
    def load_file(self, biological_source="Cell_Line", gene=None):
        """Load input file"""
        path = self._folder_in + biological_source + ut.SEP
        input_file = f"Protein_Expression_{gene}.csv"
        list_files = [x for x in os.listdir(path) if "#" not in x]
        input_file = check_input_file(input_file=input_file,
                                      list_files=list_files)
        df = pd.read_csv(path + input_file, sep=";")
        return df

    @staticmethod
    def heatmap(df, title=None, annot=True, max_annot=20, biological_source="Cell_Line"):
        """Plot heatmap for df_nx"""
        if len(df) == 0:
            raise ValueError("Empty df")
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
        y_label = biological_source.replace("_", " ") if "_" in biological_source else biological_source
        plt.ylabel(y_label, size=12, weight="bold")
        if title is None:
            title = "Protein Expression [log10 normalized iBAO intensity]"
            plt.title(title)
        return ax

    def expression_proteins(self, genes_target=None, genes_filter=None, min_exp=1, biological_source="Cell_Line"):
        """Get protein expression for all tissues. Expression is given for all target genes and filtered
        based on 'min_exp' threshold for list 'gene_filter'"""
        genes_target, genes_filter = check_genes(genes_target=genes_target, genes_filter=genes_filter)
        cols = ["Tissue", "Average Normalized Intensity"]
        dict_gene_df = {gene: self.load_file(gene=gene, biological_source=biological_source)[cols]
                        for gene in genes_target}
        list_df = []
        for gene in dict_gene_df:
            df = dict_gene_df[gene]
            df.set_index("Tissue", inplace=True)
            list_df.append(df)
        df_res = pd.concat(list_df, axis=1)
        df_res.columns = genes_target
        df_res.fillna(0, inplace=True)
        # Filter
        mask = (df_res[genes_filter] >= min_exp).all(axis=1)
        df_res = df_res[mask]
        return df_res
