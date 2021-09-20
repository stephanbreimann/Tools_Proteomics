"""
This is a script for retrieving data from the cell surface protein atlas
Bausch-Fluck D, Hofmann A, Bock T, Frei AP, Cerciello F, et al. (2015)
A Mass Spectrometric-Derived Cell Surface Protein Atlas. PLoS One 10: e0121314.
"""
import pandas as pd

import proteomics_tools._utils as ut


# I Helper Function
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


def check_genes_in_df(list_genes=None, df=None, col=None):
    """"""
    missing_genes = []
    for gene in list_genes:
        if gene not in list(df[col]):
            missing_genes.append(gene)
    if len(missing_genes) > 0:
        raise ValueError("Following genes are not in Cell Surface Atlas: {}".format(missing_genes))


# II Main Functions
class CellSurfaceAtlas:
    """Class to retrieve data from the Cell Surface Protein Atlas (https://wlab.ethz.ch/cspa/)
    Bausch-Fluck D, Hofmann A, Bock T, Frei AP, Cerciello F, et al. (2015)
    A Mass Spectrometric-Derived Cell Surface Protein Atlas. PLoS One 10: e0121314."""
    def __init__(self, folder_in=None, organism="hsa"):
        if folder_in is None:
            folder_in = ut.FOLDER_DATA + "CellSurfaceAtlas" + ut.SEP
        self.folder_in = folder_in
        self.organism = organism
        self.dict_organism_go_up = {"hsa": "Human",
                                    "mmu": "Mouse"}
        self.dict_organism_up_go = {self.dict_organism_go_up[key]: key for key in self.dict_organism_go_up}
        self.df_cell_line = self._get_df_cell_line()

    # Helper methods
    def _check_organism(self, organism, up=True):
        """Check if chosen organism is suitable"""
        list_org = ["hsa", "mmu", "Human", "Mouse"]
        if organism not in list_org:
            raise ValueError("'organism' must be one of {}".format(list_org))
        if up:
            return self.dict_organism_go_up.get(organism, organism)
        else:
            return self.dict_organism_up_go.get(organism, organism)

    def _get_df_cell_line(self):
        """Get df with cell line information"""
        file = self.folder_in + "cell_surface_atlas_description.tab"
        df = pd.read_csv(file, sep="\t")
        return df

    def _get_df_org(self):
        """Load df of given organism"""
        organism = self._check_organism(self.organism, up=False)
        file_name = "cell_surface_atlas_{}".format(organism)
        file = self.folder_in + file_name + ".tab"
        df = pd.read_csv(file, sep="\t")
        return df

    # Main methods
    def cell_line_info(self, primary=None, filter_df=True):
        """Filter information of cell lines"""
        # Check input
        organism = self._check_organism(self.organism, up=True)
        df = self.df_cell_line.copy()
        df_filtered = df[(df["organism"] == organism)]
        if not filter_df:
            return df_filtered
        else:
            # Filter
            if primary is not None:
                if primary:
                    cell_line_primary = "primary"
                else:
                    cell_line_primary = "line"
                df_filtered = df_filtered[df_filtered["cell_line_or_primary"] == cell_line_primary]
            return df_filtered

    def get_proteins(self, cell_line="HEK"):
        """Get list of all expressed proteins for given cell line"""
        # Load data
        df = self._get_df_org()
        # Check input
        list_cell_lines = [x.upper().replace("-", "") for x in list(df) if x != "ACC"]
        df.columns = ["ACC"] + list_cell_lines
        cell_line = cell_line.upper()
        if "-" in cell_line:

            cell_line = cell_line.replace("-", "")
        if cell_line not in list_cell_lines:
            raise ValueError("'cell_line' must be one of {}".format(list_cell_lines))
        list_acc = df[df[cell_line] == 1]["ACC"].tolist()   # list with protein ac
        return list_acc

    def expression_surface(self, genes_target=None, genes_filter=None, acc=True):
        """Get cell lines where all given proteins are expressed"""
        genes_target, genes_filter = check_genes(genes_target=genes_target, genes_filter=genes_filter)
        cols = ["ACC", "Gene_Name"]
        if acc:
            col = "ACC"
        else:
            col = "Gene_Name"
        df = self._get_df_org()
        check_genes_in_df(list_genes=genes_target, df=df, col=col)
        df_filtered = df[df[col].isin(genes_filter)]
        # Filter
        list_cell_lines = [x for x in list(df_filtered.dropna(axis=1)) if x not in cols]
        list_col = cols + list_cell_lines
        df = df[df[col].isin(genes_target)][list_col]
        df.reset_index(drop=True, inplace=True)
        return df
