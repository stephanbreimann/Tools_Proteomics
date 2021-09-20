"""
This is a script for creating a dict with genes associated with neuronal cell types
    for Sharma et al., 2015 and Tüshaus et al., 2020
Sharma et al., 2015 (Cell type- and brain region-resolved mouse brain proteome)
    Results: https://www.nature.com/articles/nn.4160#Sec29 (Supplementary)
Tüshaus et al., 2020 (An optimized quantitative proteomcis method establishes the
    cell type-resolved mouse brain secretome)
    Results: https://www.embopress.org/doi/full/10.15252/embj.2020105693 (Supplementary EV4)
"""
import os
import pandas as pd
import numpy as np
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt

import proteomics_tools._utils as ut


# I Helper Functions
class _CellTypeFilter:
    """Class for filtering cell type specific (proteomic) datasets to use as references"""
    def __init__(self, df=None, col_fc="Log2 fold change", col_p="-Log10 p val",
                 col_id="Majority protein IDs", list_col_info=None, sep="_",
                 p_max=0.05, log2_fc_min=1):
        self.log2_fc_min = log2_fc_min
        self.p_max = p_max
        self.df = df.copy()
        self.sep = sep
        self.col_fc = col_fc
        self.col_p = col_p
        self.col_id = col_id
        # Basic info to add to output df
        if list_col_info is None:
            list_col_info = list(df)[0:2]
        self.list_col_info = list_col_info

    # Helper methods
    @staticmethod
    def _remove_from_str(str_in, list_remove=None):
        """Remove strings given in list from str_in"""
        for s in list_remove:
            if s in str_in:
                str_in = str_in.replace(s, "")
        return str_in

    @staticmethod
    def split_ids(dict_type_ids=None, sep=";"):
        """Split multiple major ids from lists of given dict"""
        dict_type_ids_split = {}
        for ct in dict_type_ids:
            list_ids = []
            for prot_id in dict_type_ids[ct]:
                if sep in prot_id:
                    list_ids.extend(prot_id.split(sep))
                else:
                    list_ids.append(prot_id)
            dict_type_ids_split[ct] = list_ids
        return dict_type_ids_split

    @staticmethod
    def filter_isoforms(dict_type_ids=None, marker="-"):
        """Replace protein ids of isoforms with protein of canonical form (remove if already exists)"""
        dict_type_ids_without_isoforms = {}
        for ct in dict_type_ids:
            list_ids = []
            for prot_id in dict_type_ids[ct]:
                if marker in prot_id:
                    list_ids.append(prot_id.split(marker)[0])
                else:
                    list_ids.append(prot_id)
            dict_type_ids_without_isoforms[ct] = list(set(list_ids))
        return dict_type_ids_without_isoforms

    def list_cell_types(self, match_str=None):
        """Get list with all cell types/conditions in data based on col with fold change"""
        if match_str is None:
            match_str = self.col_fc
        cell_type_list = [col.split("_")[1] for col in list(self.df) if match_str in col]
        return cell_type_list

    def plot_venn(self, list_names=None, list_set_ids=None, title=None):
        """Plot Venn diagram"""
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 3, 1)
        venn3(list_set_ids[0:2] + [list_set_ids[2]],
              list_names[0:2] + [list_names[2]], ax=ax)
        ax = fig.add_subplot(1, 3, 2)
        venn3(list_set_ids[0:2] + [list_set_ids[3]],
              list_names[0:2] + [list_names[3]], ax=ax)
        ax = fig.add_subplot(1, 3, 3)
        venn2(list_set_ids[2:], list_names[2:], ax=ax)
        if title is None:
            title = "Venn plot {} (p_max: {}, log2_fc_min: {})".format(self.col_id, self.p_max, self.log2_fc_min)
        plt.suptitle(title, size=12)
        return ax

    # Main methods
    def dict_type_df(self, list_cell_types=None):
        """Filter significant gene per cell type"""
        df = self.df
        dict_cell_df = {}
        min_log10_p = -np.log10(self.p_max)
        for cell_type in list_cell_types:
            # Get columns for filtering
            col_change = self.col_fc + self.sep + cell_type
            col_p_val = self.col_p + self.sep + cell_type
            list_cell_type_col = [col for col in list(df) if cell_type in col]
            list_col = [self.col_id] + self.list_col_info + list_cell_type_col
            # Filter
            df_filtered = df[(df[col_change] >= self.log2_fc_min) & (df[col_p_val] >= min_log10_p)]
            dict_cell_df[cell_type] = df_filtered[list_col]
        return dict_cell_df

    def dict_type_ids(self, dict_type_df=None, plot_venn=False):
        """Get dict for cell types to list of major ids"""
        list_id_sets = [set(dict_type_df[ct][self.col_id]) for ct in dict_type_df]
        # Save ids for each cell type (including overlapping ids)
        dict_type_ids = {ct: list(dict_type_df[ct][self.col_id]) for ct in dict_type_df}
        list_names = list(dict_type_df.keys())
        if plot_venn:
            self.plot_venn(list_names=list_names,
                           list_set_ids=list_id_sets)
        return dict_type_ids

    @staticmethod
    def dict_id_type(dict_type_ids=None):
        """Get dict with separated major ids to significant up regulated cell type"""
        # Save cell type for each id (excluding overlapping ids)
        dict_id_type = {}
        for ct_a in dict_type_ids:
            list_ids_other_types = []
            # Get list with all ids of other cell types
            for ct_b in dict_type_ids:
                if ct_b != ct_a:
                    list_ids_other_types.extend(dict_type_ids[ct_b])
            dict_id_type.update({i: ct_a for i in dict_type_ids[ct_a]
                                 if i not in list_ids_other_types})
        return dict_id_type


# II Main Functions
class ProteomicsReferences:
    """Class for retrieving reference data sets"""

    @staticmethod
    def dict_sharma(p_max=0.05, log2_fc_min=2, key_genes=False, plot_venn=False,
                    split_ids=True, filter_isoforms=True):
        """Get dictionary with enriched (not specified) genes of cell type. Show as venn diagram.
        Ref.: Sharma et al., 2015 (Cell type- and brain region-resolved mouse brain proteome)"""
        file = ut.FOLDER_DATA + "Proteomics_References" + ut.SEP + "Sharma_neuro.xlsx"
        df = pd.read_excel(file)
        ctf = _CellTypeFilter(df=df, p_max=p_max, log2_fc_min=log2_fc_min)
        list_cell_types = ctf.list_cell_types()
        dict_type_df = ctf.dict_type_df(list_cell_types=list_cell_types)
        dict_type_ids = ctf.dict_type_ids(dict_type_df=dict_type_df, plot_venn=plot_venn)
        if split_ids:
            dict_type_ids = ctf.split_ids(dict_type_ids=dict_type_ids, sep=";")
        if filter_isoforms:
            dict_type_ids = ctf.filter_isoforms(dict_type_ids=dict_type_ids)
        dict_id_type = ctf.dict_id_type(dict_type_ids=dict_type_ids)
        if key_genes:
            return dict_id_type
        else:
            return dict_type_ids

    @staticmethod
    def dict_tueshaus(key_genes=True, plot_venn=False):
        """Results of mouse brain secretome as dictionary
        Ref.: Tüshaus et al., 2020 (An optimized quantitative proteomcis method establishes the
        cell type-resolved mouse brain secretome)"""
        file = ut.FOLDER_DATA + "Proteomics_References" + ut.SEP + "Tüshaus_neuro_secretome.xlsx"
        df = pd.read_excel(file)
        col_id = "uniprot_acc"
        ctf = _CellTypeFilter(df=df, col_id=col_id)
        list_cell_types = ctf.list_cell_types(match_str="specific")
        dict_id_type = {}
        dict_type_ids = {}
        for cell_type in list_cell_types:
            col_type = "specific_"
            col = col_type + cell_type
            list_up_ids = list(df[df[col] == "+"][col_id])
            # Modification of cell types to match with Sharma
            cell_type = cell_type.capitalize()
            if cell_type == "Neuron":
                cell_type += "s"
            dict_id_type.update({up_id: cell_type for up_id in list_up_ids})
            dict_type_ids[cell_type] = list_up_ids
        list_names = list(dict_type_ids.keys())
        list_set_ids = [set(x) for x in dict_type_ids.values()]
        if plot_venn:
            title = "Venn plot protein IDs"
            ctf.plot_venn(list_names=list_names,
                          list_set_ids=list_set_ids,
                          title=title)
        if key_genes:
            return dict_id_type
        else:
            return dict_type_ids

    @staticmethod
    def expression_cell_types(acc_target=None, dict_ref=None):
        """"""
        # TODO
        list_ct = sorted(dict_ref.keys())
        list_hits = []
        for acc in acc_target:
            hits = []
            for ct in dict_ref:
                print(ct, len(dict_ref[ct]))
                if acc in dict_ref[ct]:
                    hits.append(1)
                else:
                    hits.append(0)
            list_hits.append(hits)
        return list_hits
