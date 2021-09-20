"""
This is a script for ...
"""
import os
import platform
import pandas as pd
import numpy as np
import functools
import warnings
import math
from scipy.stats import gmean
from statsmodels.stats.multitest import multipletests

from proteomics_tools._utils import FOLDER_DATA, SEP

# I Helper Functions
def get_organism(organism=None, out="GO"):
    """Get right organism for output"""
    list_out = ["Lower", "Upper", "UP", "GO"]
    if out not in list_out:
        raise ValueError("out must be one of {}".format(list_out))
    dict_org = {"HUMAN": ["HUMAN", "hsa", "human", "Homo sapiens (Human)"],
                "MOUSE": ["mmu", "mouse", "Mus musculus (Mouse)"],
                "RAT":  ["RAT", "rno", "rat", "Rattus norvegicus (Rat)"],
                "CHICK": ["CHICKEN", "gga", "chicken", "Gallus gallus (Chicken)"]}

    dict_org_out = {"HUMAN": {"GO": "hsa", "UP": "Homo sapiens (Human)", "Lower": "human", "Upper": "HUMAN"},
                    "MOUSE": {"GO": "mmu", "UP": "Mus musculus (Mouse)", "Lower": "mouse", "Upper": "MOUSE"},
                    "RAT":  {"GO": "rno", "UP": "Rattus norvegicus (Rat)", "Lower": "rat", "Upper": "RAT"},
                    "CHICK": {"GO": "gga", "UP": "Gallus gallus (Chicken)", "Lower": "chicken", "Upper": "CHICKEN"}}
    key = ""
    for org in dict_org:
        if organism in dict_org[org]:
            key = org
    return dict_org_out[key][out]


# Warning decorator
def ignore_warning(simplefilter=True, category=RuntimeWarning):
    """Ignore user warning just for function"""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Do something before
            if simplefilter:
                warnings.simplefilter("ignore", UserWarning)
            else:
                warnings.filterwarnings("ignore", category=category)
            # Func
            value = func(*args, **kwargs)
            # Do something after
            if simplefilter:
                warnings.simplefilter("default", UserWarning)
            else:
                warnings.filterwarnings("default", category=category)
            return value
        return wrapper
    return decorator


# II Main Functions
class Utils:
    """Helper class"""
    def __init__(self):
        pass

    # Common check method
    @staticmethod
    def check_p_method(p_method=None):
        """Check if p_method is valid"""
        list_p_method = ["bonferroni", "sidak", "holm", "fdr"]
        # Set p_method to all correction methods if None
        if p_method is None:
            p_method = list_p_method
        if type(p_method) is str:
            if p_method not in list_p_method:
                raise ValueError("'p_method' must be in {}".format(list_p_method))
            p_method = [p_method]
        else:
            for m in p_method:
                if m not in list_p_method:
                    raise ValueError("'p_method' must be in {}".format(list_p_method))
        return p_method

    # Data retrieval
    @staticmethod
    def keyword_df(df_id_entry=None, keyword="growth", column="DB_Object_Name"):
        """Get df for all GO terms which have keyword in given column"""
        df = df_id_entry[[True if x != 0 else False for x in df_id_entry[column].str.count(keyword)]]
        return df

    # Column/List converting
    @staticmethod
    def _get_ratio(df=None, col=None, as_values=False):
        """Convert str to ratio"""
        if as_values:
            return [[int(x.split("/")[0]), int(x.split("/")[1])] for x in df[col]]
        else:
            return [int(x.split("/")[0]) / int(x.split("/")[1]) for x in df[col]]

    @staticmethod
    def list_to_str(list_input):
        """Remove list syntax from list_input"""
        if type(list_input) is list:
            if len(list_input) > 0:
                return str(list_input).replace("[", "").replace("]", "").replace("'", "")
            else:
                return np.NAN
        elif type(list_input) == str and list_input == "":
            return np.NAN
        else:
            return list_input

    @staticmethod
    def set_to_str(set_input):
        """Remove list syntax from list_input"""
        if type(set_input) is set:
            if len(set_input) > 0:
                return str(set_input).replace("{", "").replace("}", "").replace("'", "")
            else:
                return np.NAN
        elif type(set_input) == str and set_input == "":
            return np.NAN
        else:
            return set_input

    @staticmethod
    def str_to_list(str_input, str_split=", "):
        """Convert string of a list into list"""
        if str_input is None or str_input is np.NAN or str(str_input) == "nan":
            return []
        if str_split in str_input:
            list_str = str_input.split(str_split)
        else:
            list_str = [str_input]
        return list_str

    def add_infos(self, row=None, dict_data=None, columns=None):
        """Add information from given columns"""
        for key in columns:
            if key in dict_data.keys():
                row.append(self.list_to_str(list(dict_data[key].keys())))
            else:
                row.append(np.nan)
        return row

    # Item retrieval
    @staticmethod
    def _get_items_from_col(df=None, col="study_items"):
        """Get all items from df column"""
        all_items = []
        for x in df[col].values:
            x = str(x)
            if x != "nan":
                if "," in x:
                    all_items.extend(x.split(", "))
                else:
                    all_items.append(x)
        all_items = list(set(all_items))
        return all_items

    @staticmethod
    def _get_items_from_dict(dict_key_list=None):
        """Get list with all keys and values"""
        all_items = list(dict_key_list.keys())
        for go_id in dict_key_list:
            all_items.extend(dict_key_list[go_id])
        all_items = list(set(all_items))
        return all_items

    # Enrichment score
    @staticmethod
    def _enrichment_score(df=None):
        """Get enrichment score, i.e. -log of geometric mean of uncorrected p values for all GO ids of one cluster
        Ref: Huang et al., 2008 (Systematic and integrative analysis of large gene lists
                using DAVID bioinformatics resources)"""
        p_vals = df["p_uncorrected"].values
        enrichment_score = round(-1.0 * math.log10(gmean(p_vals)), 2) # >  1.3 equivalent to p value of 0.05
        return enrichment_score

        # P value correction

    @staticmethod
    def correct_p_values(df=None, p_method=None):
        """Correct p values by given method or methods"""
        loc_p = df.columns.get_loc("p_uncorrected")
        dict_method = {"bonferroni": "bonferroni",
                       "sidak": "sidak",
                       "holm": "holm",
                       "fdr": "fdr_bh"}
        if type(p_method) is str:
            p_vals = multipletests(list(df["p_uncorrected"]), method=dict_method[p_method])[1]
            df.insert(loc_p + 1, "p_{}".format(p_method), p_vals)
        else:
            for i, p in enumerate(p_method):
                try:
                    p_vals = multipletests(list(df["p_uncorrected"]), method=dict_method[p])[1]
                    df.insert(loc_p + 1, "p_{}".format(p), p_vals)
                except ZeroDivisionError:
                    pass
        return df

    # Filtering
    @staticmethod
    def _get_p_col(p_method=None, df=None):
        """Get p column"""
        if p_method is None:
            p_method = "uncorrected"
        elif type(p_method) is list:
            p_method = p_method[0]
        if p_method in list(df):
            col_p_method = p_method
        elif "p_" + p_method in list(df):
            col_p_method = "p_" + p_method
        else:
            raise ValueError("'p_method' must be one of {}".format([x for x in list(df) if "p_" in x]))
        return col_p_method

    @staticmethod
    def _filter_ns(df_enrich=None, n=5):
        """Filter go enrichment results"""
        list_df = []
        for ns in ["CC", "MF", "BP"]:
            df = df_enrich[df_enrich["NS"] == ns].head(n)
            list_df.append(df)
        df_enrich = pd.concat(list_df)
        return df_enrich

    def filter_df(self, df=None, p_method=None, p_max=0.05, n=None):
        """Filter df based on p value or total number of highest ranked entries"""
        col_p_method = self._get_p_col(p_method=p_method, df=df)
        # Filter
        if "study_count" in list(df):
            df = df[(df[col_p_method] < p_max) & (df["study_count"] > 0)]
        else:
            df = df[df[col_p_method] < p_max]
        if n is not None:
            df = self._filter_ns(df_enrich=df, n=n)
        df = df[df["fold_enrichment"] > 1]
        df.sort_values(by="p_uncorrected", inplace=True, ascending=True)
        df.reset_index(drop=True, inplace=True)
        return df

