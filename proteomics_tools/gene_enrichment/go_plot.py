"""
This is a script for plotting enrichment results
"""
import seaborn as sns
from matplotlib import pyplot as plt

import proteomics_tools.gene_enrichment._utils as ut


# II Main Functions
class PlotGo:
    """Class for plotting enrichment results"""
    def __init__(self, fig_format="png", verbose=False):
        self.fig_format = fig_format
        self.verbose = verbose

    # Helper methods
    @staticmethod
    def join_folder_title(title_out=None, folder_out=None):
        """Join title and folder"""
        dict_rep = {" ": "_",
                    "(": "",
                    ")": ""}
        for rep_str in dict_rep:
            if rep_str in title_out:
                title_out = title_out.replace(rep_str, dict_rep[rep_str])
        if folder_out is not None:
            out = folder_out + ut.SEP + title_out
        else:
            out = title_out
        return out

    @staticmethod
    def _shorten_terms(term=None, lim=30):
        """"""
        dict_shorter = {"negative": "neg.",
                        "positive": "pos.",
                        "activity": "act.",
                        "system": "sys.",
                        "regulation": "reg.",
                        "process": "proc.",
                        "nucleotide": "nt",
                        "protein": "prot.",
                        "multicellular": "multicel.",
                        "macromolecule": "macromol.",
                        "structural": "struct.",
                        "extracellular matrix": "ECM",
                        "inhibitor": "inhib.",
                        "molecular": "mol.",
                        "endoplasmic reticulum": "ER"}
        for shorter in dict_shorter:
            if shorter in term and len(term) > lim:
                term = term.replace(shorter, dict_shorter[shorter])
        return term

    def _split_name(self, df=None, y="NAME"):
        """"""
        names = []
        for n in list(df[y]):
            len_n = len(n)
            if " - " in n:
                n_splits = n.split(" - ")
                new_n = n_splits[0]
            elif len_n > 30:
                n = self._shorten_terms(term=n)
                if len(n) > 40:
                    n_splits = n.split(" ")
                    middle = int(len(n_splits)/2) + 1
                    new_n = " ".join(n_splits[0:middle]) + "\n" + " ".join(n_splits[middle:])
                else:
                    new_n = n
            else:
                new_n = n
            names.append(new_n)
        df[y] = names
        return df

    @staticmethod
    def _get_significance_label(df=None, p_col=None):
        """Map p values to R significance symbols"""
        # https://www.rdocumentation.org/packages/gtools/versions/3.5.0/topics/stars.pval
        p_values = list(df[p_col])
        list_stars = []
        for p in p_values:
            if p <= 0.001:
                list_stars.append("***")
            elif p <= 0.01:
                list_stars.append("**")
            elif p <= 0.05:
                list_stars.append("*")
            else:
                list_stars.append("")
        return list_stars

    @staticmethod
    def _get_fontsize(mode="long"):
        """Get fontsize based on length of table"""
        if mode == "short":
            size = 12
        else:
            size = 11
        return size

    @staticmethod
    def _get_add_y(significance_label=None):
        """Get y based on length of table"""
        add_y = len(significance_label) / (len(significance_label) * 1.25)
        return add_y

    @staticmethod
    def _get_x(x_max=10):
        """Get x value based on maximum"""
        if x_max < 10:
            x = 0.1
        elif x_max < 50:
            x = 0.5
        elif x_max < 100:
            x = 1
        else:
            x = 5
        return x

    def _add_significance_symbols(self, ax=None, rects=None, significance_label=None, x_max=10, mode=None):
        """Add significance symbols to plot"""
        # https://stackoverflow.com/questions/63521320/how-to-add-text-values-in-bar-plot-seaborn-python
        # https://stackoverflow.com/questions/37579284/matplotlib-get-coordinates-of-top-and-bottom-of-horizontal-bar-edges
        x = self._get_x(x_max=x_max)
        add_y = self._get_add_y(significance_label=significance_label)
        size = self._get_fontsize(mode=mode)
        i = 0
        for rect in rects:
            if str(rect.get_width()) != "nan":
                y = rect.get_y() + add_y
                ax.text(x, y, significance_label[i],
                        ha='left', va='bottom', rotation=0, color='lightgray', fontweight="bold", size=size)
                i += 1

    @staticmethod
    def _get_df(df=None, y=None, x=None, hue=None, p_col=None):
        """Get df based on arguments"""
        if p_col in list(df):
            df = df[[hue, y, x, p_col]]
        elif p_col is None:
            df = df[[hue, y, x]]
        else:
            raise ValueError("'p_col' must be one of {}".format([x for x in list(df) if "p_" in x]))
        return df

    @staticmethod
    def _get_mode(df=None):
        """Get short or long mode for figure"""
        if len(df) > 10:
            mode = "long"
        else:
            mode = "short"
        return mode

    def enrichment_plot(self, df=None, hue="NS", y="NAME", x="fold_enrichment", title="Fold Enrichment",
                        split_name=False, p_col=None):
        """Plot fold enrichment or enrichment score"""
        df = self._get_df(df=df, x=x, y=y, hue=hue, p_col=p_col)
        mode = self._get_mode(df=df)
        if len(df) > 0:
            sns.set_theme(style="whitegrid")
            x_max = max(df[x])
            df = df.sort_values(by=[hue, x], ascending=False)
            if split_name:
                df = self._split_name(df=df, y=y)
            if mode == "long":
                fig, ax = plt.subplots(figsize=(8, 6))
            else:
                fig, ax = plt.subplots(figsize=(8, 4))
            sns.barplot(data=df, y=y, x=x, hue=hue, orient="h", dodge=False, ci=None, ax=ax)
            plt.xlim(0, int(x_max * 1.1) + 1)   # Increase right limit of plot by 20%
            plt.legend(loc=4)
            plt.title(title, fontdict={"fontweight": "bold", "fontsize": 12})
            plt.ylabel("")
            plt.xlabel(x.replace("_", " ").capitalize(), size=12)
            if p_col is not None and p_col in list(df):
                significance_label = self._get_significance_label(df=df, p_col=p_col)
                self._add_significance_symbols(ax=ax,
                                               rects=ax.patches,
                                               significance_label=significance_label,
                                               x_max=x_max,
                                               mode=mode)
        else:
            if self.verbose:
                print("No entries to show in enrichment plot")

