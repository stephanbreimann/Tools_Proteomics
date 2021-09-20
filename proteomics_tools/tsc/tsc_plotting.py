"""
This is a script for plotting of time series clustering results
"""
import math
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib_venn import venn3
from matplotlib.lines import Line2D

from proteomics_tools.tsc import TimeSeriesClustering
from proteomics_tools.retrievers import ProteomicsReferences

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# TODO Refactor (make new)
# I Helper Functions
def get_colors(list_labels=None):
    """"""
    colors = ["tab:blue", "tab:orange", "tab:green",  "tab:gray",  "tab:olive",
              "tab:purple", "tab:brown",  "tab:cyan", "tab:red", "tab:pink"]
    unique_labels = list(set(list_labels))
    dict_label_colors = {label: color for label, color in zip(unique_labels, colors[0:len(unique_labels)])}
    list_colors = [dict_label_colors[label] for label in list_labels]
    return list_colors


class TSCPlottingSettings:
    """Settings for plotting"""
    def __init__(self):
        self._dict_nrows_x = {2: 0.95, 3: 0.925, 4: 0.90, 5: 0.875, 6: 0.85, 7: 0.8}
        self._dict_nrows_legend = {2: 20, 3: 18, 4: 16, 5: 8, 6: 4, 7: 2}
        self._dict_nrows_size = {2: 11, 3: 11, 4: 11, 5: 9, 6: 9, 7: 7}
        # https://matplotlib.org/3.1.0/gallery/color/named_colors.html
        self.colors_a = ['purple', 'green', 'steelblue', "orange", "cyan",
                         "gray",  "lightblue", "yellow", "salmon",  "red",
                         "peru", "navy", "gold", "orchid", "skyblue"]
        self.colors_b = ["tab:blue", "tab:orange", "tab:green",  "tab:gray",  "tab:olive",
                         "tab:purple", "tab:brown",  "tab:cyan", "tab:red", "tab:pink"]
        self.dict_colors_a = {i: color for i, color in enumerate(self.colors_a)}
        self.dict_colors_b = {i: color for i, color in enumerate(self.colors_b)}

    @staticmethod
    def _plot_settings():
        """Settings for plotting"""
        plt.tight_layout()
        plt.style.use("default")

    @staticmethod
    def nrows_ncol(n_clusters=None):
        """Get n rows and n columns"""
        nrows = math.ceil(n_clusters / 2) + 1
        ncol = 2
        return nrows, ncol

    def _set_min_max(self, df=None):
        """Set valuesfor plotting"""
        self.min_val = min(df.min()) - 0.5
        self.max_val = max(df.max()) + 0.5


class TSCPlottingData(TSCPlottingSettings):
    """Class for organizing data used in tsc"""
    def __init__(self, sharma_max_p=0.05, sharma_min_log2_fold_change=1):
        TSCPlottingSettings.__init__(self)
        # Sharma dict for plotting
        pr = ProteomicsReferences()
        self.dict_sharma = pr.dict_sharma(p_max=sharma_max_p,
                                          log2_fc_min=sharma_min_log2_fold_change,
                                          key_genes=True)
        self.dict_sharma_colors = self._dict_sharma_colors()

    # Helper methods
    def sharma_labels(self, list_genes=None):
        """Get list of labels from sharma df"""
        if list_genes is None:
            list_cell_types = self.dict_sharma.values()
        else:
            list_cell_types = []
            for gene in list_genes:
                if gene in self.dict_sharma.keys():
                    list_cell_types.append(self.dict_sharma[gene])
        list_labels = list(set(list_cell_types)) + ["Not Enriched"]
        list_labels.sort()
        return list_labels

    def sharma_legend_colors(self, list_ids=None):
        """Get list of colors for legend"""
        return [self.dict_sharma_colors[label] for label in self.sharma_labels(list_genes=list_ids)]

    def _dict_sharma_colors(self):
        """Get dict for sharma cell types to colors"""
        dict_colors = dict(zip(self.sharma_labels(list_genes=None), self.colors_b))
        return dict_colors


# II Main Functions
class TSCSinglePlots(TSCPlottingSettings):
    """Class for evaluation and optimization of time series clustering"""
    def __init__(self):
        TSCPlottingSettings.__init__(self)

    # Time series plot
    def time_series_plot(self, df_val=None, title=None, labels=None, x_labels=None):
        """Line plot for given df"""
        self._plot_settings()
        self._set_min_max(df=df_val)
        font_arg = {"fontweight": "bold", "fontsize": 12}
        df_val.transpose().plot(kind="line", color=labels)
        plt.style.use("seaborn")
        plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
        plt.ylim(self.min_val, self.max_val)
        # Labels
        plt.title(title, **font_arg)
        plt.ylabel("log2 ratio of mean LFQs", **font_arg)
        plt.xlabel("Time point", **font_arg)
        plt.xticks(ticks=range(0, len(x_labels)), labels=x_labels)

    # PCA plot
    @staticmethod
    def pca_plot(df=None, group_label="labels", dict_colors=None, reverse=False,
                 title=None, fig=None, nrows=None, ncol=None, i=None, pc_ratio=None):
        """Scatter plot (reversed order possible)"""
        ax = fig.add_subplot(nrows, ncol, i)
        if reverse:
            df_grouped = reversed(tuple(df.groupby(group_label)))
        else:
            df_grouped = df.groupby(group_label)
        for key, df_group in df_grouped:
            df_group.plot(ax=ax, kind='scatter', x='component 1', y='component 2',
                          edgecolor='none', alpha=0.75,
                          label=key, color=dict_colors[key], s=50)
        if pc_ratio is not None:
            xlabel = "PC1 ({}%)".format(int(pc_ratio[0] * 100))
            ylabel = "PC2 ({}%)".format(int(pc_ratio[1] * 100))
        else:
            xlabel = "PC1"
            ylabel = "PC2"
        plt.xlabel(xlabel, size=14)
        plt.ylabel(ylabel, size=14)
        if title is not None:
            plt.title(title, size=11)
        plt.tight_layout()

    # Venn Diagram
    @staticmethod
    def venn_diagramm(df=None, min_abs_log2_lfq=0.5, max_p=0.05, fig=None, nrows=None, ncol=None, i=None):
        """Venn diagramm for different days"""
        ax = fig.add_subplot(nrows, ncol, i)
        d3 = set(df[abs(df["log2 LFQ d03"]) >= min_abs_log2_lfq].index)
        d7 = set(df[abs(df["log2 LFQ d07"]) >= min_abs_log2_lfq].index)
        d14 = set(df[abs(df["log2 LFQ d14"]) >= min_abs_log2_lfq].index)
        venn3([d3, d7, d14], ["d03", "d07", "d14"], ax=ax)
        title = "Significant Genes (p_max: {}, abs_fc_min: {})".format(max_p, min_abs_log2_lfq)
        plt.style.use("ggplot")
        plt.title(title, size=12)
        plt.tight_layout()


class TSCSubPlots(TimeSeriesClustering, TSCPlottingData):
    """Class for evaluation and optimization of time series clustering"""
    def __init__(self, max_p=0.05, min_abs_log2_fc=0.5, list_max_p=None, list_min_abs_log2_fc=None,
                 sharma_max_p=0.05, sharma_min_log2_fc=1):
        TimeSeriesClustering.__init__(self)
        TSCPlottingData.__init__(self, sharma_max_p=sharma_max_p,
                                 sharma_min_log2_fold_change=sharma_min_log2_fc)
        self.list_th = self._get_th_list(max_p=max_p,
                                         min_abs_log2_lfq=min_abs_log2_fc,
                                         list_max_p=list_max_p,
                                         list_min_abs_log2_lfq=list_min_abs_log2_fc)

    # Helper methods
    @staticmethod
    def _check_legend_type(legend_type=None):
        """Check if legend type valid"""
        list_legend_types = ["cell_type", "cluster", "genes", None]
        if legend_type not in list_legend_types:
            raise ValueError("'legend_type' wrong. Select from {}".format(list_legend_types))

    # Settings
    @staticmethod
    def _get_th_list(max_p=0.05, min_abs_log2_lfq=0.5, list_max_p=None, list_min_abs_log2_lfq=None):
        """Get list of thresholds for fdr filtering
        In: a) """
        if list_max_p is not None:
            if list_min_abs_log2_lfq is not None:
                list_th = [[p, l] for p in list_max_p for l in list_min_abs_log2_lfq]
            else:
                list_th = [[p, min_abs_log2_lfq] for p in list_max_p]
        else:
            if list_min_abs_log2_lfq is not None:
                list_th = [[max_p, l] for l in list_min_abs_log2_lfq]
            else:
                list_th = [[max_p, min_abs_log2_lfq]]
        return list_th

    def _fig_n(self):
        """Get figure depending on selected thresholds"""
        n = len(self.list_th)
        dict_fig = {1: [5, 5],
                    2: [10, 5],
                    3: [10, 10],
                    4: [10, 10],
                    5: [10, 15],
                    6: [10, 15]}
        return dict_fig[n]

    def _subfig_n(self):
        """Get figure depending on selected thresholds"""
        n = len(self.list_th)
        dict_n = {1: {"nrows": 1, "ncol": 1},
                  2: {"nrows": 1, "ncol": 2},
                  3: {"nrows": 1, "ncol": 2},
                  4: {"nrows": 2, "ncol": 2},
                  5: {"nrows": 3, "ncol": 2},
                  6: {"nrows": 3, "ncol": 2}}
        return dict_n[n]

    # TODO refactor
    # Time series subplots
    def _get_colors(self, df_clusters=None, color_type=None, cluster=None):
        """Get colors for plot with all lines grouped by clusters or by cell types"""
        # Filter if specific cluster is selected
        if cluster is not None:
            df_clusters = df_clusters[df_clusters["cluster"] == cluster]
        # Select color based on legend_type
        if color_type == "cell_type":
            dict_color = self.dict_sharma_colors
            dict_color["nan"] = dict_color["Not Enriched"]
            colors = [dict_color[str(i)] for i in list(df_clusters["Identified Cell Type"])]
        elif color_type == "cluster":
            colors = [self.dict_colors_a[i] for i in list(df_clusters["cluster"])]
        else:
            colors = None
        return colors

    def _legend_colors(self, colors=None, legend_type=None, list_genes=None):
        """Get legend information"""
        if legend_type == "cell_type":
            custom_lines = [Line2D([0], [0], color=color, lw=2) for color
                            in self.sharma_legend_colors(list_ids=list_genes)]
            labels = self.sharma_labels(list_genes=list_genes)
        else:
            label_range = range(0, len(set(colors)))
            custom_lines = [Line2D([0], [0], color=self.colors_a[i], lw=2) for i in label_range]
            if legend_type == "cluster":
                labels = ["cluster {}".format(i) for i in label_range]
            else:
                labels = None
        return custom_lines, labels

    # TODO refactor
    def _time_series_add_subplot(self, df=None, title=None, ax=None, nrows=2, colors=None, legend_type=None, alpha=1.0):
        """Line plot for given df as subplot"""
        self._check_legend_type(legend_type=legend_type)
        y = 0.80
        x = self._dict_nrows_x[nrows]
        n = self._dict_nrows_legend[nrows]
        if len(df) > n or legend_type is None:
            legend = False
        else:
            legend = True
        font_arg = {"fontweight": "medium", "fontsize": self._dict_nrows_size[nrows]}
        df.transpose().plot(kind="line", ax=ax, legend=legend, color=colors, alpha=alpha)
        plt.ylim(self.min_val, self.max_val)
        plt.xticks(ticks=[0, 1, 2, 3], labels=["d00", "d03", "d07", "d14"])
        plt.ylabel("log2 ratio of mean LFQs", **font_arg)
        plt.xlabel("Time point", **font_arg)

        ax.text(y, x, title,
                horizontalalignment='center',
                transform=ax.transAxes)
        if legend_type in ["cell_type", "cluster"]:
            list_genes = list(self.add_col(df=df.copy())["gene_names"])
            custom_lines, labels = self._legend_colors(colors=colors,
                                                       legend_type=legend_type,
                                                       list_genes=list_genes)
            ax.legend(custom_lines, labels, loc="upper left")
        elif legend:
            plt.legend(loc="upper left")
            if n/2 < len(df) <= n:
                plt.legend(ncol=2, loc="upper left")

    @staticmethod
    def _add_cell_type_bars(df_clusters=None, ax=None, percentage=False, group_col=None,
                            label="cluster", colors=None):
        """Add overview with cell types per cluster"""
        sns.set_palette("deep")
        if percentage:
            df_clusters.replace(np.NAN, "Not Enriched", inplace=True)
        df = df_clusters.groupby([group_col])[label].value_counts().unstack(0).replace(np.NAN, 0)
        if percentage:
            df = df.apply(lambda x: x/x.sum()*100, axis=1)
            ylabel = "% Identified Cell Type"
            plt.ylim(0, 100)
        else:
            ylabel = "Count"
        df.plot(kind="bar", ax=ax, color=colors)
        ax.legend().set_title = ""
        plt.legend(ncol=2, loc="upper left", framealpha=0.5)
        plt.xlabel("Cluster")
        plt.ylabel(ylabel)
        sns.set_palette("bright")

    def time_series_subplots(self, df_clusters=None, n_clusters=4, title=None, group_col=None,
                             coloring=None, percentage=False):
        """Plot time series and clusters of time series with mean p value"""
        nrows, ncol = self.nrows_ncol(n_clusters=n_clusters)
        fig = plt.figure(figsize=(10, 10))
        plt.title(title, fontdict={"fontweight": "bold", "fontsize": 14})
        sns.despine(top=True, bottom=True, left=True, right=True)
        plt.axis('off')
        ax = fig.add_subplot(nrows, ncol, 1)
        title = "All genes"
        colors = self._get_colors(df_clusters=df_clusters, color_type=coloring)
        self._time_series_add_subplot(df=df_clusters[self.list_group_col],
                                      ax=ax,
                                      title=title,
                                      nrows=nrows,
                                      colors=colors,
                                      legend_type=coloring,
                                      alpha=0.7)
        for i in list(range(0, n_clusters)):
            df_clust = df_clusters[df_clusters["cluster"] == i][self.list_group_col]
            df_clust.index = [self.dict_ind_id[i] for i in df_clust.index]
            ax = fig.add_subplot(nrows, ncol, i + 2)
            title = "Cluster {} (n: {})".format(i, len(df_clust))
            colors = self._get_colors(df_clusters=df_clusters, color_type=coloring, cluster=i)
            if coloring in ["cell_type", "cluster"]:
                legend_type = None
            else:
                legend_type = coloring
            self._time_series_add_subplot(df=df_clust, title=title, ax=ax, nrows=nrows,
                                          colors=colors,
                                          legend_type=legend_type,
                                          alpha=0.7)
            plt.tight_layout()
        ax = fig.add_subplot(nrows, ncol, n_clusters + 2)
        list_ids = list(self.add_col(df=df_clusters.copy())["protein_id"])
        colors = self.sharma_legend_colors(list_ids=list_ids)
        self._add_cell_type_bars(df_clusters=df_clusters, ax=ax, group_col=group_col,
                                 colors=colors, percentage=percentage)
        plt.tight_layout()
        return plt.gca()



class TSCPlotting(TSCSubPlots, TSCSinglePlots):
    """Class for plotting time series data
    In: I Perseus arguments
        a1) df: df with raw values
        a2) conditions: conditions to select from df columns for processing
        a3) lfq_str: string to select lfq columns
        a4) group_rep: str part of group representative
        II Filtering thresholds
        b1) max_p: threshold of maximum
        b2) min_abs_log2_lfq: threshold of minimum absolute log2 lfq value
        b3) list_max_p: list of max_p values (for comparing multiple conditions)
        b4) list_min_abs_log2_lfq: list of min_abs_log2_lfq values (for comparing multiple conditions)
        III Filtering data
        c1) df_log2_lfq: df with log2 lfq values genes
        c2) df_p_vals: df with p values for genes
        c3) filter_conditions: list with columns to use for filtering"""
    def __init__(self, max_p=0.05, min_abs_log2_fc=0.5, list_max_p=None, list_min_abs_log2_fc=None,
                 df_log2_lfq=None, df_p_vals=None, filter_conditions=None, fig_format="png", n_max=12,
                 sharma_max_p=0.05, sharma_min_log2_fc=0.5):
        TSCSubPlots.__init__(self, max_p=max_p, min_abs_log2_fc=min_abs_log2_fc,
                             list_max_p=list_max_p, list_min_abs_log2_fc=list_min_abs_log2_fc,
                             sharma_max_p=sharma_max_p, sharma_min_log2_fc=sharma_min_log2_fc)
        TSCSinglePlots.__init__(self)
        # Filtering data
        self.n_max = n_max
        self.fig_format = fig_format
        self.df_log2_lfq = df_log2_lfq
        self.df_p_vals = df_p_vals
        self.filter_conditions = filter_conditions
        # Settings
        self._set_min_max(df=df_log2_lfq)

    # Plots
    def pca_clusters(self, df_ratio_p_val=None, model=None, model_arg=None, show=True, n_clusters=None,
                     metric="euclidean", david=False, title_out=None):
        """PCA plot for time series clusters"""
        # Clustering (with PCA plot)
        fig_n = self._fig_n()
        subfig_n = self._subfig_n()
        fig = plt.figure(figsize=fig_n)
        for i, p_lfq in enumerate(self.list_th):
            max_p = p_lfq[0]
            min_abs_log2_lfq = p_lfq[1]
            # Filtering
            df_fdr = self.get_df_fdr(df_ratio_pval=df_ratio_p_val,
                                     min_abs_log2_fc=min_abs_log2_lfq)
            # Clustering
            df_clusters = self.optimized_clustering(df_val=df_fdr,
                                                    model=model,
                                                    model_arg=model_arg,
                                                    n_min=4,
                                                    n_max=self.n_max,
                                                    n_clusters=n_clusters,
                                                    metric=metric)

            # David output
            if david:
                self.get_df_cluster_ids(df_clusters=df_clusters, max_p=max_p, min_abs_log2_lfq=min_abs_log2_lfq)
            # PCA
            title = "PCA for LFQ ratios (p_max: {}, abs_fc_min: {}, n: {})".format(max_p, min_abs_log2_lfq, len(df_fdr))
            labels = ["cluster {}".format(i) for i in list(set(df_clusters["cluster"]))]
            dict_colors = dict(zip(labels, self.colors_a))
            list_labels = ["cluster {}".format(i) for i in list(df_clusters["cluster"])]

            df_projected, pc_ratio = self.pca_on_distance(df=df_fdr, n_components=2, list_labels=list_labels)
            self.pca_plot(df=df_projected,
                          pc_ratio=pc_ratio,
                          group_label="labels",
                          dict_colors=dict_colors,
                          reverse=False,
                          title=title,
                          fig=fig, **subfig_n, i=i + 1)
        return plt.gca()

    def pca_cell_type(self, show=True, title_out=None):
        """PCA plot with different cell types"""
        # PCA cell type
        fig_n = self._fig_n()
        subfig_n = self._subfig_n()
        fig = plt.figure(figsize=fig_n)
        for i, p_lfq in enumerate(self.list_th):
            max_p = p_lfq[0]
            min_abs_log2_lfq = p_lfq[1]
            # Filtering
            df_fdr = self.get_df_fdr(df_log2_lfq=self.df_log2_lfq,
                                     df_p_vals=self.df_p_vals,
                                     filter_conditions=self.filter_conditions,
                                     max_p=max_p,
                                     min_abs_log2_fc=min_abs_log2_lfq)
            # PCA
            list_ids = [self.dict_ind_id[i] for i in df_fdr.index]
            list_labels = [self.dict_sharma[up_id] if up_id in self.dict_sharma.keys()
                           else "Not Enriched" for up_id in list_ids]
            title = "PCA for LFQ ratios (p_max: {}, abs_fc_min: {}, n: {})".format(max_p, min_abs_log2_lfq, len(df_fdr))
            df_projected, pc_ratio = self.pca_on_distance(df=df_fdr, n_components=2, list_labels=list_labels)
            self.pca_plot(df=df_projected,
                          pc_ratio=pc_ratio,
                          group_label="labels",
                          dict_colors=self.dict_sharma_colors,
                          reverse=True,
                          title=title,
                          fig=fig, **subfig_n, i=i+1)
        if title_out is None:
            title_out = "PCA_cell_types"
        return plt.gca()

    def venn_diagram_conditions(self, show=True, title_out=None):
        """Venn diagram for data sets"""
        fig_n = self._fig_n()
        subfig_n = self._subfig_n()
        fig = plt.figure(figsize=fig_n)
        for i, p_lfq in enumerate(self.list_th):
            max_p = p_lfq[0]
            min_abs_log2_lfq = p_lfq[1]
            df_fdr = self.get_df_fdr(df_log2_lfq=self.df_log2_lfq,
                                     df_p_vals=self.df_p_vals,
                                     filter_conditions=self.filter_conditions,
                                     max_p=max_p,
                                     min_abs_log2_fc=min_abs_log2_lfq)
            self.venn_diagramm(df=df_fdr, max_p=max_p, min_abs_log2_lfq=min_abs_log2_lfq,
                               fig=fig, **subfig_n, i=i+1)
        if title_out is None:
            title_out = "Venn_diagram"
        self._show_or_save(show=show, title_out=title_out)

    # TODO adjust
    def tsc_plots(self, model=None, model_arg=None, show=True, n_clusters=None, metric="euclidean", coloring="genes",
                  title_out=None, percentage=False):
        """Subplots for time series clustering
        In: a) coloring: 'cluster', 'cell_type', 'genes' """
        # Time Series Clustering
        for i, p_lfq in enumerate(self.list_th):
            max_p = p_lfq[0]
            min_abs_log2_lfq = p_lfq[1]
            df_fdr = self.get_df_fdr(df_log2_lfq=self.df_log2_lfq.copy(),
                                     df_p_vals=self.df_p_vals,
                                     filter_conditions=self.filter_conditions,
                                     max_p=max_p,
                                     min_abs_log2_fc=min_abs_log2_lfq)
            df_clusters = self.optimized_clustering(df_val=df_fdr,
                                                    model=model,
                                                    model_arg=model_arg,
                                                    n_min=4,
                                                    n_max=self.n_max,
                                                    n_clusters=n_clusters,
                                                    metric=metric)
            n_clusters = len(set(df_clusters["cluster"]))
            group_col = "Identified Cell Type"
            list_ids = [self.dict_ind_id[i] for i in df_fdr.index]
            df_clusters[group_col] = [self.dict_sharma[i] if i in list(self.dict_sharma.keys())
                                      else np.NAN for i in list_ids]
            sns.set_style("whitegrid")
            title = "Time series of LFQ ratios (p_max: {}, abs_fc_min: {}, n: {})".format(max_p,
                                                                                          min_abs_log2_lfq,
                                                                                          len(df_fdr))
            self.time_series_subplots(df_clusters=df_clusters,
                                      n_clusters=n_clusters,
                                      title=title,
                                      group_col=group_col,
                                      coloring=coloring,
                                      percentage=percentage)
            return plt.gca()


# III Test/Caller Functions


# IV Main



