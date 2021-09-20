"""
This is a script for clustering optimization and time series clustering
"""

import numpy as np
import pandas as pd
from sklearn.metrics.cluster import calinski_harabasz_score, silhouette_score, davies_bouldin_score
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

from proteomics_tools.tsc.tsc_base import TimeSeriesClusteringBase


# I Helper Functions
def check_col_in_df(df=None, test_cols=None):
    """Check if test_col not None and in columns of given df"""
    if test_cols is None:
        raise ValueError("Column should not be None")
    list_cols = list(df)
    for col in test_cols:
        if col not in list_cols:
            raise ValueError("{} not in following df columns: {}".format(col, list_cols))


class ClusteringOptimization:
    """Optimization of clustering using evaluation metrics"""
    def __init__(self):
        pass

    # Helper methods
    @staticmethod
    def _check_score(score=None):
        """Check if score valid"""
        list_scoring = ["SS", "CHS", "DBS"]     # Silhouette score, Calinksi Harabasz score, Davies-Bouldin Index
        if score not in list_scoring:
            raise ValueError("'scoring' must be one of following: ", list_scoring)

    @staticmethod
    def _evaluation_metrics(df_val=None, labels=None):
        """Evaluation of clustering algorithms using silhouette_score, calinski_harabasz_score, davies_bouldin_score
        Q:  a) https://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation"""
        # Silhouette score: from -1 (incorrect clustering) over 0 (overlapping clustering) to 1 (highly dense)
        sil_score = silhouette_score(df_val, labels=labels)
        sil_score = (sil_score + 1) / 2  # 0 incorrect, 1 highly dense
        # Calinksi Harabasz score (Variance Ratio Criterion): sum of between-clusters d/ inter-cluster d (all clusters)
        # where d dispersion is defined as the sum of distances squared
        # The higher the score, the more dense are the clusters (ranges from 0 to max 1000 for data set)
        cH_score = calinski_harabasz_score(df_val, labels=labels)
        cH_score = cH_score / 1000  # (1000 - c) / 1000
        # Davies-Bouldin Index (average ‘similarity’ between clusters)
        # The smaller, the lower the similarity between clusters, the higher the partition
        db_score = davies_bouldin_score(df_val, labels=labels)
        db_score = 1 - db_score
        scores = [round(val, 2) for val in [sil_score, cH_score, db_score]]
        return scores

    # Main methods
    @staticmethod
    def clustering_on_distance(df_val=None, n_clusters=7, metric="euclidean", model=None, model_arg=None):
        """Perform clustering on pairwise distances
        In: a) metric: metric for pairwise distances ['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan']
        Q:  a) https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html"""
        df_val = df_val.copy()
        p = pairwise_distances(df_val, metric=metric)
        clust_model = model(**model_arg, n_clusters=n_clusters)
        labels = clust_model.fit_predict(p)
        return labels

    # Clustering optimization
    def clustering_evaluation(self, df_val, model=None, model_arg=None, rank=True,
                              n_min=4, n_max=12, metric="euclidean"):
        """Get df with number of clusters and evaluation scores"""
        dict_n_scores = {}
        for i in list(range(n_min, n_max)):
            labels = self.clustering_on_distance(df_val=df_val,
                                                 n_clusters=i,
                                                 model=model,
                                                 model_arg=model_arg,
                                                 metric=metric)
            list_evaluation = self._evaluation_metrics(df_val=df_val, labels=labels)
            unique, counts = np.unique(labels, return_counts=True)
            list_evaluation.extend([min(counts), max(counts)])
            dict_n_scores[i] = list_evaluation
        df_eval = pd.DataFrame.from_dict(dict_n_scores, orient="index", columns=["SS", "CHS", "DBS", "min", "max"])
        if rank:
            df_eval = df_eval.rank(axis=0, ascending=False)
        return df_eval

    def select_optimum_n_clusters(self, df_eval=None, score=None, rank=True):
        """Select best number of clusters based on ranking of clusters"""
        if score is None:
            mean_rank = df_eval.mean(axis=1)
            n_clusters = mean_rank.idxmin()
        else:
            self._check_score(score=score)
            if rank:
                df_top = pd.concat([df_eval.idxmin(), df_eval.min()], axis=1)
            else:
                df_top = pd.concat([df_eval.idxmax(), df_eval.max()], axis=1)
            df_top.columns = ["index", "score"]
            n_clusters = int(df_top.loc[score]["index"])
        return n_clusters


# II Main Functions
class TimeSeriesClustering(TimeSeriesClusteringBase):
    """Class for times series clustering with following steps:
    1. Preprocessing (PerseusPipeline)
    2. Statistical tests (ttest and correlation)
    3. Pairwise distance ('cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan')
    4. Evaluation via calinski_harabasz_score, davies_bouldin_score, silhouette_score
    Q: https://scikit-learn.org/stable/modules/classes.html#module-sklearn.cluster"""
    def __init__(self):
        TimeSeriesClusteringBase.__init__(self)

    # Helper function
    @staticmethod
    def get_max_cluster(df_clusters=None, col_label="cluster"):
        """Get cluster with most data points from df_clust"""
        max_cluster = df_clusters[col_label].value_counts().index.tolist()[0]
        return max_cluster

    # Main methods
    @staticmethod
    def pca_on_distance(df_val=None, metric="euclidean", n_components=2):
        """Perform PCA on pairwise distance matrix
        Q:  a) https://stats.stackexchange.com/questions/87681/performing-pca-with-only-a-distance-matrix
        !:  a) Just on euclidean distance"""
        p = pairwise_distances(df_val, metric=metric)
        pca = PCA(n_components=n_components)
        projected = pca.fit_transform(p)
        pc_ratio = pca.explained_variance_ratio_
        df_projected = pd.DataFrame(projected, columns=["component {}".format(i) for i in range(1, n_components + 1)])
        return df_projected, pc_ratio

    @staticmethod
    def optimized_clustering(df_val=None, model=None, model_arg=None, n_min=4, n_max=12,
                             score="SS", metric="euclidean", n_clusters=None):
        """Clustering with optimization of number of clusters
        In: a) df_val: df with values to perform clustering on
            b1) model: clustering model
            b2) model_arg: arguments of clustering model,
            b3) n_min, n_max: range of minimum to maximum number of clusters
            b4) score: clustering optimization score {'SS': Silhouette, 'CHS': Calinksi Harabasz, 'DBS': Davies-Bouldin}
            b5) metric: distance metric for df_val
            c) n_clusters: number of clusters -> no optimization
        Out:a) labels: list of clustering labels for each row of df_val
        """
        co = ClusteringOptimization()
        if n_clusters is None:
            df_eval = co.clustering_evaluation(df_val=df_val,
                                               model=model,
                                               model_arg=model_arg,
                                               n_min=n_min,
                                               n_max=n_max,
                                               metric=metric)
            n_clusters = co.select_optimum_n_clusters(df_eval=df_eval,
                                                      score=score)
        labels = co.clustering_on_distance(df_val=df_val,
                                           model=model,
                                           model_arg=model_arg,
                                           n_clusters=n_clusters,
                                           metric=metric)
        return labels

    @staticmethod
    def get_df_cluster_ids(df_clusters=None, col_id="ACC", col_label="cluster"):
        """Create gene list for each clustering (for David analysis)"""
        dict_clust_list = {"cluster {}".format(i): [] for i in set(list(df_clusters[col_label]))}
        for i, row in df_clusters.iterrows():
            entry = row[col_id]
            clust = row[col_label]
            dict_clust_list["cluster {}".format(clust)].append(entry)
        len_max = max([len(dict_clust_list[i]) for i in dict_clust_list])
        for i in dict_clust_list:
            ids = dict_clust_list[i]
            list_empty = [""] * (len_max - len(ids))
            dict_clust_list[i].extend(list_empty)
        df_clusters_ids = pd.DataFrame(dict_clust_list)
        return df_clusters_ids

    @staticmethod
    def get_df_venn(df_fdr=None, cols_group=None, labels_group=None, col_id="ACC", min_abs_log2_lfq=0.5):
        """Get protein ids as given in venn diagram"""
        # Check input
        if cols_group is None or labels_group is None or len(cols_group) != len(labels_group):
            raise ValueError("Size of 'cols_group' and 'labels_group' must match")
        check_col_in_df(df=df_fdr, test_cols=cols_group)
        # Create dict with information for venn diagram
        dict_group_ids = {}
        for col, label in zip(cols_group, labels_group):
            list_ids = list(df_fdr[abs(df_fdr[col]) >= min_abs_log2_lfq][col_id])
            dict_group_ids[label] = list_ids
        dict_venn = {x: "" for x in list(df_fdr[col_id])}
        for group in dict_group_ids:
            for id_ in dict_group_ids[group]:
                dict_venn[id_] += group
        list_venn = [dict_venn[p] for p in list(df_fdr[col_id])]
        return list_venn


class TSCSaver:
    """Save results of time series clustering"""
    def __init__(self):
        pass

    @staticmethod
    def save_cluster(df_clusters=None, max_p=0.05, min_abs_log2_fc=0.5):
        """Save df cluster"""
        title = "df_clusters_p_max_{}_min_abs_log2_fc_{}".format(max_p, min_abs_log2_fc)
        df_clusters = df_clusters.sort_value(by="cluster")
        return title, df_clusters

    @staticmethod
    def save_gene_list(df_clusters_ids=None, max_p=0.05, min_abs_log2_fc=0.5):
        """Save df cluster ids containing lists of ids for each cluster"""
        title = "df_gene_list_p_max_{}_min_abs_log2_fc_{}".format(max_p, min_abs_log2_fc)
        return title, df_clusters_ids

