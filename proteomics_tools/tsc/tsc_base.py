"""
This is a script for time series clustering
"""
import numpy as np


# I Helper Functions
def check_col_in_df(df=None, test_cols=None):
    """Check if test_col not None and in columns of given df"""
    if test_cols is None:
        raise ValueError("Column should not be None")
    list_cols = list(df)
    for col in test_cols:
        if col not in list_cols:
            raise ValueError("{} not in following df columns: {}".format(col, list_cols))


# II Main Functions
class TimeSeriesClusteringBase:
    """Class for times series clustering with following steps:
    1. Preprocessing (PerseusPipeline)
    2. Statistical tests (ttest and correlation)
    3. Pairwise distance ('cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan')
    5. Evaluation via calinski_harabasz_score, davies_bouldin_score, silhouette_score
    Q: https://scikit-learn.org/stable/modules/classes.html#module-sklearn.cluster"""

    def __init__(self):
        pass

    @staticmethod
    def get_df_fdr(df_ratio_pval=None, min_abs_log2_fc=0.5, max_pval=0.05, cols_ratio=None, cols_pval=None,
                   p_val_log10_in=True, drop_na=True):
        """Filter df_fc_pval by threshold for minimum absolute log2 fold change/ratio [fc] and maximum p value
        for given columns of fc and pval"""
        check_col_in_df(df=df_ratio_pval, test_cols=cols_ratio)
        check_col_in_df(df=df_ratio_pval, test_cols=cols_pval)
        df = df_ratio_pval.copy()
        # Filter for fold change/ratio
        df = df[[x.any() for x in np.array(abs(df[cols_ratio])) >= min_abs_log2_fc]]
        # Filter by p value
        if p_val_log10_in:
            min_pval = -np.log10(max_pval)
            df_fdr = df[[x.any() for x in np.array(df[cols_pval]) >= min_pval]]
        else:
            df_fdr = df[[x.any() for x in np.array(df[cols_pval]) <= max_pval]]
        if drop_na:
            df_fdr = df_fdr.dropna()
        df_fdr.reset_index(drop=True, inplace=True)
        return df_fdr

    @staticmethod
    def normalize_col(df=None, col_ref=None, cols=None, log_scale=True):
        """Get normalized for log2 lfq for first column in df"""
        df = df.copy()
        df_ref = df[col_ref]
        for col in cols:
            if log_scale:
                df[col] = df[col] - df_ref
            else:
                df[col] = df[col] / df_ref
        return df

