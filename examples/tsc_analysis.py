"""
This is a script for ...
"""
import time
import pandas as pd
import matplotlib.pyplot as plt

from proteomics_tools.tsc import TimeSeriesClustering, TSCPlotting

from sklearn.cluster import AgglomerativeClustering
from proteomics_tools.perseuspy import PerseusPipeline, get_dict_groups
import examples._config as conf


# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def tsc_test():
    """"""
    df = pd.read_excel(conf.FOLDER_DATA + "Demyelination_lfq.xlsx")
    # 1. Perseus Analysis
    groups = ["d00", "d03", "d07", "d14"]
    dict_col_group = get_dict_groups(df=df, groups=groups)
    pp = PerseusPipeline(df=df, dict_col_group=dict_col_group)
    df_ratio_pval = pp.run()
    # 1.2 Filter columns
    list_cols = [x for x in list(df_ratio_pval) if "/d00" in x]
    list_ratio_col = [x for x in list_cols if "ratio" in x]
    list_pval_col = [x for x in list_cols if "p value" in x]

    # 2. TSC
    tsc = TimeSeriesClustering()
    df_fdr = tsc.get_df_fdr(df_ratio_pval=df_ratio_pval, cols_ratio=list_ratio_col, cols_pval=list_pval_col)
    print(df_fdr)
    model = AgglomerativeClustering
    model_arg = {"linkage": "average", "affinity": "precomputed"}
    df_clusters = df_fdr.copy()
    # 2.1 Clustering
    labels = tsc.optimized_clustering(df_val=df_fdr[list_ratio_col],
                                      model=model,
                                      model_arg=model_arg)
    df_clusters["cluster"] = labels
    max_cluster = tsc.get_max_cluster(df_clusters=df_clusters)
    # 2.2 Set analysis (Venn diagram)
    list_venn = tsc.get_df_venn(df_fdr=df_fdr, cols_group=list_ratio_col, labels_group=["d00", "d03", "d07", "d14"])
    df_clusters["venn"] = list_venn

    # 2.3 Gene lists
    df_cluster_ids = tsc.get_df_cluster_ids(df_clusters=df_clusters)
    df_pop = df.dropna(how="all", subset=[x for x in list(df) if x not in ["Protein ID", "Gene Names", "Description"]])
    tsc_plot = TSCPlotting()
    # Plotting
    plt.show()


# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    tsc_test()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()