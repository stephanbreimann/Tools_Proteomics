"""
This is a script for ...
"""
import time
import pandas as pd
import numpy as np

import examples._config as conf

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe
#warnings.simplefilter(action="ignore", category=UserWarning)

# I Helper Functions


# II Main Functions
def perseus_analysis():
    """"""
    # Load data
    file = conf.FOLDER_DATASETS + "Demyelination_lfq.xlsx"
    df = pd.read_excel(file)
    df_sub = pd.read_excel(conf.FOLDER_DATASETS + "df_tmd_pred_ss_type_I_hsa_mmu_SB_21_06_14.xlsx")
    vals = []
    for i in set(df_sub["class"]):
        d = df_sub
        d = d[d["class"] == i]
        start_mean = np.mean(d["start"])
        stop_mean = np.mean(d["stop"])
        vals.append([start_mean, stop_mean])
    df = pd.DataFrame(data=vals, columns=["start", "stop"])
    df["test"] = 1
    print(df)


    """
    # 1. Perseus Analysis
    groups = ["d00", "d03", "d07", "d14"]
    dict_col_group = get_dict_groups(df=df, groups=groups)
    pp = PerseusPipeline(df=df, dict_col_group=dict_col_group)
    df_ratio_pval = pp.run(gmean=False)
    print(df_ratio_pval)
    # 1.3 Plotting
    fig_format = "svg"
    gene_list = ["Sez6"]
    pp.volcano_plot(df_ratio_pval=df_ratio_pval,
                    gene_list=gene_list,
                    col_pval="-log10 p value (d07/d00)",
                    col_ratio="log2 ratio (d07/d00)",
                    th_filter=(0.05, 0.5),
                    th_text=(0.05, -1, 2.5),
                    force=(0.9, 0.50, 0.25),
                    avoid_conflict=0.3,
                    precision=0.01,
                    box=False,
                    verbose=True,
                    label_bold=False,
                    label_size=10,
                    filled_circle=False,
                    title="Test_Volcano_mean",
                    fig_format=fig_format)
    plt.show()
    # 2. Gene Enrichment
    organism = "mouse"
    p_method = ["fdr", "bonferroni", "holm"]
    ge = GeneEnrichment(organism=organism)
    dict_go, dict_org_go, df_org_go = ge.load_go()
    dict_kegg = ge.load_kegg()
    pop = df_ratio_pval["ACC"].tolist()
    study = df_ratio_pval[(df_ratio_pval["-log10 p value (d07/d00)"] >= -np.log(0.05)) &
                          (abs(df_ratio_pval["log2 ratio (d07/d00)"]) <= 0.5)]["ACC"].tolist()
    # 2.1 GO enrichment analysis
    study, pop, df_enrich = ge.go_enrichment(study=study,
                                             pop=pop,
                                             df_org_go=df_org_go,
                                             dict_go=dict_go,
                                             p_method=p_method)
    df_enrich_filtered = ge.filter_df(df=df_enrich,
                                      p_max=0.05,
                                      p_method="p_fdr")
    ge.plot_enrichment_go(df_enrich_filtered=df_enrich_filtered,
                          p_col="p_fdr",
                          n=5,
                          split_name=True)
    plt.show()
    # 2.2 Semantic clustering (via REVIGO)
    df_sem_clust = ge.semantic_clustering(df_enrich=df_enrich,
                                          df_enrich_filtered=df_enrich_filtered,
                                          dict_go=dict_go,
                                          dict_org_go=dict_org_go,
                                          n=5)
    ge.plot_revigo_clustering(df_sem_clust=df_sem_clust)
    plt.show()
    # 2.3 KEGG enrichment analysis
    df_pathway, df_disease, df_study = ge.kegg_analysis(df_enrich_filtered=df_enrich_filtered,
                                                        study=study,
                                                        pop=pop,
                                                        p_method=p_method)
    ge.plot_enrichment_kegg(df_pathway=df_pathway, df_disease=df_disease)
    plt.show()
    """

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    perseus_analysis()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
