"""
This is a script for ...
"""
import time
import pandas as pd
import matplotlib.pyplot as plt

import examples._config as conf
from proteomics_tools.gene_enrichment import GeneEnrichment

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions
def go(pop=None, study=None, organism="human", folder_out=None):
    """"""
    p_method = ["fdr", "bonferroni", "holm"]
    ge = GeneEnrichment(organism=organism, verbose=True)
    dict_go, dict_org_go, df_org_go = ge.load_go()
    dict_kegg = ge.load_kegg()
    # GO enrichment analysis
    study, pop, df_enrich = ge.go_enrichment(study=study,
                                             pop=pop,
                                             df_org_go=df_org_go,
                                             dict_go=dict_go,
                                             p_method=p_method)
    print(len(study))
    df_enrich_filtered = ge.filter_df(df=df_enrich,
                                      p_max=0.05,
                                      p_method="p_fdr")
    ge.plot_enrichment_go(df_enrich_filtered=df_enrich_filtered,
                          p_col="p_fdr",
                          n=5,
                          split_name=True)
    #plt.show()
    # 2 Semantic clustering (via REVIGO)
    """
    df_sem_clust = ge.semantic_clustering(df_enrich=df_enrich,
                                          df_enrich_filtered=df_enrich_filtered,
                                          dict_go=dict_go,
                                          dict_org_go=dict_org_go,
                                          n=5)
    #print(df_sem_clust)
    #ge.plot_revigo_clustering(df_sem_clust=df_sem_clust)
    #plt.show()
    plt.tight_layout()

    plt.savefig(folder_out + "plot_revigo_clustering.png")
    """

    # 2.3 KEGG enrichment analysis
    df_pathway, df_disease, df_study = ge.kegg_analysis(df_enrich_filtered=df_enrich_filtered,
                                                        study=study,
                                                        pop=pop,
                                                        p_method=p_method,
                                                        dict_kegg=dict_kegg,)
    ge.plot_enrichment_kegg(df_pathway=df_pathway, df_disease=df_disease)
    plt.show()
    df_enrich_filtered.to_excel(folder_out + "df_enrich.xlsx", index=False)
    #df_sem_clust.to_excel(folder_out + "df_revigo.xlsx", index=False)
    df_pathway.to_excel(folder_out + "df_pathway.xlsx", index=False)
    df_disease.to_excel(folder_out + "df_disease.xlsx", index=False)
    df_study.to_excel(folder_out + "df_study.xlsx", index=False)

# II Main Functions


# III Test/Caller Functions
def go_analysis():
    """"""
    folder_in = conf.FOLDER_DATA + "datasets" + conf.SEP
    f = conf.FOLDER_RESULTS + "substrates" + conf.SEP
    df = pd.read_excel(folder_in + "df_tmd_pred_ss_type_I_hsa_mmu_SB_21_06_14.xlsx")
    print(df)
    df.insert(2, "organism", [x.split("_")[1] for x in df["name"]])
    i = list(df).index("len_c_terminal")
    df.insert(i+1, "len_seq", [len(s) for s in df["sequence"]])
    df.insert(i, "len_tmd", [len(x) for x in df["tmd"]])
    df.insert(i, "len_n_terminal", [int(x)-1 for x in df["start"]])
    df_mmu = pd.read_excel(f + "df_study_mmu_candidates.xlsx")
    df_hsa = pd.read_excel(f + "df_study_hsa_candidates.xlsx")


    """
    
    df_go = pd.concat([df_hsa, df_mmu], axis=0)
    df_go = df_go.set_index("STUDY")
    df_go.index.name = "entry"
    df_go = df_go.drop("NAME", axis=1)
    list_ids = list(df_mmu["STUDY"]) + list(df_hsa["STUDY"])
    df_candidates = df[df["entry"].isin(list_ids)]
    df_candidates = df_candidates.set_index("entry")
    df_candidates = df_candidates.join(df_go)
    df = df_candidates.sort_values(by="len_n_terminal")
    df.to_excel(f + "df_candidates_21_06_23_SB.xlsx")
    """
    """
    df["NAME"] = [x.split("_")[0] for x in df["name"]]
    
    list_name = list(df[df["class"].isin(["SUBEXPERT", "SUBLIT"])]["NAME"])
    df = df[~df["NAME"].isin(list_name)]
    df.drop("NAME", inplace=True, axis=1)
    df.to_excel(f + "df_all.xlsx", index=False)
    
    df_mmu = df[df["organism"] == "MOUSE"]
    df_hsa = df[df["organism"] == "HUMAN"]
    df_mmu.to_excel(f + "df_mmu.xlsx", index=False)
    df_hsa.to_excel(f + "df_hsa.xlsx", index=False)
    folder_out = f + "mmu" + conf.sep
    
    go(pop=df_mmu["entry"], study=df_mmu[df_mmu["MEAN"] > 60]["entry"], organism="mouse", folder_out=folder_out)
    folder_out = f + "hsa" + conf.sep
    go(pop=df_hsa["entry"], study=df_hsa[df_hsa["MEAN"] > 60]["entry"], organism="human", folder_out=folder_out)
    """


# IV Main
def main():
    t0 = time.time()
    go_analysis()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
