"""
This is a script for ...
"""
import time
import pandas as pd

# Settings
import examples._config as conf
from proteomics_tools.api_loaders import KEGGDB

pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def kegg_caller():
    """"""
    folder_out = conf.FOLDER_DATA + "GO_KEGG/"
    for organism in ["mmu", "hsa"]:
        kd = KEGGDB(organism=organism)
        # 1. Get KEGG pathway via Biopython REST API and save
        df_pathway = kd.get_pathway(verbose=True)
        """
        file_pathway = folder_out + "kegg_pathway_{}".format(organism)
        df_pathway.to_excel(file_pathway + ".xlsx", index=False)
        df_pathway.to_csv(file_pathway + ".tab", index=False, sep="\t")
        """
        # 1. Get KEGG disease via Biopython REST API and save
        df_disease = kd.get_disease(df_pathway=df_pathway)
        print(df_disease)
        """
        file_disease = folder_out + "kegg_disease_{}".format(organism)
        df_disease.to_excel(file_disease + ".xlsx", index=False)
        df_disease.to_csv(file_disease + ".tab", index=False, sep="\t")
        """

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    kegg_caller()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()


