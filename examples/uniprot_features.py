"""
This is a script for ...
"""
import time
import pandas as pd

from proteomics_tools.uniprot_features import UniprotFeatures
import examples._config as conf

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


def mapping_org():
    """"""
    folder_mapper = conf.FOLDER_DATA + "UniProt/"
    df_human = pd.read_excel(folder_mapper + "mapping_db_hsa_reviewed.xlsx")
    df_mouse = pd.read_excel(folder_mapper + "mapping_db_mmu_reviewed.xlsx")
    df_mapper_org = UniprotFeatures().mapper_org(df_org1=df_human, df_org2=df_mouse,
                                                 str_org1="hsa", str_org2="mma")
    file_mapper = folder_mapper + "mapping_hsa_to_mmu.xlsx"
    df_mapper_org.to_excel(file_mapper, index=False)


def up_features():
    """"""
    folder_uniprot = conf.FOLDER_DATA + "UniProt/"
    list_name = ["uniprot_{}".format(w) for w in ["hsa", "hsa_reviewed", "mmu", "mmu_reviewed"]]
    list_name = ["uniprot_{}".format(w) for w in ["rat", "rat_reviewed"]]
    uf = UniprotFeatures()
    # Raw files in archieve
    for name in list_name:
        print(name)
        df = pd.read_csv(folder_uniprot + "raw_{}.tab".format(name), sep="\t")
        print(df.head(10))
        df = uf.run_features(df.copy())
        print(df.head(10))
        df.to_excel(folder_uniprot + "{}.xlsx".format(name), index=False)
        df.to_csv(folder_uniprot + "{}.tab".format(name), index=False, sep="\t")
        df = uf.run_one_hot_encoding(df)
        df.to_excel(folder_uniprot + "encoded_{}.xlsx".format(name), index=False)
        df.to_csv(folder_uniprot + "encoded_{}.tab".format(name), index=False, sep="\t")


def dummy_encode_data():
    """"""
    folder_uniprot = conf.FOLDER_DATA + "UniProt/"
    name = "uniprot_mmu_reviewed"
    df = pd.read_csv(folder_uniprot + "{}.tab".format(name), sep="\t")
    uf = UniprotFeatures()
    df = uf.run_one_hot_encoding(df)
    df.to_excel(folder_uniprot + "encoded_{}.xlsx".format(name), index=False)
    df.to_csv(folder_uniprot + "encoded_{}.text".format(name), index=False, sep="\t")


# IV Main
def main():
    t0 = time.time()
    #mapping_org()
    up_features()
    #dummy_encode_data()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
