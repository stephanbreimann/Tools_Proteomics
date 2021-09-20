"""
This is a script for ...
"""
import time
import pandas as pd
import examples._config as conf
import numpy as np


# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
class ParserSubstrates:
    """Class for parsing substrates from source to uniprot accession numbers"""
    def __int__(self):
        pass

    @staticmethod
    def _adjust_mapper(df):
        """Adjust of mapper df"""
        df["Protein names"] = df["Protein names"].str.upper()
        df["Gene names"] = df["Gene names"].str.upper()
        df.dropna(subset=["Gene_Name", "Protein names"], inplace=True)
        df.reset_index(drop=True, inplace=True)
        return df

    def one_to_one_acc(self, df, df_mapper=None, col_name="Gene names"):
        """One to one mapper"""
        df_mapper = self._adjust_mapper(df_mapper)
        list_acc = []
        for i, row in df.iterrows():
            gene = row[col_name]
            acc = row["ACC"]
            if str(acc) == "nan":
                col_ = df_mapper[df_mapper["Gene_Name"] == gene.upper()]
                # Exact matching
                if len(col_) < 1:
                    col_ = df_mapper[df_mapper["Protein names"] == gene.upper()]
                if len(col_) < 1:
                    col_ = df_mapper[df_mapper["ID"].str.split("_") == gene.upper()]
                # Substring matching
                if len(col_) < 1:
                    # Replace nan with False
                    mask = [True if x is True else False for x in df_mapper["Gene_Name"].str.contains(gene.upper())]
                    col_ = df_mapper[mask]
                if len(col_) < 1:
                    mask = [True if x is True else False for x in df_mapper["Protein names"].str.contains(gene.upper())]
                    col_ = df_mapper[mask]
                # Except if only one row matches
                if len(col_) == 1:
                    acc = col_["ACC"].values[0]
                else:
                    acc = np.NAN
            list_acc.append(acc)
        df["ACC"] = list_acc
        return df

    def one_to_may_acc(self, df, df_mapper=None, col_name="Gene names"):
        """One to many mapper. Results are printed so that selection can be performed manually on pre-filtered data"""
        df_mapper = self._adjust_mapper(df_mapper)
        list_no_match = []
        for i, row in df.iterrows():
            gene = row[col_name]
            acc = row["ACC"]
            if str(acc) == "nan":
                mask = [True if x is True else False for x in df_mapper["Gene_Name"].str.contains(gene.upper())]
                col_ = df_mapper[mask]
                if len(col_) < 1:
                    mask = [True if x is True else False for x in df_mapper["Protein names"].str.contains(gene.upper())]
                    col_ = df_mapper[mask]
                if len(col_) != 0:
                    print(gene)
                    print(col_)
                else:
                    list_no_match.append(gene)
        print("No matches for following {} genes".format(len(set(list_no_match))))
        for gene in set(list_no_match):
            print(gene)

    @staticmethod
    def add_mapper_info(df, df_mapper):
        """Merge df and adjust column naming"""
        df_filtered = df_mapper[df_mapper["ACC"].isin(set(df["ACC"]))][["ACC", "ID", "Gene_Name", "Organism"]]
        df_merged = df.merge(df_filtered, how="left", on="ACC")
        cols_x = [x for x in list(df_merged) if "_x" in x]
        cols_y = [y for y in list(df_merged) if "_y" in y]
        df_merged.drop(cols_x, axis=1, inplace=True)
        df_merged.rename({y: y.replace("_y", "") for y in cols_y}, axis=1, inplace=True)
        return df_merged

    @staticmethod
    def list_bace_substrates(file_hemming=None, min_hl_ratio_hemming=3, sum_min=1):
        """Get list of BACE substrates based on filtering thresholds"""
        df_bace = pd.read_excel(file_hemming, sheet_name="Summary")
        df = pd.read_excel(file_hemming, sheet_name="Hemming 2010")
        df_bace_hemming = df[df["Ratio (H/L)"] >= min_hl_ratio_hemming]
        list_genes_hemming = list(df_bace_hemming["Gene"])
        list_genes = df_bace["Gene"]
        df_bace["Hemming"] = [1 if x in list_genes_hemming else 0 for x in list_genes]
        df_bace.index = list_genes
        col_datasets = [x for x in list(df_bace) if x not in ["Gene", "Sum"]]
        df_bace["Sum"] = list(df_bace[col_datasets].transpose().sum())
        list_bace_sub = df_bace[df_bace["Sum"] >= sum_min].index.tolist()
        return list_bace_sub

    def df_bace_substrates(self, file_hemming=None, df_mapper=None):
        """"""
        list_bace1_sub = self.list_bace_substrates(file_hemming=file_hemming)
        df_bace = df_mapper[df_mapper["Gene_Name"].isin([x.upper() for x in list_bace1_sub])]
        df_bace = df_bace[["ACC", "ID", "Gene_Name", "Organism"]]
        df_bace.insert(0, "Hemming_2010", df_bace["Gene_Name"])
        list_missing_sub = list(set([x.upper() for x in list_bace1_sub]).difference(set(df_bace["Gene_Name"])))
        list_missing_sub.sort()
        dict_missing = {"Hemming_2010": list_missing_sub,
                        "ACC": [0] * len(list_missing_sub),
                        "ID": [0] * len(list_missing_sub),
                        "Gene_Name": [0] * len(list_missing_sub),
                        "Organism": [0] * len(list_missing_sub)}
        df_missing = pd.DataFrame(data=dict_missing)
        df_bace = pd.concat([df_bace, df_missing], axis=0)
        df_bace.reset_index(drop=True, inplace=True)
        return df_bace


# III Test/Caller Functions
def parse_adam_substrates():
    """"""
    folder_mapper = conf.FOLDER_DATA + "UniProt/"
    folder_adam = conf.FOLDER_DATA + "Substrates/ADAM/"
    file_adam = folder_adam + "ADAM_Substrates.xlsx"
    df_adam = pd.read_excel(file_adam)
    #df_adam.to_excel(file_adam, index=False)
    df_human = pd.read_excel(folder_mapper + "mapping_db_hsa.xlsx")
    #ps = ParserSubstrates()
    #df = ps.one_to_one_acc(df_adam, df_mapper=df_human, col_name="Lambrecht_2018")
    #ps.one_to_may_acc(df_adam, df_mapper=df_human, col_name="Lambrecht_2018")
    #df_merged = ps.add_info(df_adam, df_mapper=df_human)
    #df_merged.to_excel(folder_adam + "test.xlsx")


def adjust_gamma_sec():
    """"""
    folder_in = conf.FOLDER_DATA + "Substrates/"
    """
    df = pd.read_excel(folder_in + "substrates_GAMMA-SEC.xlsx")
    print(df)
    df_mapper = pd.read_csv(conf.folder_data + "UniProt/" + "uniprot_hsa.tab", sep="\t")
    list_genes = df["Gene_Name"].tolist()
    
    df_filtered = df_mapper[df_mapper["Gene_Name"].isin(list_genes)]
    dict_gene_acc = dict(zip(df_filtered["Gene_Name"], df_filtered["ACC"]))
    list_acc_hsa = []
    for i, row in df.iterrows():
        gene = row["Gene_Name"]
        if gene in dict_gene_acc.keys():
            list_acc_hsa.append(dict_gene_acc[gene])
        else:
            list_acc_hsa.append(np.NAN)
    df["ACC_HUMAN"] = list_acc_hsa
    """
    df = pd.read_excel(folder_in + "non-substrates_GAMMA-SEC.xlsx")
    print(df)
    df_sub = pd.read_excel(folder_in + "substrates_GAMMA-SEC.xlsx")
    print(df_sub[df_sub["class"] == "SUBEXPERT"])
    """
    loc_id = df.columns.get_loc("ID")
    df.insert(loc_id + 1, "Organism", [dict_org[x.split("_")[1]]["UP"] for x in df["ID"]])
    df.insert(loc_id + 1, "Gene_Name", [x.split("_")[0] for x in df["ID"]])
    df.to_excel(folder_in + "Substrates_GAMMA-SEC.xlsx", index=False)
    """


def adjust_substrates():
    """"""
    file = "non-substrates_PARL.xlsx"
    folder_in = conf.FOLDER_DATA + "Substrates/PARL/"
    folder_up = conf.FOLDER_DATA + "UniProt/"

    df = pd.read_excel(folder_in + file)
    df_mapper = pd.read_csv(folder_up + "uniprot_hsa.tab", sep="\t")
    ps = ParserSubstrates()
    ps.one_to_one_acc()
    df = ps.add_mapper_info(df, df_mapper)
    df.to_excel(conf.FOLDER_DATA + "Substrates/" + file, index=False)
    print(df)
    """
    loc_id = df.columns.get_loc("ACC")
    df.insert(loc_id + 1, "Organism", [dict_org[x.split("_")[1]]["UP"] for x in df["ID"]])
    df.insert(loc_id + 1, "Gene_Name", [x.split("_")[0] for x in df["ID"]])
    #ps = ParserSubstrates()
    """


def sppl_substrates():
    """"""
    file = "substrates_SPPL3.xlsx"
    folder_in = conf.FOLDER_DATA + "Substrates/"
    folder_up = conf.FOLDER_DATA + "UniProt/"
    df = pd.read_excel(folder_in + file)
    df["ID"] = [x.split("(")[1].replace(")", "").strip() for x in df["ACC"]]
    df["ACC"] = [x.strip().split(" ")[2] for x in df["ACC"]]
    loc_id = df.columns.get_loc("ACC")
    df.insert(loc_id + 1, "Gene_Name", [x.split("_")[0] for x in df["ID"]])
    df.to_excel(folder_in + file, index=False)
    print(df)
    """
    #df = TmdParser().read_df(df=df)
    loc_id = df.columns.get_loc("ACC")
    df.insert(loc_id + 1, "Organism", [dict_org[x.split("_")[1]]["UP"] for x in df["ID"]])
    df.insert(loc_id + 1, "Gene_Name", [x.split("_")[0] for x in df["ID"]])
    df.to_excel(folder_in + file, index=False)
    print(df)
    """


# IV Main
def main():
    t0 = time.time()
    #parse_adam_substrates()
    #adjust_gamma_sec()
    #adjust_substrates()
    sppl_substrates()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
