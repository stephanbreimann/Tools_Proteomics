"""
This is a script for adding uniprot features (e.g., glycosylation or subcellular location) to data frame
"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder


# I Helper Functions


# II Main Functions
class UniprotFeatures:
    """Class for organizing features from uniprot"""
    def __init__(self):
        self.dict_cc = {"Nucleus": "Nucleus",
                        "Endoplasmic reticulum": "ER",
                        "Cytoplasm": "Cytoplasm",
                        "Cytoskeleton": "Cytoskeleton",
                        "Centrosome": "Centrosome",
                        "Mitochondrion": "Mitochondrion",
                        "Golgi": "Golgi",
                        "Endosome": "Endosome",
                        "Lysosome": "Lysosome",
                        "Cell membrane": "PM",  # Plasma Membrane
                        "Peripheral membrane protein": "PMPs",
                        "Membrane": "Membrane",
                        "Extracellular matrix": "ECM",
                        "Cell junction": "Cell junction",
                        "Cell projection": "Cell projection",
                        "Secreted": "Secreted"}
        self.dict_tm1_type = {"Type I ": "1",
                              "Type II ": "2",
                              "Type III ": "3",
                              "Type IV ": "4"}
        self.dict_lip = {"palmito": "Palmitoylation",
                         "geranylgeranyl": "Geranylgeranylation",
                         "myristoyl": "Myristoylation",
                         "GPI-anchor": "GPI-anchor",
                         "farnesyl": "Farnesylation"}
        self.dict_gly = {"N-linked": "N-linked",
                         "O-linked": "O-linked"}
        self.dict_ptm = {"Phospho": "Phosphorylation",
                         "Ubiqui": "Ubiquitination"}
        self.dict_col = {"Subcellular location [CC]": list(self.dict_cc.values()),
                         "Glycosylation": list(self.dict_gly.values()),
                         "Lipidation": list(self.dict_lip.values()),
                         "PTM": list(self.dict_ptm.values())}

    @staticmethod
    def mapper_org(df_org1=None, df_org2=None, str_org1="hsa", str_org2="mma"):
        """Match ACCs of tow organism based on IDs"""
        list_org1 = [x.split("_")[0] for x in df_org1["ID"]]
        list_org2 = [x.split("_")[0] for x in df_org2["ID"]]
        df_org1["Gene_Name"] = [x.split("_")[0] for x in df_org1["ID"]]
        df_org2["Gene_Name"] = [x.split("_")[0] for x in df_org2["ID"]]
        df_1 = df_org2[df_org2["Gene_Name"].isin(list_org1)].sort_values(by="ID")
        df_2 = df_org1[df_org1["Gene_Name"].isin(list_org2)].sort_values(by="ID")
        col_id_org1 = "ID_{}".format(str_org1)
        col_acc_org1 = "ACC_{}".format(str_org1)
        col_id_org2 = "ID_{}".format(str_org2)
        col_acc_org2 = "ACC_{}".format(str_org2)
        dict_mapper_org = {"Gene_Name": df_2["Gene_Name"].tolist(),
                           col_id_org1: df_2["ID"].tolist(),
                           col_acc_org1: df_2["ACC"].tolist(),
                           col_id_org2: df_1["ID"].tolist(),
                           col_acc_org2: df_1["ACC"].tolist()}
        df_mapper_org = pd.DataFrame(data=dict_mapper_org)
        return df_mapper_org

    # Parse uniprot features
    @staticmethod
    def _parse_column(df, col_name="Subcellular location [CC]", dict_items=None, out_list=False):
        """Extract information given in dict_items from column"""
        if dict_items is None:
            raise ValueError("'dict_items' should be given by attributes (e.g. self.dict_cc)")
        list_values = []
        for i, row in df.iterrows():
            val = row[col_name]
            if type(val) is str:
                list_val = []
                for c in val.replace(".", ";").split(";"):
                    for w in dict_items.keys():
                        if w.lower() in c.lower():
                            list_val.append(dict_items[w])
                list_val = list(set(list_val))
                list_val.sort()
                str_val = ";".join(list_val)
                if len(str_val) == 0:
                    str_val = np.NaN
            else:
                str_val = np.NaN
            list_values.append(str_val)
        if out_list:
            return list_values
        else:
            df[col_name] = list_values
            return df

    @staticmethod # TODO remove if proven to be wrong
    def _get_tm_type(df):
        """Add type of single span transmembrane protein"""
        dict_type = {"i": 1,
                     "ii": 2,
                     "iii": 3,
                     "iv": 4}
        list_type = []
        for i, row in df.iterrows():
            tm_type = np.NaN
            val = row["Transmembrane"]
            if type(val) is str:
                count_val = 0
                for c in val.replace(".", ";").split(";"):
                    if "TRANSMEM" in c:
                        count_val += 1
                if count_val == 1:
                    for n in ["ii", "iii", "iv"]:
                        str_type = "type {} ".format(n)
                        if str_type in val.lower():
                            tm_type = dict_type[n]
                    if tm_type is np.NaN:
                        tm_type = 1     # If not specified assumed to be type I (minority of TPs)
            list_type.append(tm_type)
        return list_type

    @staticmethod
    def _count_str(df, col_name="Transmembrane", str_count="TRANSMEM"):
        """Count TM"""
        list_values = []
        for i, row in df.iterrows():
            val = row[col_name]
            if type(val) is str:
                count_val = 0
                for c in val.replace(".", ";").split(";"):
                    if str_count in c:
                        count_val += 1
            else:
                count_val = 0
            list_values.append(count_val)
        return list_values

    @staticmethod
    def _adjust_col_names(df, dict_name=None):
        """Adjust column names to standardized form"""
        df.rename(dict_name, inplace=True, axis=1)
        df["Gene_Name"] = [x.split("_")[0] for x in df["ID"]]
        list_col = list(df)
        for col in ["Status", "Gene_Name"]:
            list_col.remove(col)
            list_col.insert(2, col)
        for col in ["Lipidation", "Glycosylation"]:
            list_col.remove(col)
            list_col.insert(len(list_col), col)
        return df[list_col]

    # One hot encode selected features
    def one_hot_encoding_categories(self, df, col_names=None):
        """Dummy encode features for given columns"""
        df.fillna({c: "None" for c in col_names}, inplace=True)
        enc_label = LabelEncoder()
        enc_onehot = OneHotEncoder(sparse=False)
        list_encoded = []
        list_labels = []
        list_col = []
        for col in col_names:
            cols = sorted(self.dict_col[col])   # single labels
            list_col.extend(cols)
            # integer encoded
            integer_encoded = enc_label.fit_transform(df[col].to_numpy())
            labels = enc_label.classes_.tolist()    # all labels (single & multi labels
            integer_encoded = integer_encoded.reshape((len(df), 1))
            # binary encoded
            onehot_encoded = enc_onehot.fit_transform(integer_encoded)
            # Add columns if not all labels are given
            col_dif = set(cols).difference(set(labels))
            n_dif = len(col_dif)
            if n_dif != 0:
                labels.extend(list(col_dif))
                dif_dummy = np.zeros((len(integer_encoded), n_dif))
                onehot_encoded = np.concatenate([onehot_encoded, dif_dummy], axis=1)
            list_labels.extend(labels)
            list_encoded.append(onehot_encoded)
        data = np.concatenate(list_encoded, axis=1)
        df = pd.DataFrame(data=data, columns=list_labels)
        df.drop(columns="None", inplace=True)
        # Add values from multi-labels to single labels (columns)
        for col in list(df):
            if ";" in col:
                a = np.array(df[col]).reshape(len(df), 1)
                df[col.split(";")] += a
        return df[list_col]

    @staticmethod
    def one_hot_encoding_integers(df, col_name=None, list_int=None, label_int="TM"):
        """Dummy encode features for given columns"""
        # Labels
        list_labels = ["{}{}".format(label_int, x) for x in list_int]
        # Encoding
        list_encoded = []
        for val in df[col_name]:
            list_val = []
            for x in list_int:
                if val == str(x):
                    list_val.append(1)
                else:
                    list_val.append(0)
            list_encoded.append(list_val)
        df = pd.DataFrame(data=list_encoded, columns=list_labels)
        return df

    # Main function
    def run_features(self, df):
        """Run feature modification pipeline"""
        # Transmembrane information
        list_tm_count = self._count_str(df, col_name="Transmembrane", str_count="TRANSMEM")
        list_tm_type = self._parse_column(df,
                                          col_name="Subcellular location [CC]",
                                          dict_items=self.dict_tm1_type,
                                          out_list=True)
        # Modify categorical data
        df = self._parse_column(df,
                                col_name="Subcellular location [CC]",
                                dict_items=self.dict_cc)
        df = self._parse_column(df,
                                col_name="Lipidation",
                                dict_items=self.dict_lip)
        df = self._parse_column(df,
                                col_name="Glycosylation",
                                dict_items=self.dict_gly)
        df = self._parse_column(df,
                                col_name="Post-translational modification",
                                dict_items=self.dict_ptm)
        # Modify continuous data
        loc_tm = df.columns.get_loc("Transmembrane")
        df.insert(loc_tm, "TM1_Type", list_tm_type)
        df.insert(loc_tm, "TM_n", list_tm_count)
        # Adjust column names
        df.drop(labels=["Transmembrane"], axis=1, inplace=True)
        dict_name = {"Entry": "ACC",
                     "Entry name": "ID",
                     "Post-translational modification": "PTM"}
        df = self._adjust_col_names(df, dict_name=dict_name)
        return df

    def run_one_hot_encoding(self, df):
        """One hot encode categorical features from df"""
        drop_list = ["Glycosylation", "Lipidation", "PTM", "Subcellular location [CC]"]
        list_cat = ["Glycosylation", "Lipidation", "PTM"]
        df_cc = self.one_hot_encoding_categories(df, col_names=["Subcellular location [CC]"])
        df_cat = self.one_hot_encoding_categories(df, col_names=list_cat)
        df_tm = self.one_hot_encoding_integers(df, col_name="TM1_Type",
                                               list_int=[1, 2, 3, 4],
                                               label_int="TM1_Type_")
        df_tm["Multi_TM"] = [1 if n > 1 else 0 for n in df["TM_n"]]
        # Ad multi TM "Multi_TM"
        df.drop(labels=drop_list, axis=1, inplace=True)
        df.drop(labels=["TM_n", "TM1_Type"], axis=1, inplace=True)
        list_df = [df, df_cc, df_tm, df_cat]
        df = pd.concat(list_df, axis=1)
        return df



