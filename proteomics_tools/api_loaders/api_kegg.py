"""
This is a script for creating a pathway and disease database from KEGG via REST API
"""
import io
import pandas as pd
import numpy as np
from bioservices import KEGG
from Bio.KEGG import REST

import proteomics_tools.api_loaders._utils as ut
from proteomics_tools.mapper import Mapper


# II Main Functions
class _KEGGDB(Mapper, ut.Utils):
    """Class for creating pathway and disease database for given organism"""
    def __init__(self, organism="mmu", reviewed=False):
        Mapper.__init__(self)
        ut.Utils.__init__(self)
        self.organism = organism
        self.reviewed = reviewed
        self.k = KEGG(verbose=False, )
        self.k.organism = self.organism
        self.columns = ["PATHWAY", "NAME", "ko_pathway", "n_proteins", "n_rel_pathways",
                        "protein_ids", "rel_pathways", "disease_ids", "drug_ids", "compound_ids", "module_ids"]
        self._id_columns = ["protein_ids", "disease_ids", "drug_ids", "compound_ids"]

    # Helper methods
    @staticmethod
    def _to_df(result):
        """Return a df, given tabular text"""
        return pd.read_table(io.StringIO(result), header=None)

    def _get_df_pathway_via_api(self):
        """Get pathway data for given organism via biopyghon REST API"""
        result = REST.kegg_list("pathway", self.organism).read()
        df_pathway = self._to_df(result)
        df_pathway.columns = ["pathways", "name"]
        return df_pathway

    # Convert list of strings into list
    def _df_str_to_list(self, df=None, out="list", columns=None):
        """Convert string list in df to lists or sets"""
        if columns is None:
            columns = self._id_columns
        if type(columns) is str:
            columns = [columns]
        list_out = ["list", "set"]
        if out not in list_out:
            raise ValueError("'out' must be in {}".format(list_out))
        if out == "list":
            for col in columns:
                df[col] = [self.str_to_list(i) if str(i) != "nan" else [] for i in df[col]]
        else:
            for col in columns:
                df[col] = [set(self.str_to_list(i)) if str(i) != "nan" else set() for i in df[col]]
        return df

    def _ko_to_up_ids(self, list_kos=None, dict_keggid_upid=None):
        """Convert list of ko into uniprot ids"""
        org_str = self.organism.upper()
        list_upids = []
        for ko in list_kos:
            dict_res = self.k.parse(self.k.get(ko))
            if "GENES" in dict_res.keys():
                if org_str in dict_res["GENES"].keys():
                    kegg_id = "{}:{}".format(self.organism, dict_res["GENES"][org_str])
                    if "(" in kegg_id:
                        kegg_id = kegg_id.split("(")[0]
                    try:
                        list_upids.append(dict_keggid_upid[kegg_id])
                    except KeyError:
                        pass
        return list_upids

    @staticmethod
    def _get_all_feat_items(df_pathway=None, feat_name=None):
        """Get all items for given feature from df pathway"""

        feat_items = []
        for i, row in df_pathway.iterrows():
            feat_item = row["{}_ids".format(feat_name)]
            if str(feat_item) != "nan":
                feat_items.extend(feat_item)
        feat_items = list(set(feat_items))
        feat_items.sort()
        return feat_items

    # 1. Get pathway with related pathways
    def get_df_pathways(self, verbose=False):
        """Get complex pathway df with information on related pathways"""
        # Perform the query
        df_pathway = self._get_df_pathway_via_api()
        list_pathways = list(df_pathway["pathways"])
        dict_keggid_upid = self.mapping_dict(fr="KEGG_ID", to="ACC", reviewed=self.reviewed)
        kegg_org = "{}:".format(self.organism)
        list_row = []
        for i, p in enumerate(list_pathways):
            if verbose:
                print(f"\rInformation for {i+1} out of {len(list_pathways)} pathways are retrieved", end="")
            result = REST.kegg_get(p).read()
            dict_result = self.k.parse(result)
            name_full = dict_result["NAME"][0].split(" - ")
            name = " - ".join(name_full[:-1])
            # Add genes
            if "GENE" in dict_result.keys():
                genes = []
                for i in dict_result["GENE"].keys():
                    try:
                        up_id = dict_keggid_upid[kegg_org + i]
                        genes.append(up_id)
                    except KeyError:
                        pass
                n_genes = len(genes)
            else:
                genes = np.NAN
                n_genes = 0
            ko_pathway = dict_result["KO_PATHWAY"]
            result = REST.kegg_get(ko_pathway).read()
            dict_result_ko = self.k.parse(result)
            # Get related pathways
            if "REL_PATHWAY" in dict_result.keys():
                list_related_pathways = [i for i in self.str_to_list(dict_result["REL_PATHWAY"], " ")
                                         if self.organism in i]
                n_related_pathways = len(list_related_pathways)
            else:
                list_related_pathways = np.NAN
                n_related_pathways = 0
            row = [p, name, ko_pathway, n_genes, n_related_pathways,
                   self.list_to_str(genes),
                   self.list_to_str(list_related_pathways)]
            row = self.add_infos(row=row, dict_data=dict_result_ko, columns=["DISEASE", "DRUG", "COMPOUND", "MODULE"])
            list_row.append(row)
        df = pd.DataFrame(list_row, columns=self.columns)
        df.sort_values(["n_rel_pathways", "n_proteins"], ascending=False, inplace=True)
        if verbose:
            print("\n")
        return df

    # 2. Merge proteins of related pathways
    def get_modules(self, df_pathway=None):
        """Get df for each module and list of protein ids for given organism"""
        dict_keggid_upid = self.mapping_dict(fr="KEGG_ID", to="ACC", reviewed=self.reviewed)
        # Get List of modules
        list_modules = []
        for m in df_pathway["module_ids"]:
            if str(m) != "nan":
                list_modules.extend(self.str_to_list(m))
        list_modules = list(set(list_modules))
        list_modules.sort()
        # Get list of genes for each module
        list_row = []
        for m in list_modules:
            dict_res = self.k.parse(self.k.get(m))
            str_kos = dict_res["DEFINITION"]
            str_kos = str_kos.replace("(", "").replace(")", "").replace(",", " ").replace("+", " ").replace("-", " ")
            list_kos = str_kos.split(" ")
            list_upids = self._ko_to_up_ids(list_kos=list_kos,
                                            dict_keggid_upid=dict_keggid_upid)
            str_up_ids = self.list_to_str(list_upids)
            str_kos = self.list_to_str(list_kos)
            n_protein_ids = len(list_upids)
            n_ko_ids = len(list_kos)
            list_row.append([m, n_protein_ids, n_ko_ids, str_up_ids, str_kos])
        columns = ["MODULE", "n_protein_ids", "n_ko_ids", "protein_ids", "ko_ids"]
        df_pathway = pd.DataFrame(list_row, columns=columns)
        df_pathway.sort_values(by="MODULE", inplace=True)
        return df_pathway

    def modules_to_pathway(self, df_pathway=None, df_module=None):
        """Join protein ids from modules into pathways"""
        df_pathway = self._df_str_to_list(df=df_pathway, out="set", columns="protein_ids")
        df_pathway = self._df_str_to_list(df=df_pathway, out="list", columns="module_ids")
        df_module = self._df_str_to_list(df=df_module, out="set", columns="protein_ids")
        list_proteins = []
        list_n_proteins = []
        for i, row in df_pathway.iterrows():
            proteins = row["protein_ids"]
            modules = row["module_ids"]
            for m in modules:
                proteins_module = list(df_module[df_module["MODULE"] == m]["protein_ids"])[0]
                proteins = proteins.union(proteins_module)
            list_proteins.append(self.list_to_str(list(proteins)))
            list_n_proteins.append(len(proteins))
        df_pathway["protein_ids"] = list_proteins
        df_pathway["n_proteins"] = list_n_proteins
        df_pathway["module_ids"] = [self.list_to_str(i) for i in df_pathway["module_ids"]]
        return df_pathway

    def get_df_feat(self, df_pathway=None, feat_name="disease", verbose=False):
        """Get df disease or drug with associated uniprot accession numbers"""
        dict_keggid_upid = self.mapping_dict(fr="KEGG_ID", to="ACC",  reviewed=self.reviewed),
        df_pathway = self._df_str_to_list(df=df_pathway, out="set", columns="{}_ids".format(feat_name))
        feat_items = self._get_all_feat_items(df_pathway=df_pathway, feat_name=feat_name)
        list_rows = []
        for i, item in enumerate(feat_items):
            if verbose:
                print(f"\rInformation for {i+1} out of {len(feat_items)} {feat_name} are retrieved", end="")
            result = REST.kegg_get(item).read()
            dict_result = self.k.parse(result)
            name = dict_result["NAME"][0]
            if ";" in name:
                name = name.replace(";", "")
            category = dict_result.get("CATEGORY", np.NAN)
            if "GENE" in dict_result.keys():
                gene_list = list(dict_result["GENE"].values())
                list_kos = [i.split("KO:")[1].replace("]", "") for i in gene_list if "KO:" in i]
                list_upids = self._ko_to_up_ids(list_kos=list_kos,
                                                dict_keggid_upid=dict_keggid_upid)
            else:
                list_kos = []
                list_upids = []
            n_protein_ids = len(list_upids)
            n_ko_ids = len(list_kos)
            list_rows.append([item, name, category,
                              n_protein_ids, n_ko_ids,
                              self.list_to_str(list_upids),
                              self.list_to_str(list_kos)])
        columns = [feat_name.upper(), "NAME", "CATEGORY", "n_protein_ids", "n_ko_ids", "protein_ids", "ko_ids"]
        df = pd.DataFrame(list_rows, columns=columns)
        if verbose:
            print("\n")
        return df


class KEGGDB:
    """Class to retrieve pathway and disease information from KEGG via REST API and
    parse it into df format"""
    def __init__(self, organism="mmu", reviewed=False):
        self._organism = organism
        self._reviewed = reviewed

    def get_pathway(self, include_ids_from_module=True, verbose=False):
        """Get KEGG pathway via Biopython REST API and save as df"""
        kd = _KEGGDB(organism=self._organism, reviewed=self._reviewed)
        df_pathway = kd.get_df_pathways(verbose=verbose)
        if include_ids_from_module:
            df_module = kd.get_modules(df_pathway=df_pathway)
            df_pathway = kd.modules_to_pathway(df_pathway=df_pathway,
                                               df_module=df_module)
        return df_pathway

    def get_disease(self, df_pathway=None, verbose=False):
        """Get KEGG disease with associated protein ids via Biopython REST API and save as df"""
        kd = _KEGGDB(organism=self._organism, reviewed=self._reviewed)
        df_disease = kd.get_df_feat(df_pathway=df_pathway, feat_name="disease", verbose=verbose)
        return df_disease
