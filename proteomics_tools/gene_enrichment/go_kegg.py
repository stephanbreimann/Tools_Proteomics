"""
This is a script for KEGG analysis (check enrichment of genes in pathway or disease)
"""
import numpy as np
import pandas as pd
from bioservices import KEGG
from scipy.stats import fisher_exact

import proteomics_tools.gene_enrichment._utils as ut
from proteomics_tools.mapper import Mapper


# I Helper Functions
def _retriever(df=None, query="PATHWAY"):
    """Abstract retiever function for dict_pathway and dict_disease"""
    if query in list(df):
        df = df.set_index(query)
    dict_query = df.to_dict(orient="index")
    # Convert id strings into id lists
    for key_outer in dict_query:
        for key in dict_query[key_outer]:
            if "ids" in key and "n_" not in key:
                dict_query[key_outer][key] = ut.Utils().str_to_list(dict_query[key_outer][key])
    return dict_query


def _dict_upid_retriever(df=None, query="PATHWAY"):
    """Get dict for UniProt id to pathway or disease"""
    list_up_ids = []
    dict_query_up_ids = {}
    for i, row in df.iterrows():
        up_ids_str = row["protein_ids"]
        query_val = row[query]
        if str(up_ids_str) != "nan":
            up_ids = ut.Utils().str_to_list(up_ids_str)
            dict_query_up_ids[query_val] = up_ids
            list_up_ids.extend(up_ids)
        else:
            dict_query_up_ids[query_val] = []
    list_up_ids = list(set(list_up_ids))
    dict_upid_query = {up_id: [] for up_id in list_up_ids}
    for query_val in dict_query_up_ids:
        for up_ids in dict_query_up_ids[query_val]:
            dict_upid_query[up_ids].append(query_val)
    return dict_upid_query


# II Main Functions
# 1. Preload Kegg
class KeggLoader:
    """Load all data for KEGG analysis"""
    def __init__(self, organism="hsa", reviewed=True, folder_data=None):
        self.reviewed = reviewed
        self.organism = organism
        if folder_data is None:
            folder_data = ut.FOLDER_DATA + "GO_KEGG/"
        self.folder_data = folder_data

    # Set general dictionaries (for pathway and disease)
    def run(self):
        """Load kegg_pathway and kegg_disease"""
        file_disease = self.folder_data + "kegg_disease_{}".format(self.organism)
        file_pathway = self.folder_data + "kegg_pathway_{}".format(self.organism)
        kegg_disease = pd.read_csv(file_disease + ".tab", sep="\t")
        kegg_pathway = pd.read_csv(file_pathway + ".tab", sep="\t")
        # Dict for data from database
        dict_pathway = _retriever(df=kegg_pathway, query="PATHWAY")
        dict_disease = _retriever(df=kegg_disease, query="DISEASE")
        # Dict for upid id and ID fir pathway and disease
        dict_upid_pathway = _dict_upid_retriever(df=kegg_pathway, query="PATHWAY")
        dict_upid_disease = _dict_upid_retriever(df=kegg_disease, query="DISEASE")
        # Mappers
        mp = Mapper(organism=self.organism)
        dict_keggid_upid = mp.mapping_dict(fr="KEGG_ID", to="ACC", reviewed=self.reviewed)
        dict_upid_keggid = mp.mapping_dict(fr="ACC", to="KEGG_ID", reviewed=self.reviewed)
        dict_upid_name = mp.mapping_dict(fr="ACC", to="Gene_Name", reviewed=self.reviewed)
        dict_kegg = dict(kegg_pathway=kegg_pathway,
                         kegg_disease=kegg_disease,
                         dict_pathway=dict_pathway,
                         dict_disease=dict_disease,
                         dict_upid_pathway=dict_upid_pathway,
                         dict_upid_disease=dict_upid_disease,
                         dict_keggid_upid=dict_keggid_upid,
                         dict_upid_keggid=dict_upid_keggid,
                         dict_upid_name=dict_upid_name)
        return dict_kegg


# 2. Kegg Analysis
class KeggSetter:
    """Class for setting dictionaries used during processing"""
    def __init__(self, organism="hsa", reviewed=True, folder_data=None, dict_kegg=None):
        self.reviewed = reviewed
        self.organism = organism
        if folder_data is None:
            folder_data = ut.FOLDER_DATA + "GO_KEGG/"
        if dict_kegg is None:
            dict_kegg = dict()
        self.folder_data = folder_data
        # Database for pathway and disease
        self.kegg_pathway = dict_kegg.get("kegg_pathway", None)
        self.kegg_disease = dict_kegg.get("kegg_disease", None)
        self._set_kegg()
        # Dict for data from database
        self.dict_pathway = dict_kegg.get("dict_pathway", None)
        self.dict_disease = dict_kegg.get("dict_disease", None)
        # Dict for upid id and ID fir pathway and disease
        self.dict_upid_pathway = dict_kegg.get("dict_upid_pathway", None)
        self.dict_upid_disease = dict_kegg.get("dict_upid_disease", None)
        # Mappers
        self.dict_keggid_upid = dict_kegg.get("dict_keggid_upid", None)
        self.dict_upid_keggid = dict_kegg.get("dict_upid_keggid", None)
        self.dict_upid_name = dict_kegg.get("dict_upid_name", None)

    # Set general dictionaries (for pathway and disease)
    def _set_kegg(self):
        """Load kegg_pathway and kegg_disease"""
        file_disease = self.folder_data + "kegg_disease_{}".format(self.organism)
        file_pathway = self.folder_data + "kegg_pathway_{}".format(self.organism)
        if self.kegg_pathway is None:
            self.kegg_pathway = pd.read_csv(file_pathway + ".tab", sep="\t")
        if self.kegg_disease is None:
            self.kegg_disease = pd.read_csv(file_disease + ".tab", sep="\t")

    def _set_dict_pathway(self):
        """Set dict for pathway information"""
        if self.dict_pathway is None:
            self.dict_pathway = _retriever(df=self.kegg_pathway, query="PATHWAY")

    def _set_dict_disease(self):
        """Set dict for disease information"""
        if self.dict_disease is None:
            self.dict_disease = _retriever(df=self.kegg_disease, query="DISEASE")

    # Dict upid to pathway and disease ids
    def _set_dict_upid_pathway(self):
        """Set dict for UniProt id to pathway"""
        if self.dict_upid_pathway is None:
            self.dict_upid_pathway = _dict_upid_retriever(df=self.kegg_pathway, query="PATHWAY")

    def _set_dict_upid_disease(self):
        """Set dict for UniProt id to pathway"""
        if self.dict_upid_disease is None:
            self.dict_upid_disease = _dict_upid_retriever(df=self.kegg_disease, query="DISEASE")

    # Mappers
    def _set_dict_keggid_upid(self):
        """Set dict for KEGG id for proteins to uniprot id (accession number)"""
        if self.dict_keggid_upid is None:
            mp = Mapper(organism=self.organism)
            self.dict_keggid_upid = mp.mapping_dict(fr="KEGG_ID", to="ACC",
                                                    reviewed=self.reviewed)

    def _set_dict_upid_keggid(self):
        """Set dict for uniprot id (accession number) to KEGG id for proteins"""
        if self.dict_upid_keggid is None:
            mp = Mapper(organism=self.organism)
            self.dict_upid_keggid = mp.mapping_dict(fr="ACC", to="KEGG_ID",
                                                    reviewed=self.reviewed)

    def _set_dict_upid_name(self):
        """Set dict for uniprot id (accession number) to gene name"""
        if self.dict_upid_name is None:
            mp = Mapper(organism=self.organism)
            self.dict_upid_name = mp.mapping_dict(fr="ACC", to="Gene_Name",
                                                  reviewed=self.reviewed)


class KeggGoBase(KeggSetter):
    """Class for KEGG pathway analysis using bioservice package"""

    def __init__(self, organism="hsa", dict_kegg=None):
        KeggSetter.__init__(self, organism=organism, dict_kegg=dict_kegg)
        self.kegg = KEGG(verbose=False)
        self.kegg.organism = self.organism
        self.list_dbs = ["PATHWAY", "DISEASE"]

    # Helper methods
    def _get_study_name(self, study_item=None):
        """Get protein name of study_name"""
        self._set_dict_upid_name()
        try:
            study_name = self.dict_upid_name[study_item]
        except KeyError:
            try:
                mp = Mapper(organism=self.organism)
                study_name = mp.fr_to(fr="ACC", to="Protein_name", query=study_item)
            except:
                study_name = np.NAN
        return study_name

    # General enrichment method
    def _get_list_query(self, study=None, query="PATHWAY"):
        """Get list of pathways in which study items occur"""
        if query == "PATHWAY":
            self._set_dict_upid_pathway()
            dict_upid_query = self.dict_upid_pathway
        else:
            self._set_dict_upid_disease()
            dict_upid_query = self.dict_upid_disease
        list_query = []
        for up_id in study:
            if up_id in dict_upid_query.keys():
                pathways = dict_upid_query[up_id]
                list_query.extend(pathways)
        list_query = list(set(list_query))
        return list_query

    @staticmethod
    def _enrichment(study=None, pop=None, dict_data=None):
        """Enrichment analysis for diseases and pathways"""
        study_n = len(study)
        pop_n = len(pop)
        genes = ut.Utils().str_to_list(dict_data["protein_ids"])
        if len(genes) > 0:
            genes_count = len(genes)
            study_items = list(set(study).intersection(set(genes)))
            study_count = len(study_items)
            pop_count = len(set(pop).intersection(set(genes)))
            ratio_study = "{}/{}".format(study_count, study_n)
            ratio_pop = "{}/{}".format(pop_count, pop_n)
            if pop_count != 0:
                fold_enrichment = round((study_count/study_n)/(pop_count/pop_n), 2)
            else:
                fold_enrichment = np.NAN
            p_uncorrected = fisher_exact([[study_count, study_n], [pop_count, pop_n]])[1]
            freq_in_study = round(study_count / study_n * 100, 1)
            freq_study_in_total = round(study_count / genes_count * 100, 1)
        else:
            genes = study_items = np.NAN
            study_count = pop_count = genes_count = ratio_study = ratio_pop = np.NAN
            fold_enrichment = p_uncorrected = freq_in_study = freq_study_in_total = np.NAN
        if study_items is not np.NAN and len(study_items) == 0:
            study_items = np.NAN
        result = [ratio_study, ratio_pop, fold_enrichment, p_uncorrected,
                  freq_in_study, freq_study_in_total,
                  study_count, pop_count, genes_count,
                  ut.Utils().list_to_str(study_items), ut.Utils().list_to_str(genes)]
        return result

    # Pathway enrichment
    def _enrichment_df(self, study=None, pop=None, query="PATHWAY"):
        """Create df with results of enrichment analysis for pathway"""
        if query not in self.list_dbs:
            raise ValueError("query must be in {}".format(self.list_dbs))
        if query == "PATHWAY":
            self._set_dict_pathway()
            dict_query = self.dict_pathway
            add_info = "DISEASE"
        else:
            self._set_dict_disease()
            dict_query = self.dict_disease
            add_info = "PATHWAY"
        add_info_col = [add_info, "DRUG"]
        str_add_info_ids = "{}_ids".format(add_info.lower())
        list_query = self._get_list_query(study=study, query=query)
        list_rows = []
        for query_val in list_query:
            dict_data = dict_query[query_val]
            name = dict_data["NAME"]
            row = [query_val, name]
            # Add enrichment results
            res = self._enrichment(dict_data=dict_data, study=study, pop=pop)
            row.extend(res)
            # Add data for drugs and diseases if available
            row = ut.Utils().add_infos(row=row, dict_data=dict_data, columns=add_info_col)
            list_rows.append(row)
        columns = [query, "NAME",
                   "ratio_in_study", "ratio_in_pop",
                   "fold_enrichment", "p_uncorrected",
                   "study_%", "study_in_proteins_%",
                   "study_count", "pop_count", "protein_count",
                   "study_items", "protein_ids", str_add_info_ids, "drug_ids"]
        df = pd.DataFrame(list_rows, columns=columns)
        return df


class KeggGo(KeggGoBase):
    """Class for KEGG pathway analysis using bioservice package for retrieving information from the KEGG
     database via REST API
     """

    def __init__(self, organism="hsa", dict_kegg=None, verbose=False):
        KeggGoBase.__init__(self, organism=organism, dict_kegg=dict_kegg)
        self.verbose = verbose

    # Load all necessary data
    def filter_study_pop(self, study=None, pop=None):
        """Remove all study and pop items that are not in KEGG"""
        df = self.kegg_pathway
        protein_ids = []
        for i in df["protein_ids"]:
            protein_ids.extend(ut.Utils().str_to_list(i))
        protein_ids = list(set(protein_ids))
        study_kegg = [i for i in study if i in protein_ids]
        pop_kegg = [i for i in pop if i in protein_ids]
        # Give information about ratios in KEGG via UserWarning
        ratio_study_in_kegg = round(len(study_kegg)/len(study) * 100, 2)
        ratio_pop_in_kegg = round(len(pop_kegg)/len(pop) * 100, 2)
        if self.verbose:
            warn_str = "\n{}% of study items are in KEGG (n={})".format(ratio_study_in_kegg,
                                                                        len(study_kegg))
            print(warn_str)
            warn_str = "\n{}% of population items (background) are in KEGG (n={})".format(ratio_pop_in_kegg,
                                                                                          len(pop_kegg))
            print(warn_str)
        return study_kegg, pop_kegg

    # Pathway enrichment
    def pathway_enrichment(self, study=None, pop=None, p_method="bonferroni"):
        """Enrichment analysis for pathways"""
        p_method = ut.Utils().check_p_method(p_method=p_method)
        study, pop = self.filter_study_pop(study=study, pop=pop)
        df = self._enrichment_df(study=study, pop=pop, query="PATHWAY")
        # Correct p value
        df = ut.Utils().correct_p_values(df=df, p_method=p_method)
        return df

    # Disease enrichment
    def disease_enrichment(self, study=None, pop=None, p_method="bonferroni"):
        """Enrichment analysis for diseases"""
        p_method = ut.Utils().check_p_method(p_method=p_method)
        study, pop = self.filter_study_pop(study=study, pop=pop)
        df = self._enrichment_df(study=study, pop=pop, query="DISEASE")
        # Correct p value
        df = ut.Utils().correct_p_values(df=df, p_method=p_method)
        return df

    # Study items
    def study_items_enrichment(self, study=None, df_enrich=None, df_pathway=None, df_disease=None):
        """Summary of enrichment analysis for given study items"""
        list_ids = ["GO", "PATHWAY", "DISEASE"]
        list_rows = []
        for study_item in study:
            study_name = self._get_study_name(study_item=study_item)
            row = [study_item, study_name]
            for i, df in enumerate([df_enrich, df_pathway, df_disease]):
                df_filtered = df[df["study_items"].str.contains(study_item)]
                count = len(df_filtered)
                ids = ut.Utils().list_to_str(list(df_filtered[list_ids[i]]))
                names = ut.Utils().list_to_str(list(df_filtered["NAME"]))
                row.extend([count, ids, names])
            list_rows.append(row)
        columns = ["STUDY", "NAME", "GO_count", "GO_ids", "GO_names",
                   "PATHWAY_count", "PATHWAY_ids", "PATHWAY_names",
                   "DISEASE_count", "DISEASE_ids", "DISEASE_names"]
        df = pd.DataFrame(list_rows, columns=columns)
        df = df[["STUDY", "NAME", "GO_count", "PATHWAY_count", "DISEASE_count",
                 "GO_ids", "PATHWAY_ids", "DISEASE_ids",
                 "GO_names", "PATHWAY_names", "DISEASE_names"]]
        df = df.sort_values(by=["PATHWAY_count", "DISEASE_count", "GO_count"], ascending=False).reset_index(drop=True)
        return df

    def show_pathways(self, df_pathway=None, n=5):
        """Show best n pathways in KEGG"""
        self._set_dict_upid_keggid()
        for i, row in df_pathway.head(n).iterrows():
            study_items = ut.Utils().str_to_list(row["study_items"])
            dict_keggid_color = {self.dict_upid_keggid[i]: "red,black" for i in study_items}
            pathway = row["KEGG"]
            self.kegg.show_pathway(pathway, keggid=dict_keggid_color)
        # TODO Adjust coloring according to enrichment

    def run(self, study=None, pop=None, p_method=None, df_enrich_filtered=None):
        """Complete KEGG enrichment analysis using bioservice package for retrieving information from the KEGG
        database via REST API for pathways and diseases.
        In: a) study: list of study items (UniProt accession numbers)
            b) pop: list of population items (UniProt accession numbers
            c) p_method: method for p-value correction ["bonferroni", "sidak", "holm", "fdr"]
            d) df_enrich_filtered: df with enriched proteins
        Out:a) df_pathway: df with pathways for which study items are enriched
            b) df_disease: df with diseases for which study items are enriched
            c) df_study: df with study items and all associated information from KEGG
        References
        ----------
        Cokelaer et al. BioServices: a common Python package to access biological Web Services
            programmatically Bioinformatics (2013) 29 (24): 3241-3242
        """
        dict_enrich = {"study": study, "pop": pop, "p_method": p_method}
        df_pathway = self.pathway_enrichment(**dict_enrich)
        df_disease = self.disease_enrichment(**dict_enrich)
        df_study = self.study_items_enrichment(study=study,
                                               df_enrich=df_enrich_filtered,
                                               df_pathway=df_pathway,
                                               df_disease=df_disease)
        return df_pathway, df_disease, df_study

