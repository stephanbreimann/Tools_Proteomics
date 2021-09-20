"""
This is a script for gene enrichment using python package goatools.
Ref.: https://github.com/tanghaibao/goatools
"""
import os
import pandas as pd
from goatools.go_enrichment import GOEnrichmentStudy
import sys

from proteomics_tools.gene_enrichment.go_base import BaseGO
import proteomics_tools.gene_enrichment._utils as ut


# I Helper Functions
class EnrichmentBase(BaseGO):
    """Class for loading GO DB and to query GO terms"""

    def __init__(self, dict_go=None, df_org_go=None):
        BaseGO.__init__(self, dict_go=dict_go)
        self.df_org_go = df_org_go

    # Helper methods
    @staticmethod
    def _block_print():
        sys.stdout = open(os.devnull, 'w')

    @staticmethod
    def _enable_print():
        sys.stdout = sys.__stdout__

    @staticmethod
    def _get_up_ids(df=None, up_col="DB_Object_ID"):
        """Get uniprot ids """
        if up_col in list(df):
            list_up_ids = list(set(df[up_col]))
        else:
            raise ValueError("'up_col' missing in df")
        list_up_ids = [i for i in list_up_ids if i is not None]
        return list_up_ids

    def _get_dict_upid_goid(self, df_org_go=None):
        """Get uniprot and go ids from index (u_id|go_id) and convert into dict"""
        list_up_ids = self._get_up_ids(df=df_org_go)
        dict_upid_goid = {x: [] for x in list_up_ids}  # Uniprot ID to GO id
        for u_id_go_id in list(df_org_go.index):
            u_id, go_id = u_id_go_id.split("|")
            dict_upid_goid[u_id].append(go_id)
        dict_upid_goid = {x: set(dict_upid_goid[x]) for x in dict_upid_goid}
        return dict_upid_goid

    # Enrichment
    def _go_enrichment(self, study=None, pop=None, block_print=True):
        """GO term enrichment using goatools and implemented as in DAVID (https://david.ncifcrf.gov/home.jsp)"""
        dict_upid_goid = self._get_dict_upid_goid(df_org_go=self.df_org_go)
        # GO enrichment study (Run fishers exact test)
        if block_print:
            self._block_print()
        g = GOEnrichmentStudy(pop, dict_upid_goid, self.dict_go,
                              propagate_counts=True,
                              alpha=0.05,
                              methods=["bonferroni"])
        g_res = g.run_study(study)
        # Save temporary file
        file_tmp = ut.FOLDER_DATA + "tmp_gene_enrichment.xlsx"
        g.wr_xlsx(file_tmp, goea_results=g_res)
        df_go_enrich = pd.read_excel(file_tmp)
        if block_print:
            self._enable_print()
        # Add fold enrichment from ratio information to df (as calculated by DAVID)
        loc_ratio_in_pop = df_go_enrich.columns.get_loc("ratio_in_pop")
        list_ratio_in_study = self._get_ratio(df=df_go_enrich, col="ratio_in_study")
        list_ratio_in_pop = self._get_ratio(df=df_go_enrich, col="ratio_in_pop")
        fold_enrichment = [round(s / list_ratio_in_pop[i], 2) for i, s in enumerate(list_ratio_in_study)]
        df_go_enrich.insert(loc_ratio_in_pop + 1, "fold_enrichment", fold_enrichment)
        os.remove(file_tmp)
        # Filter
        df_go_enrich.drop("p_bonferroni", inplace=True, axis=1)
        df_go_enrich.dropna(axis=0, subset=["GO"], inplace=True)
        df_go_enrich.rename(columns={"name": "NAME"}, inplace=True)
        return df_go_enrich


# II Main Functions
class EnrichmentGo(EnrichmentBase):
    """Class for loading GO DB and to query GO terms using python package goatools and
    computation of the the fold enrichment as in DAVID"""
    def __init__(self, dict_go=None, df_org_go=None, verbose=False):
        EnrichmentBase.__init__(self, dict_go=dict_go, df_org_go=df_org_go)
        self.verbose = verbose

    # Quality check of study and pop sets
    def filter_study_pop(self, study=None, pop=None, up_col="DB_Object_ID"):
        """Check if study is subset of pop is subset of up ids in df_org_go"""
        list_up_ids = self._get_up_ids(df=self.df_org_go, up_col=up_col)
        if not set(study).issubset(set(pop)):
            pop = list(set(study).union(set(pop)))
            if self.verbose:
                print("study set is not subset of pop set")
        if not set(pop).issubset(set(list_up_ids)):
            dif = list(set(pop).difference(set(list_up_ids)))
            if self.verbose:
                warn_str = "Following up ids are removed from pop and study because no GO term exists\n: {}".format(dif)
                print(warn_str)
            pop = [p for p in pop if p not in dif]
            study = [x for x in study if x in pop]
        return study, pop

    # Main method
    def run(self, study=None, pop=None, p_method=None, block_print=True):
        """Enrichment analysis, where e (enriched) and p (purified) means  significantly higher resp. lower.
        GO terms enrichment is performed using the python package goatools. Additionally, the fold enrichment
        is computed as in DAVID.
        In: a) study: list with uniprot ids for study set
            b) pop: list with uniprot ids for population (background)
            c) df_org_go: df with all uniprot ids that are associated with a GO term for given organism
            d) p_method: str or list with str for p value corrections ["bonferroni", "sidak", "holm", "fdr"]
        Out:a) df_go_enrich: result of enrichment analysis
        References
        ----------
        Klopfenstein DV, Zhang L, Pedersen BS, ... Tang H GOATOOLS:
        A Python library for Gene Ontology analyses Scientific reports | (2018) 8:10872 | DOI:10.1038/s41598-018-28948-z
        Nature Protocols 2009; 4(1):44 & Nucleic Acids Res. 2009;37(1):1
        """
        p_method = self.check_p_method(p_method=p_method)
        study, pop = self.filter_study_pop(study=study, pop=pop)
        df = self._go_enrichment(study=study, pop=pop, block_print=block_print)
        # Correct p value
        df = self.correct_p_values(df=df, p_method=p_method)
        df_go_enrichment = df
        return df_go_enrichment
