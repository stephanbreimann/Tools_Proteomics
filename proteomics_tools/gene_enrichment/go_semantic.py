"""
This is a script for a semantic analysis of GO terms using the recursive clustering algorithm of REVIGO
Supek et al., 2011 (REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms)
            https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0021800
"""
import math
import numpy as np
import pandas as pd

from proteomics_tools.gene_enrichment.go_base import BaseGO


# I Helper Functions
class SemanticGoBase(BaseGO):
    """Class for loading GO DB and to query GO terms"""

    def __init__(self, dict_go=None, termcounts=None):
        BaseGO.__init__(self, dict_go=dict_go)
        self.termcounts = termcounts

    # Helper methods
    @staticmethod
    def _check_namespace(ns=None):
        """Check if namespace valid"""
        list_ns = ["BP", "MF", "CC"]
        if ns not in list_ns:
            raise ValueError("'ns' must be in: {}".format(list_ns))

    @staticmethod
    def _check_single_p_method(p_method=None, df=None):
        """Check if p_method is valid"""
        list_p_method = ["bonferroni", "sidak", "holm", "fdr"]
        if p_method is None:
            p_method = "bonferroni"
        if type(p_method) is list:
            p_method = p_method[0]
        elif type(p_method) is not str or p_method not in list_p_method:
            raise ValueError("'p_method' must be single value from {}".format(list_p_method))
        col_p_method = "p_" + p_method
        if col_p_method not in list(df):
            raise ValueError("'p_method' not in given df")
        return p_method

    # Semantic similarity
    def _min_branch_length(self, go_ids=None):
        """Finds the minimum branch length between two terms in the GO DAG"""
        # First get the deepest common ancestor
        ncp = self.nearest_common_parent(go_ids=go_ids)
        # Then get the distance from the DCA to each term
        ncp_depth = self.dict_go[ncp].depth
        d1 = self.dict_go[go_ids[0]].depth - ncp_depth
        d2 = self.dict_go[go_ids[1]].depth - ncp_depth
        # Return the total distance - i.e., to the deepest common ancestor and back.
        return d1 + d2

    def semantic_distance_edge_based(self, go_ids=None):
        """Finds the semantic distance (minimum number of connecting branches)
        between two GO terms."""
        sem_dist = self._min_branch_length(go_ids=go_ids)
        return sem_dist

    def semantic_similarity_edge_based(self, go_ids=None):
        """Finds the semantic similarity (inverse of the semantic distance) between two GO terms."""
        sem_sim = 1.0 / float(self.semantic_distance_edge_based(go_ids=go_ids))
        return sem_sim

    # Four different semantic similarity measures
    # Implementation adjusted from NLTK: https://www.nltk.org/_modules/nltk/corpus/reader/wordnet.html
    def _information_content(self, go_id=None):
        """Calculates the information content (i.e.-log("freq of GO term")) of a GO term."""
        # Get the observed frequency of the GO term
        freq = self.termcounts.get_term_freq(go_id=go_id)
        ic = -1.0 * math.log(freq)  # Calculate the information content (-log(freq)) to base e
        return ic

    def _res_similarity(self, go_ids=None):
        """Resnik's similarity measure (i.e. information content of nearest common parent)
        E.g. resnik similarity for 'biological process' = 0 because it occurs in every go term as parent
        Ref: Resnik, 1995 (Using Information Content to Evaluate Semantic Similarity in a Taxonomy)"""
        ncp = self.nearest_common_parent(go_ids=go_ids)     # called least common subsumer (lcs) in WordNet
        ic_ncp = self._information_content(go_id=ncp)
        res_sim = ic_ncp
        return res_sim

    def _jcn_similarity(self, go_ids=None):
        """Compute Jiang Conrath similarity based on resnik similarity. Considers the information content of
        nearst common parent (nsc) and the two compared concepts to calculate the distance between the two concepts.
        Ref: Jiang & Conrath, 1997 (Semantic Similarity Based on Corpus Statistics and Lexical Taxonomy)"""
        ic_a = self._information_content(go_id=go_ids[0])
        ic_b = self._information_content(go_id=go_ids[1])
        ic_ncp = self._res_similarity(go_ids=go_ids)
        if ic_a == 0 or ic_b == 0:
            return 0
        # Distance function (as in https://www.nltk.org/_modules/nltk/corpus/reader/wordnet.html)
        ic_difference = ic_a + ic_b - 2 * ic_ncp
        if ic_difference == 0:
            return math.inf
        jcn_sim = 1/ic_difference   # similarity as reciprocal of distance
        return jcn_sim

    def _lin_similarity(self, go_ids=None):
        """Compute Lin similarity based on resnik similarity. Considers the information content of
        lowest common subsumer (lcs) and the two compared concepts.
        Ref: Lin, 1998 (An Information-Theoretic Definition of Similarity)"""
        ic_a = self._information_content(go_id=go_ids[0])
        ic_b = self._information_content(go_id=go_ids[1])
        ic_ncp = self._res_similarity(go_ids=go_ids)
        lin_sim = (2.0 * ic_ncp) / (ic_a + ic_b)
        return lin_sim

    def _sim_rel(self, go_ids=None):
        """Compute semantic similarity according to Pesquita. Based on weighted Lin similarity.
        Ref: Pesquita et al., 2009 (Semantic Similarity in Biomedical Ontologies)"""
        ncp = self.nearest_common_parent(go_ids=go_ids)
        freq_ncp = self.termcounts.get_term_freq(go_id=ncp)
        lin_sim = self._lin_similarity(go_ids=go_ids)
        sim_rel = lin_sim * (1-freq_ncp)
        return sim_rel

    @staticmethod
    def _check_similarity_measure(sim_measure="res", go_ids=None):
        """Check similarity measure"""
        list_sim_measures = ["res", "lin", "jcn", "sim_rel"]
        if sim_measure not in list_sim_measures:
            raise ValueError("'sim_measure' not in {}".format(list_sim_measures))
        if len(go_ids) != 2 and sim_measure in list_sim_measures[1:]:
            raise ValueError("'{}' just accepts two go_ids".format(sim_measure))

    def _semantic_similarity(self, go_ids=None, sim_measure="res"):
        """Check similarity measure (node-based) for given set of go_ids. Resnik's measure accepts more than two go_ids.
        The other measures (Lin, Jiang & Conrath, and sim_rel) accept just two go_ids."""
        self._check_similarity_measure(sim_measure=sim_measure, go_ids=go_ids)
        if sim_measure == "res":
            return self._res_similarity(go_ids=go_ids)
        elif sim_measure == "lin":
            return self._lin_similarity(go_ids=go_ids)
        elif sim_measure == "jcn":
            return self._jcn_similarity(go_ids=go_ids)
        else:
            return self._sim_rel(go_ids=go_ids)


# II Main Functions
class SemanticGo(SemanticGoBase):
    """Class for loading GO DB and to query GO terms"""

    def __init__(self, dict_go=None, termcounts=None):
        SemanticGoBase.__init__(self, dict_go=dict_go, termcounts=termcounts)

    # Helper methods
    @staticmethod
    def _p_value_filtering(df_enrich=None, ns="BP", p_max=0.05):
        """Go ids from df_go_enrich based on p value"""
        go_ids = list(df_enrich[(df_enrich["p_uncorrected"] < p_max) & (df_enrich["NS"] == ns)]["GO"])
        if len(go_ids) < 5:
            raise ValueError("p value filtering settings are to strict for semantic analysis")
        p_val = list(df_enrich[df_enrich["p_uncorrected"] < p_max]["p_uncorrected"])
        dict_go_id_p_val = dict(zip(go_ids, p_val))
        return dict_go_id_p_val

    def _pairwise_similarity(self, go_ids=None, sim_measure="jcn"):
        """Compute pairwise similarity for all given go_ids"""
        m_sim = np.zeros((len(go_ids), len(go_ids)))
        for i, g1 in enumerate(go_ids):
            for j, g2 in enumerate(go_ids):
                if i >= j:
                    m_sim[i, j] = 0
                else:
                    m_sim[i, j] = self._semantic_similarity(go_ids=[g1, g2], sim_measure=sim_measure)
        # Normalize similarity
        max_val = np.max(m_sim)
        m_sim = m_sim / max_val
        return m_sim

    @staticmethod
    def _merge_clusters(dict_clusters=None, go_id_keep=None, go_id_remove=None):
        """Merge clusters by removing go_idb from cluster and append add to go_ida list including all ids of cluster"""
        cluster = dict_clusters[go_id_remove]
        cluster.append(go_id_remove)
        del dict_clusters[go_id_remove]
        dict_clusters[go_id_keep].extend(cluster)
        return dict_clusters

    # 1. Revigo clustering
    def _revigo_clustering(self, sim_min=0.7, df_sim=None, dict_go_id_p_val=None, dict_clusters=None):
        """Recursive clustering algorithm to select non-redundant set of GO terms (go_ids) based on
        semantic similarity (min_similarity) for given min_similarity and significant (dict_go_id_p_val)
        In: a) min_similarity: C threshold for filtering 0.9 (large), 0.7 (medium), 0.5 (small). 0.4 (tiny)"""
        # 1. Filter similarity df for remaining go_ids
        go_ids = dict_clusters.keys()
        df_sim = df_sim[go_ids].loc[go_ids]
        # 2. Find most similar pair t1 & t2 (Highest Resnik's similarity)
        max_val = df_sim.max().max()
        # 3. Are t1 & t2 less similar than user-specified cutoff min_similarity (yes -> finish)
        if max_val < sim_min:    # Stop condition
            return dict_clusters
        else:
            # 4. Remove either t1 or t2 depending on criteria (listed in order of priority)
            go_ida = df_sim.max(axis=0).idxmax()    # col_max
            go_idb = df_sim.max(axis=1).idxmax()   # index max
            arg_keep_b = {"dict_clusters": dict_clusters, "go_id_keep": go_idb, "go_id_remove": go_ida}
            arg_keep_a = {"dict_clusters": dict_clusters, "go_id_keep": go_ida, "go_id_remove": go_idb}
            # a) One term has a very broad interpretation (frequency >5%) -> reject less specific
            ic_a = self._information_content(go_id=go_ida)
            ic_b = self._information_content(go_id=go_idb)
            if ic_a > ic_b:
                dict_clusters = self._merge_clusters(**arg_keep_b)
            elif ic_a < ic_b:
                dict_clusters = self._merge_clusters(**arg_keep_a)
            else:
                # b) One term has a less significant p-value (weaker enrichment) -> reject less significant
                if dict_go_id_p_val[go_ida] > dict_go_id_p_val[go_idb]:
                    dict_clusters = self._merge_clusters(**arg_keep_b)
                elif dict_go_id_p_val[go_ida] < dict_go_id_p_val[go_idb]:
                    dict_clusters = self._merge_clusters(**arg_keep_a)
                else:
                    # c) t1 and t2 are in parent-child relation -> reject child
                    parents_a = self.get_parents(go_id=go_ida)
                    parents_b = self.get_parents(go_id=go_idb)
                    if go_ida in list(parents_b):
                        dict_clusters = self._merge_clusters(**arg_keep_b)
                    elif go_idb in list(parents_a):
                        dict_clusters = self._merge_clusters(**arg_keep_a)
                    else:
                        # Nothing of a-c -> reject random (t1 to ensure reproducibility)
                        dict_clusters = self._merge_clusters(**arg_keep_b)
            return self._revigo_clustering(df_sim=df_sim,
                                           sim_min=sim_min,
                                           dict_go_id_p_val=dict_go_id_p_val,
                                           dict_clusters=dict_clusters)

    def revigo(self, df_enrich=None, p_max=0.05, ns="BP", sim_min=0.7, sim_measure="res"):
        """ Implementation of REVIGO algorithm (semantic clustering) for removing redundant GO terms
        Ref: Supek et al., 2011 (REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms)
            https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0021800
        In: a) df_enrich: df from enrichment
            b1) p_method: method for p value correction ("bonferroni", "sidak", "holm", "fdr")
            b2) p_max: threshold for p value
            c) ns: namespace for filtering 'BP' (bio. process) , 'MF' (mol function), 'CC' (cellular component)
            d1) sim_min: minimum of semantic similarity for REVIGO clustering (0.9 large dataset, 0.4 tiny dataset
            d2) sim_measure: similarity measure ("res", "lin", "jcn", "sim_rel")
        Out:a) dict_clusters: dict with GO clusters (keys are non-redundant set of representative GO ids)"""
        # Check input
        self._check_namespace(ns=ns)
        dict_go_id_p_val = self._p_value_filtering(df_enrich=df_enrich, p_max=p_max, ns=ns)
        go_ids = list(dict_go_id_p_val.keys())
        # REVIGO Algorithm
        # Similarity matrix for all GO ids from given GO (for normalization)
        matrix_sim = self._pairwise_similarity(go_ids=go_ids, sim_measure=sim_measure)
        # Clustering on matrix for go_ids that are associated with study terms
        df_sim = pd.DataFrame(matrix_sim, columns=go_ids, index=go_ids)
        dict_clusters = {go_id: [] for go_id in go_ids}
        dict_clusters = self._revigo_clustering(df_sim=df_sim,
                                                dict_go_id_p_val=dict_go_id_p_val,
                                                sim_min=sim_min,
                                                dict_clusters=dict_clusters)
        return dict_clusters

    # Fold enrichment
    def cluster_go_enrichment(self, df_enrich=None, dict_clusters=None):
        """Calculate fold enrichment as in DAVID (geometric mean of uncorrected p values all go terms in one cluster)"""
        all_go_ids = self._get_items_from_dict(dict_key_list=dict_clusters)
        all_study_items = self._get_items_from_col(df=df_enrich[df_enrich["GO"].isin(all_go_ids)], col="study_items")
        list_enrichment_items = []
        for go_id in dict_clusters:
            # Get go ids
            go_ids_in_cluster = dict_clusters[go_id]
            list_go_ids = go_ids_in_cluster + [go_id]
            df = df_enrich[df_enrich["GO"].isin(list_go_ids)]
            ns = df[df["GO"] == go_id]["NS"].values[0]
            depth = df[df["GO"] == go_id]["depth"].values[0]
            study_items = self._get_items_from_col(df=df)
            name = df[df["GO"] == go_id]["NAME"].values[0]
            enrichment_score = self._enrichment_score(df=df)
            # Count and percentage
            str_go_ids_in_cluster = self.list_to_str(list_go_ids)
            str_study_items = self.list_to_str(study_items)
            count_study_items = len(study_items)
            count_go_ids = len(list_go_ids)
            freq_study_items = round(count_study_items / len(all_study_items) * 100, 2)
            freq_go_ids = round(count_go_ids / len(all_go_ids) * 100, 2)
            list_enrichment_items.append([go_id, ns, name, enrichment_score, count_go_ids, count_study_items,
                                          freq_go_ids, freq_study_items, depth, str_go_ids_in_cluster, str_study_items])
        list_col = ["GO", "NS", "NAME", "enrichment_score", "GO_count", "study_count",
                    "GO_%", "study_%", "depth", "GO_cluster", "study_items"]
        df_enrich_clust = pd.DataFrame(list_enrichment_items, columns=list_col)
        df_enrich_clust.sort_values(by=["enrichment_score", "study_%", "depth"], inplace=True, ascending=False)
        df_enrich_clust.reset_index(drop=True, inplace=True)
        return df_enrich_clust

    def run(self, df_enrich=None, df_enrich_filtered=None, sim_min=0.4, n=4):
        """Run REVIGO algorithm (semantic clustering) for removing redundant GO terms and sum results up in df
        Ref: Supek et al., 2011 (REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms)
            https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0021800
        In: a1) df_enrich: df from enrichment analysis
            a2) df_enrich_filtered: p value filtered df from enrichment analysis
            b) sim_min: minimum of semantic similarity for REVIGO clustering (0.9 large dataset, 0.4 tiny dataset
            c) n: number of maximal chosen clusters for sub-ontology/namespace
        Out:a) df_sem_clust: df with GO clusters
        Note
        ----
        Resnik's semantic similarity measure (i.e. information content of nearest common parent) used
        since it is only measure accounting for more than two GO terms.
        E.g. resnik similarity for 'biological process' = 0 because it occurs in every go term as parent
        Ref: Resnik, 1995 (Using Information Content to Evaluate Semantic Similarity in a Taxonomy)
        """
        # sim_min (similarity minimum): 0.4 (tiny), 0.5 (medium), 0.7 (default), 0.9 (large)
        # na (name space): 'BP' (biological process), 'MF' (molecular function), 'CC' (cellular component)
        list_df_enrich_clust = []
        for ns in ["BP", "CC", "MF"]:
            if len(df_enrich_filtered[df_enrich_filtered["NS"] == ns]) > 3:
                dict_clusters = self.revigo(df_enrich=df_enrich,
                                            sim_min=sim_min,
                                            ns=ns,
                                            sim_measure="res")
                df_go_enrich_clust = self.cluster_go_enrichment(df_enrich=df_enrich,
                                                                dict_clusters=dict_clusters)
                if df_go_enrich_clust is not None:
                    list_df_enrich_clust.append(df_go_enrich_clust.head(n))
        if len(list_df_enrich_clust) > 0:
            df_sem_clust = pd.concat(list_df_enrich_clust)
        else:
            df_sem_clust = pd.DataFrame()
        return df_sem_clust
