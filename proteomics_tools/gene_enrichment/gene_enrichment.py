"""
This is a script for merging all GO classes into one class
"""
import pandas as pd

import proteomics_tools.gene_enrichment._utils as ut
from proteomics_tools.gene_enrichment.term_count import TermCounts
from proteomics_tools.gene_enrichment.go_load import LoadGo
from proteomics_tools.gene_enrichment.go_enrichment import EnrichmentGo
from proteomics_tools.gene_enrichment.go_kegg import KeggGo, KeggLoader
from proteomics_tools.gene_enrichment.go_semantic import SemanticGo
from proteomics_tools.gene_enrichment.go_plot import PlotGo

# I Helper Functions


# II Main Functions
class GeneEnrichment:
    """Class for performing gene enrichment analysis as in DAVID and REVIGO using"""
    def __init__(self, organism="mouse", verbose=False):
        self.organism = organism
        self.verbose = verbose

    def load_go(self):
        """Download current file from official GOA (https://www.ebi.ac.uk/GOA/downloads) by EMBL-EBI
        Out:a) dict_go: Basic GO ontology (i.e., GO terms presented in major organisms like mouse or human)
            b) dict_org_go: dict of GO ontology terms for selected organism
            c) df_org_go: df of GO ontology terms for selected organism
        References
        ----------
        Huntley RP, Sawford T, Mutowo-Meullenet P, Shypitsyna A, Bonilla C, Martin MJ, O’Donovan C
        The GOA database: Gene Ontology annotation updates for 2015. Nucleic Acids Res. 2015 Jan; 43:D1057-63
        """
        lgo = LoadGo()
        dict_go = lgo.basic()
        df_org_go = lgo.gaf_for_organism(out_dict="df", organism=self.organism)
        dict_org_go = lgo.gaf_for_organism(out_dict="dict", organism=self.organism)
        return dict_go, dict_org_go, df_org_go

    def load_kegg(self, reviewed=True):
        """Preload all KEGG """
        go_organism = ut.get_organism(organism=self.organism, out="GO")
        kl = KeggLoader(organism=go_organism, reviewed=reviewed)
        kegg_db = kl.run()
        return kegg_db

    @staticmethod
    def filter_df(df=None, p_max=0.05, p_method=None, n=None):
        """Filter df based on p value or total number of highest ranked entries"""
        df_filtered = ut.Utils().filter_df(df=df, p_max=p_max, p_method=p_method, n=n)
        return df_filtered

    # Analysis pipeline
    def go_enrichment(self, study=None, pop=None, dict_go=None, df_org_go=None, p_method=None):
        """Enrichment analysis, where e (enriched) and p (purified) means  significantly higher resp. lower.
        GO terms enrichment is performed using the python package goatools. Additionally, the fold enrichment
        is computed as in DAVID.
        In: a1) study: list with uniprot ids for study set
            a2) pop: list with uniprot ids for population (background)
            b1) dict_go: Basic GO ontology (i.e., GO terms presented in major organisms like mouse or human)
            c) df_org_go: df of GO ontology terms for selected organism
            d) p_method: str or list with str for p value corrections ["bonferroni", "sidak", "holm", "fdr"]
        Out:a1) study: filtered study list (ids that have no associated GO term are removed)
            a2) pop: filtered pop list (ids that have no associated GO term are removed)
            b) df_enrich: result of enrichment analysis
        References
        ----------
        Klopfenstein DV, Zhang L, Pedersen BS, ... Tang H GOATOOLS:
            A Python library for Gene Ontology analyses Scientific reports, 2018
        DAVID: Nature Protocols 2009; 4(1):44 & Nucleic Acids Res. 2009;37(1):1
        """
        # TODO include n of genes associated with GO term, % in study, & % in background
        ego = EnrichmentGo(dict_go=dict_go, df_org_go=df_org_go, verbose=self.verbose)
        study, pop = ego.filter_study_pop(study=study, pop=pop)
        df_enrich = ego.run(study=study,
                            pop=pop,
                            p_method=p_method)
        return study, pop, df_enrich

    @staticmethod
    def semantic_clustering(df_enrich_filtered=None, df_enrich=None, dict_go=None, dict_org_go=None, n=5, sim_min=0.4):
        """Run REVIGO algorithm (semantic clustering) for removing redundant GO terms and sum results up in df
        In: a1) df_enrich: df from enrichment analysis
            a2) df_enrich_filtered: p value filtered df from enrichment analysis
            b) sim_min: minimum of semantic similarity for REVIGO clustering (0.9 large dataset, 0.4 tiny dataset)
            c) n: number of maximal chosen clusters for sub-ontology/namespace
        Out:a) df_sem_clust: df with GO clusters
        Notes
        -----
        Resnik's semantic similarity measure (i.e. information content of nearest common parent) used
        since it is only measure accounting for more than two GO terms.
        E.g. resnik similarity for 'biological process' = 0 because it occurs in every go term as parent
        References
        ----------
        Supek et al., 2011 (REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms)
            https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0021800
        Resnik, 1995 (Using Information Content to Evaluate Semantic Similarity in a Taxonomy)
        """
        termcounts = TermCounts(dict_go=dict_go,
                                dict_org_go=dict_org_go)
        sgo = SemanticGo(dict_go=dict_go,
                         termcounts=termcounts)
        df_sem_clust = sgo.run(df_enrich_filtered=df_enrich_filtered,
                               df_enrich=df_enrich,
                               n=n, sim_min=sim_min)
        return df_sem_clust

    def kegg_analysis(self, df_enrich_filtered=None, study=None, pop=None, p_method=None, dict_kegg=None):
        """Complete KEGG enrichment analysis using bioservice package for retrieving information from the KEGG
        database via REST API for pathways and diseases.
        In: a) study: list of study items (UniProt accession numbers)
            b) pop: list of population items (UniProt accession numbers
            c) p_method: method for p-value correction ["bonferroni", "sidak", "holm", "fdr", None (default value)]
            d) df_enrich_filtered: df with enriched proteins
        Out:a) df_pathway: df with pathways for which study items are enriched
            b) df_disease: df with diseases for which study items are enriched
            c) df_study: df with study items and all associated information from KEGG
        References
        ----------
        Cokelaer et al. BioServices: a common Python package to access biological Web Services
            programmatically Bioinformatics (2013) 29 (24): 3241-3242
        """
        go_organism = ut.get_organism(organism=self.organism, out="GO")
        kgo = KeggGo(organism=go_organism, dict_kegg=dict_kegg, verbose=self.verbose)
        df_pathway, df_disease, df_study = kgo.run(study=study,
                                                   pop=pop,
                                                   p_method=p_method,
                                                   df_enrich_filtered=df_enrich_filtered.copy())
        return df_pathway, df_disease, df_study

    # Plotting
    def plot_enrichment_go(self, df_enrich_filtered=None, title="Fold Enrichment GO Terms",
                           split_name=True, p_col="p_uncorrected", n=5):
        """Plot GP Enrichment
        In: a) df_enrich_filtered: df with enriched proteins
            b) title: plot title
            c) split_name: boolean indicating whether feature names can be split in two rows
            d) p_col: p value column indicating the significance of GO terms
            e) n: number of entries per names space (MF, CC, BP)
        """
        pgo = PlotGo(verbose=self.verbose)
        df_enrich_filtered = self.filter_df(df=df_enrich_filtered.copy(), n=n)
        pgo.enrichment_plot(df=df_enrich_filtered, hue="NS", title=title, split_name=split_name, p_col=p_col)

    def plot_revigo_clustering(self, df_sem_clust=None, title="Enrichment Score REVIGO Clustering", split_name=True):
        """Plot enrichment score from REVIGO clustering"""
        pgo = PlotGo(verbose=self.verbose)
        pgo.enrichment_plot(df=df_sem_clust.copy(), x="enrichment_score", title=title, split_name=split_name)

    def plot_enrichment_kegg(self, df_pathway=None, df_disease=None, title="Fold Enrichment KEGG Analysis",
                             split_name=True, n=7, p_col="p_uncorrected"):
        """Plot fold enrichment from KEGG Analysis"""
        df_p = df_pathway.copy()
        df_d = df_disease.copy()
        pgo = PlotGo(verbose=self.verbose)
        df_p.sort_values(by=["p_uncorrected", "NAME"], inplace=True)
        df_d.sort_values(by=["p_uncorrected", "NAME"], inplace=True)
        df_p["NS"] = "Pathway"
        df_p.drop("PATHWAY", inplace=True, axis=1)
        df_d["NS"] = "Disease"
        df_d.drop("DISEASE", inplace=True, axis=1)
        df = pd.concat([df_p.head(n), df_d.head(n)], axis=0)
        pgo.enrichment_plot(df=df, title=title, p_col=p_col, split_name=split_name)
