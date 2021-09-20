"""
This is a script with a mapping class for ids of bioservices
A) Mapper for uniprot to GO term
Q:  a) https://bioservices.readthedocs.io/en/master/
    b) https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html
"""
import pandas as pd
import time
from bioservices import UniProt
import warnings
import numpy as np

import proteomics_tools._utils as ut


# Settings
# Valid IDs for mapping dictionary (otherwise direct retrieving via rest API of bio services)
LIST_FR_TO = ["ACC", "ID", "KEGG_ID", "GENECARDS_ID", "Gene_Name"]

# Valid IDs for bioservices REST API (used in GO_KEGG/Gene Ontology)
# https://www.genome.jp/kegg/catalog/org_list.html
MAP_BIOSERVICE = {"UniProtKB AC/ID": "ACC+ID",
                  "UniProtKB AC": "ACC",   # Modified
                  "UniProtKB ID": "ID",    # Modified
                  "UniParc": "UPARC",
                  "UniRef50": "NF50",
                  "UniRef90": "NF90",
                  "UniRef100": "NF100",
                  "EMBL/GenBank/DDBJ": "EMBL_ID",
                  "EMBL/GenBank/DDBJ CDS": "EMBL",
                  "PIR": "PIR",
                  "UniGene": "UNIGENE_ID",
                  "Entrez Gene (GeneID)": "P_ENTREZGENEID",
                  "GI number*": "P_GI",
                  "IPI": "P_IPI",
                  "RefSeq Protein": "P_REFSEQ_AC",
                  "RefSeq Nucleotide": "REFSEQ_NT_ID",
                  "PDB": "PDB_ID",
                  "DisProt": "DISPROT_ID",
                  "HSSP": "HSSP_ID",
                  "DIP": "DIP_ID",
                  "MINT": "MINT_ID",
                  "Allergome": "ALLERGOME_ID",
                  "MEROPS": "MEROPS_ID",
                  "mycoCLAP": "MYCOCLAP_ID",
                  "PeroxiBase": "PEROXIBASE_ID",
                  "PptaseDB": "PPTASEDB_ID",
                  "REBASE": "REBASE_ID",
                  "TCDB": "TCDB_ID",
                  "PhosSite": "PHOSSITE_ID",
                  "DMDM": "DMDM_ID",
                  "Aarhus/Ghent-2DPAGE": "AARHUS_GHENT_2DPAGE_ID",
                  "World-2DPAGE": "WORLD_2DPAGE_ID",
                  "DNASU": "DNASU_ID",
                  "Ensembl": "ENSEMBL_ID",
                  "Ensembl Protein": "ENSEMBL_PRO_ID",
                  "Ensembl Transcript": "ENSEMBL_TRS_ID",
                  "Ensembl Genomes": "ENSEMBLGENOME_ID",
                  "Ensembl Genomes Protein": "ENSEMBLGENOME_PRO_ID",
                  "Ensembl Genomes Transcript": "ENSEMBLGENOME_TRS_ID",
                  "GeneID": "P_ENTREZGENEID",
                  "GenomeReviews": "GENOMEREVIEWS_ID",
                  "GO_KEGG": "KEGG_ID",
                  "PATRIC": "PATRIC_ID",
                  "UCSC": "UCSC_ID",
                  "VectorBase": "VECTORBASE_ID",
                  "AGD": "AGD_ID",
                  "ArachnoServer": "ARACHNOSERVER_ID",
                  "CGD": "CGD",
                  "ConoServer": "CONOSERVER_ID",
                  "CYGD": "CYGD_ID",
                  "dictyBase": "DICTYBASE_ID",
                  "EchoBASE": "ECHOBASE_ID",
                  "EcoGene": "ECOGENE_ID",
                  "euHCVdb": "EUHCVDB_ID",
                  "EuPathDB": "EUPATHDB_ID",
                  "FlyBase": "FLYBASE_ID",
                  "GeneCards": "GENECARDS_ID",
                  "GeneFarm": "GENEFARM_ID",
                  "GenoList": "GENOLIST_ID",
                  "H-InvDB": "H_INVDB_ID",
                  "HGNC": "HGNC_ID",
                  "HPA": "HPA_ID",
                  "LegioList": "LEGIOLIST_ID",
                  "Leproma": "LEPROMA_ID",
                  "MaizeGDB": "MAIZEGDB_ID",
                  "MIM": "MIM_ID",
                  "MGI": "MGI_ID",
                  "neXtProt": "NEXTPROT_ID",
                  "Orphanet": "ORPHANET_ID",
                  "PharmGKB": "PHARMGKB_ID",
                  "PomBase": "POMBASE_ID",
                  "PseudoCAP": "PSEUDOCAP_ID",
                  "RGD": "RGD_ID",
                  "SGD": "SGD_ID",
                  "TAIR": "TAIR_ID",
                  "TubercuList": "TUBERCULIST_ID",
                  "WormBase": "WORMBASE_ID",
                  "WormBase Transcript": "WORMBASE_TRS_ID",
                  "WormBase Protein": "WORMBASE_PRO_ID",
                  "Xenbase": "XENBASE_ID",
                  "ZFIN": "ZFIN_ID",
                  "eggNOG": "EGGNOG_ID",
                  "GeneTree": "GENETREE_ID",
                  "HOGENOM": "HOGENOM_ID",
                  "HOVERGEN": "HOVERGEN_ID",
                  "KO": "KO_ID",
                  "OMA": "OMA_ID",
                  "OrthoDB": "ORTHODB_ID",
                  "ProtClustDB": "PROTCLUSTDB_ID",
                  "BioCyc": "BIOCYC_ID",
                  "Reactome": "REACTOME_ID",
                  "UniPathWay": "UNIPATHWAY_ID",
                  "CleanEx": "CLEANEX_ID",
                  "GermOnline": "GERMONLINE_ID",
                  "ChEMBL": "CHEMBL_ID",
                  "ChiTaRS": "CHITARS_ID",
                  "DrugBank": "DRUGBANK_ID",
                  "GenomeRNAi": "GENOMERNAI_ID",
                  "NextBio": "NEXTBIO_ID"}


# I Helper Functions
def check_mapper_arg(fr="ACC", to="ID"):
    """Check if mapper arguments are valid"""
    valid_ids = list(MAP_BIOSERVICE.values())
    if fr not in valid_ids or to not in valid_ids:
        raise ValueError(f"'fr' ({fr}) and 'to' ({to}) must be on of following: {valid_ids}")
    if fr not in LIST_FR_TO or to not in LIST_FR_TO:
        warn = "'fr' ({}) and/or 'to' ({}) are not main ids ({}).\n\tEntries will be retrieved " \
               "via BioServices REST API, which can take some time.".format(fr, to, LIST_FR_TO)
        warnings.warn(warn)
        time.sleep(1)


def check_organism(organism="hsa"):
    """Check if organism valid and set to right term"""
    # https://www.genome.jp/kegg/catalog/org_list.html
    list_organism = ["hsa", "mmu", "tah"]   # Human, Mouse, Arabidopsis Thaliana
    dict_org = {"mouse": ["mmu", "mus", "mouse"],
                "human": ["hsa", "human"],
                "arabidopsis": ["tah"]}
    hit = False
    if organism in list_organism:
        return organism
    for key in dict_org:
        if organism.lower() in dict_org[key]:
            organism = dict_org[key][0]
            hit = True
    if not hit:
        raise ValueError("'organism' must be from {}".format(list_organism))
    return organism


# II Main Functions
class Mapper:
    """Mapping class between ids from list_fr_to"""
    def __init__(self, organism="hsa", folder_up=None, verbose=False):
        if folder_up is None:
            folder_up = ut.FOLDER_DATA + "UniProt" + ut.SEP
        self.verbose = verbose
        self.folder_up = folder_up
        self.organism = check_organism(organism=organism)
        self.up = UniProt(verbose=False)

    def mapping_dict(self, fr="ACC", to="ID", reviewed=True):
        """Create dictionary for mapping"""
        check_mapper_arg(fr=fr, to=to)
        if reviewed:
            df = pd.read_csv(self.folder_up + "mapping_db_{}_reviewed.text".format(self.organism), sep="\t")
        else:
            df = pd.read_csv(self.folder_up + "mapping_db_{}.text".format(self.organism), sep="\t")
        dict_fr_to = {}
        # Check if queries in df
        df.sort_values(by="Status", inplace=True)
        for i, row in df.iterrows():
            key = str(row[fr])
            val = str(row[to])
            if ";" in val:
                val = val.split(";")[0]
            if key != "nan" and key not in dict_fr_to.keys():
                if fr == "ID":
                    key = key.split("_")[0]
                if ";" in key:
                    for k in key.split(";"):
                        dict_fr_to[k] = val
                else:
                    dict_fr_to[key] = val
        return dict_fr_to

    # Mapping via bioservices (check tutorial for dict of ids)
    def _fr_to_str(self, fr="ACC", to="ID", query_str=None, mapping_dict=None):
        """Mapper function for query string"""
        out = None
        if mapping_dict is not None:
            out = mapping_dict.get(query_str, None)
        if out is None:
            try:
                out = self.up.mapping(fr=fr, to=to, query=query_str)[query_str][0].split("_")[0]
            except KeyError:
                pass
        return out

    def _fr_to_list(self, fr="ACC", to="ID", query_list=None, mapping_dict=None):
        """Mapper function for query list"""
        list_to = []
        if self.verbose:
            print(f"Mapper retrieves {to} (n={len(query_list)}) via BioServices REST API")
        for i, q in enumerate(query_list):
            if self.verbose:
                progress = min(np.round(i/ len(query_list) * 100, 2), 100)
                progress_bar = "#" * int(progress/2) + " " * (50 - int(progress/2))
                print(f"\r   |{progress_bar}| {progress}%", end="")
            list_to.append(self._fr_to_str(fr=fr, to=to, query_str=q, mapping_dict=mapping_dict))
        if self.verbose:
            progress_bar = "#" * 50
            print(f"\r   |{progress_bar}| 100.00%")
        return list_to

    def fr_to(self, fr="ACC", to="ID", query=None):
        """Mapping function for query str or query list"""
        check_mapper_arg(fr=fr, to=to)
        kwargs = {"fr": fr, "to": to}
        if fr in LIST_FR_TO and to in LIST_FR_TO:
            mapping_dict = self.mapping_dict(fr=fr, to=to, reviewed=False)
            kwargs["mapping_dict"] = mapping_dict
        if type(query) is list:
            return self._fr_to_list(**kwargs, query_list=query)
        else:
            return self._fr_to_str(**kwargs, query_str=query)

