"""
This is a script for retrieving gene ontology
References
----------
Huntley RP, Sawford T, Mutowo-Meullenet P, Shypitsyna A, Bonilla C, Martin MJ, O’Donovan C
The GOA database: Gene Ontology annotation updates for 2015. Nucleic Acids Res. 2015 Jan; 43:D1057-63
"""
import gzip
import json
import os
from ftplib import FTP
from urllib.request import urlopen
import pandas as pd
import wget
from Bio.UniProt import GOA as GOA
from goatools import obo_parser

import proteomics_tools.gene_enrichment._utils as ut


# II Main Functions
class LoadGo:
    """Class for loading GO DBs"""

    def __init__(self, folder_data=None):
        if folder_data is None:
            folder_data = ut.FOLDER_DATA + "GO_KEGG/"
        self.folder_data = folder_data

    # Helper methods
    @staticmethod
    def _check_organism(organism):
        """Check if organism in list of allowed organisms from EMBL"""
        list_dirs = ["uniprot", "human", "mouse", "rat", "arabidopsis", "zebrafish",
                     "chicken", "cow", "dog", "pig", "fly", "worm", "yeast", "pdb", "proteomes"]
        if organism.lower() not in list_dirs:
            raise ValueError("{} should be in {}".format(organism, list_dirs))

    @staticmethod
    def _check_gaf_output(out):
        """Check if output format is dict or df"""
        list_out = ["dict", "df"]
        if out not in list_out:
            raise ValueError("{} not in {}".format(out, list_out))

    # Loading methods
    def basic(self, dict_out=True):
        """Load GO DB"""
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        # Check if the file exists already
        if not os.path.isfile(self.folder_data + '/go-basic.obo'):
            go_obo = wget.download(go_obo_url, self.folder_data + '/go-basic.obo')
        else:
            go_obo = self.folder_data + '/go-basic.obo'
        if dict_out:
            dict_go = obo_parser.GODag(go_obo)
            return dict_go
        else:
            return go_obo

    # Get Go Term from QuickGo
    @staticmethod
    def quickgo(go_id=None, just_results=True):
        """This function retrieves the definition of a given Gene Ontology term, sing EMBL-EBI's QuickGO browser.
        In: a) go_id - a valid Gene Ontology ID, e.g. GO:0048527
        Out:a) OBO-XML dictionary"""
        quickgo_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + go_id
        ret = urlopen(quickgo_url)
        # Check the response
        if ret.getcode() == 200:
            go_term = json.loads(ret.read())   # OBO-XML dictionary
            if just_results:
                return go_term["results"][0]
            else:
                return go_term
        else:
            raise ValueError("Couldn't receive information from QuickGO. Check GO ID and try again.")

    # Download GAF (Gene Association File: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/)
    @staticmethod
    def _gaf_to_dict_go(go_organism_gaf=None):
        """Unzip go organism gaf and save as dictionary for ids to entries"""
        # File is a gunzip file, so we need to open it in this way
        with gzip.open(go_organism_gaf, 'rt') as go_organism_fp:
            dict_go = {}    # Up_id|GO_id to entry
            # Iterate on each function using Bio.UniProt.GOA library.
            for entry in GOA.gafiterator(go_organism_fp):
                uniprot_id = entry['DB_Object_ID']
                go_id = entry["GO_ID"]
                dict_go["{}|{}".format(uniprot_id, go_id)] = entry  # GO entry
        return dict_go

    def gaf_for_organism(self, out_dict="dict", organism="human"):
        """Download current file from official GOA (https://www.ebi.ac.uk/GOA/downloads) by EMBL-EBI
        In: a) organism
        Out:a) dict_id_entry: dict for each id (in gaf) with entry information
        References
        ----------
        Huntley RP, Sawford T, Mutowo-Meullenet P, Shypitsyna A, Bonilla C, Martin MJ, O’Donovan C
        The GOA database: Gene Ontology annotation updates for 2015. Nucleic Acids Res. 2015 Jan; 43:D1057-63
        """
        self._check_organism(organism)
        self._check_gaf_output(out_dict)
        go_organism_uri = '/pub/databases/GO/goa/{}/goa_{}.gaf.gz'.format(organism.upper(), organism.lower())
        go_organism_file = go_organism_uri.split('/')[-1]
        # Check if the file exists already
        go_organism_gaf = os.path.join(self.folder_data, go_organism_file)  # GAF: GO Annotation Format
        if not os.path.isfile(go_organism_gaf):
            # Login to FTP server
            ebi_ftp = FTP('ftp.ebi.ac.uk')
            ebi_ftp.login()  # Logs in anonymously
            # Download
            with open(go_organism_gaf, 'wb') as arab_fp:
                ebi_ftp.retrbinary('RETR {}'.format(go_organism_uri), arab_fp.write)
            # Logout from FTP server
            ebi_ftp.quit()
        dict_id_entry = self._gaf_to_dict_go(go_organism_gaf=go_organism_gaf)
        if out_dict == "dict":
            return dict_id_entry
        else:
            df_go = pd.DataFrame.from_dict(dict_id_entry, orient="index")
            df_go = df_go.applymap(lambda x: x[0] if isinstance(x, list) else x)  # Unzip lists
            return df_go
