# Proteomics Tools - Collection of bioinformatic tools for protein analysis
Proteomics Tools comprises following lightweighted tools to study protein expression:
1. API Loader: Updating underlying database (e.g., KEGG IDs for pathways and diseases)
2. Gene Enrichment: Tools to perform enrichment analysis (e.g., for GO terms or KEGG IDs)
3. Mapper: A simple mapper between different bioinformatic ids (e.g., from UniProt ACC to gene name)
4. Parsers: Collection of different sequence parsers
5. Perseus Pipeline: Perform proteomic analysis as in Perseus (e.g., creating Volcano plot)
6. Retrievers: Retrieve and compare expression information for various proteins
    based on gold-standard databases (e.g., ProteomicDB or Human Protein Atlas)
7. Time Series Analysis: Pipeline to perform time series analysis (e.g., clustering) on expression data
8. Uniprot Features: Tool to include UniProt annotation information (e.g., subcellular location) 
    to pd.DataFrame
9. Potency analysis

The underlying database is available via DropBox: https://www.dropbox.com/home/Share/Proteomics_Tools