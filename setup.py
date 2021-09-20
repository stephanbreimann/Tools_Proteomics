#!/usr/bin/env python

from setuptools import setup

setup(
    name="proteomics_tools",
    version='1.0',
    description='Proteomics analysis tool kit',
    author='Stephan Breimann',
    author_email='stephanbreimann@yahoo.de',
    url=None,
    packages=['proteomics_tools',
              'proteomics_tools.api_loaders',
              'proteomics_tools.gene_enrichment',
              'proteomics_tools.mapper',
              'proteomics_tools.parsers',
              'proteomics_tools.perseuspy',
              'proteomics_tools.potency_analysis',
              'proteomics_tools.potency_analysis',
              'proteomics_tools.retrievers',
              'proteomics_tools.tsc',
              'proteomics_tools.uniprot_features'],
    include_package_data=True,
    package_data={"": ["*.xlsx"]}
)