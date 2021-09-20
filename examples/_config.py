#! /usr/bin/python3
"""
Config with folder structure
"""
import os
import platform


# Helper function
def _folder_path(super_folder, folder_name):
    """Modification of separator (OS depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


# Folder
SEP = "\\" if platform.system() == "Windows" else "/"

FOLDER_PROJECT = os.path.dirname(os.path.abspath(__file__)).split('examples')[0].replace('/', SEP) + "examples" + SEP

FOLDER_DATA = _folder_path(FOLDER_PROJECT, "datasets")
FOLDER_RESULTS = _folder_path(FOLDER_PROJECT, 'results')
