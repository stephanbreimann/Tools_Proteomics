"""
This is a script for getting data folder
"""
import os
import platform


URL = "https://www.dropbox.com/home/Share/Proteomics_Tools"
# Helper methods
def _folder_path(super_folder, folder_name):
    """Modification of separator (OS depending)"""
    path = os.path.join(super_folder, folder_name + SEP)
    return path


def find_path_to_data_folder(dir_name=None):
    root, dirs, files = next(os.walk(dir_name))
    parent_root = SEP.join(root.split(SEP)[0:-1]) + SEP
    root, dirs, files = next(os.walk(parent_root))
    parent_root = SEP.join(root.split(SEP)[0:-1]) + SEP
    found_folder = False
    list_parents = [parent_root]
    path_data = None
    counter = 0
    while not found_folder and counter < 4:
        counter += 1
        if "data" in dirs:
            found_folder = True
            path_data = root + "data" + SEP
        else:
            root, dirs, files = next(os.walk(parent_root))
            parent_root = SEP.join(root.split(SEP)[0:-2]) + SEP
            list_parents.append(parent_root)
    if path_data is None:
        error = "Data folder not found in project. Create folder 'data' in main project level" \
                f" and save data downloaded from: {URL}"
        raise ValueError(error)
    return path_data


# Folder
SEP = "\\" if platform.system() == "Windows" else "/"
FOLDER_DATA = find_path_to_data_folder(dir_name=os.getcwd())

