"""
This is a script for ...
"""
import time
import pandas as pd

import examples._config as conf
from proteomics_tools.mapper import Mapper

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def map_ensemble_to_uniprot():
    """"""
    folder_in = conf.FOLDER_DATA + "HumanProteinAtlas" + conf.SEP
    file = "rna_celline.tsv"
    df = pd.read_csv(folder_in + file, sep="\t")
    mp = Mapper(organism="mmu", verbose=True)
    mp.mapping_dict()
    list_ensem = list(set(df["Gene"]))
    list_acc = mp.fr_to(fr="ENSEMBL_ID", to="ACC", query=list_ensem)
    print(list_acc)

# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    map_ensemble_to_uniprot()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()
