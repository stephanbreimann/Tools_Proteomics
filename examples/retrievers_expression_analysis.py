"""
This is a script for ...
"""
import time
import pandas as pd
import matplotlib.pyplot as plt

import examples._config as conf
from proteomics_tools import retrievers as rt

# Settings
pd.set_option('expand_frame_repr', False)  # Single line print for pd.Dataframe


# I Helper Functions


# II Main Functions
def cell_surface_atlas_caller():
    """"""
    organism = "hsa"
    csa = rt.CellSurfaceAtlas(organism=organism)
    df_cell_line_info = csa.cell_line_info(primary=None, filter_df=False)
    # Get list of genes
    list_genes = ["TSPAN14", "TM2D3", "TM2D1", "TM2D2"][0:3]
    #list_acc = ["Q8NG11", "Q9BRN9", "Q9BX73", "Q9BX74"]
    df_cell_lines = csa.expression_surface(genes_target=list_genes,
                                           genes_filter=list_genes[1],
                                           acc=True)
    print(df_cell_lines)


def dict_caller():
    """"""
    pr = rt.ProteomicsReferences()
    _dict_sharma = pr.dict_sharma(filter_isoforms=True,
                                  split_ids=True,
                                  log2_fc_min=0.1,
                                  p_max=0.05)
    _dict_tueshaus = pr.dict_tueshaus(key_genes=False)
    list_acc = ["P12023", "Q7TSK2", "P20917"]     # APP, SEZ6, MAG (Mouse)  ???
    pr.expression_cell_types(acc_target=list_acc, dict_ref=_dict_sharma)


def proteomics_db_caller():
    """"""
    genes_target = ["APP", "TSPAN14", "TM2D3", "TM2D2", "TM2D1"]
    pdb = rt.ProteomicsDB()
    df = pdb.expression_proteins(genes_target=genes_target, genes_filter=genes_target[1],
                                 min_exp=5, biological_source="Cell_Line")
    print(set(df.index))
    pdb.heatmap(df=df)
    #plt.tick_params(rotation=25, axis="y")
    plt.show()
    file_out = "protein_expression_{}".format("_".join(genes_target))
    #df.to_excel(conf.folder_results + file_out + "xlsx")
    #plt.savefig(conf.folder_results + file_out + "png", bbox_inches="tight")


def covid_co_expression():
    """"""
    list_genes = ["ACE2", "ADAM10", "ADAM17", "TMPRSS2",  "AVPR1B", "AGTR1"]
    co_expressed_genes = ["ADAM10", "ADAM17"] # "AVPR1B", "AGTR1"]
    folder_out = conf.FOLDER_RESULTS + "HPA_CoEpxression_" + "_".join(co_expressed_genes) + "/"

    hpa = rt.HumanProteinAtlas()
    kwargs = {"list_genes": list_genes,
              "co_expressed_genes": co_expressed_genes,
              "heatmap": True}
    files_excluded = ['rna_mouse_brain_hpa.tsv', 'rna_celline.tsv']
    hpa.tissue_type_expression(**kwargs, level="ssRNA", min_th=10)
    hpa.tissue_type_expression(**kwargs, level="Protein")
    hpa.all_rna_expression(list_genes=list_genes,
                           co_expressed_genes=co_expressed_genes,
                           min_rna=0.5,
                           max_n=15,
                           files_excluded=files_excluded)


def xao_cell_line_co_expression():
    """"""
    list_genes = ["APP", "TM2D3", "TM2D2", "TM2D1", "MS4A4A", "ATP8B4"]
    # list_genes = [x.lower().capitalize() for x in list_genes]
    file = "rna_celline"
    # Adust figure size
    hpa = rt.HumanProteinAtlas()
    df = hpa.rna_expression(file=file, list_genes=list_genes, co_expressed_genes=["MS4A4A"],
                            filter_col="NX", min_rna=5, min_rna_step=0.5, max_n=30, heatmap=True, annot=True)


def kathrin_cell_line_tspan14():
    """"""
    genes_target = ["APP", "TSPAN14", "TM2D3", "TM2D2", "TM2D1"]
    pdb = rt.ProteomicsDB()
    df = pdb.expression_proteins(genes_target=genes_target, genes_filter=genes_target[1],
                                 min_exp=5, biological_source="Cell_Line")
    pdb.heatmap(df=df)
    # plt.tick_params(rotation=25, axis="y")
    plt.show()
    file_out = "protein_expression_{}".format("_".join(genes_target))
    hpa = rt.HumanProteinAtlas()
    file = "rna_celline"
    hpa.all_rna_expression(list_genes=genes_target,
                           co_expressed_genes=genes_target,
                           min_rna=0.5,
                           max_n=15)
    #df = hpa.rna_expression(file=file, list_genes=genes_target[1], co_expressed_genes=genes_target,
    #                        filter_col="NX", min_rna=5, min_rna_step=0.5, max_n=30, heatmap=True, annot=True)


# III Test/Caller Functions


# IV Main
def main():
    t0 = time.time()
    covid_co_expression()
    #cell_surface_atlas_caller()
    #dict_caller()
    #proteomics_db_caller()
    #kathrin_cell_line_tspan14()
    t1 = time.time()
    print("Time:", t1 - t0)


if __name__ == "__main__":
    main()


