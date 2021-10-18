#! /usr/bin/python3
"""
Config with folder structure
"""
import os
import platform
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


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


def plot_settings(fig_format="pdf", verbose=False, grid=False, grid_axis="y"):
    """General plot settings"""
    if verbose:
        print(plt.rcParams.keys)  # Print all plot settings that can be modified in general
    sns.set_context("talk", font_scale=0.9)
    # plt.style.use("Solarize_Light2") # https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    # Font settings https://matplotlib.org/3.1.1/tutorials/text/text_props.html
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Arial"
    plt.rcParams["axes.labelweight"] = "bold"
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["axes.labelsize"] = 17  # 13.5
    plt.rcParams["axes.titlesize"] = 16.5  # 15
    if fig_format == "pdf":
        mpl.rcParams['pdf.fonttype'] = 42
    elif "svg" in fig_format:
        mpl.rcParams['svg.fonttype'] = 'none'
    font = {'family': 'Arial',
            "weight": "bold"}
    mpl.rc('font', **font)
    # Error bars
    plt.rcParams["errorbar.capsize"] = 10  # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.errorbar.html
    # Grid
    plt.rcParams["axes.grid.axis"] = grid_axis  # 'y', 'x', 'both'
    plt.rcParams["axes.grid"] = grid
    # Legend
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.fontsize"] = "medium"  # "x-small"
    plt.rcParams[
        "legend.loc"] = 'upper right'  # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.legend.html
