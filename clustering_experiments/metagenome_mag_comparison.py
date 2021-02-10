#!/usr/bin/env python3

import re
import os
import glob

import plotly.express as px
import pandas as pd

from evaluate_clusters import ClusterExperiment


def plot_potu_stats(otu_mat: pd.DataFrame) -> None:
    palette = px.colors.cyclical.Twilight
    labs = {}
    bar_plt = px.bar(otu_mat.groupby(["Data_type"]).sum().reset_index(),
                     x="Data_type", y="Count", color="Data_type",
                     color_discrete_sequence=palette,
                     labels=labs,
                     title="Total ORFs classified between Sakinaw Lake metagenomes and MAGs")
    bar_plt.update_traces(marker_line_color='rgb(8,48,107)',
                          marker_line_width=1.5, opacity=0.6)
    bar_plt.update_layout(xaxis={'visible': False, 'showticklabels': False})
    bar_plt.show()

    bar_plt = px.bar(otu_mat.groupby(["Data_type", "RefPkg"]).nunique().reset_index(),
                     x="RefPkg", y="Count", color="Data_type",
                     color_discrete_sequence=palette, barmode="group",
                     labels=labs,
                     title="Discrete pOTUs identified in Sakinaw Lake metagenomes and MAGs")
    bar_plt.show()
    return


def compare_cluster_occupancy(dir_path):
    matrices = glob.glob(dir_path + os.sep + "**" + os.sep + "phylotu_matrix.tsv", recursive=True)
    master_mat = None
    for mat in matrices:  # type: str
        phylotu_out_dir = mat.split('final_outputs')[0]
        ce = ClusterExperiment(directory=phylotu_out_dir)
        ce.pkg_name = re.match(r"phylotu_(\w+)_.*", mat.split(os.sep)[1]).group(1)

        # Read the matrix
        otu_mat = ce.load_otu_matrix()

        # Reformat the columns for pandas.wide_to_long()
        otu_mat.columns = [re.sub('.*sak_', '', col) for col in otu_mat.columns]
        prefixed_cols = []
        for col in otu_mat.columns:
            if col == "#OTU_ID":
                prefixed_cols.append(col)
                continue
            ws = col.split('_')
            prefixed_cols.append(ws[-1] + '_'.join(ws[:-1]))
        otu_mat.columns = prefixed_cols

        # Pivot table to long-table format
        otu_mat = pd.wide_to_long(otu_mat, ["MAGs", "metaG"], j="Sample", i="#OTU_ID", suffix=r'\w+')
        otu_mat["RefPkg"] = ce.pkg_name

        if master_mat is None:
            master_mat = otu_mat
        else:
            master_mat = master_mat.append(otu_mat)

    master_mat = master_mat.reset_index().melt(id_vars=["#OTU_ID", "Sample", "RefPkg"],
                                               var_name="Data_type",
                                               value_name="Count")
    plot_potu_stats(master_mat)

    return


if __name__ == "__main__":
    compare_cluster_occupancy("Sakinaw")
