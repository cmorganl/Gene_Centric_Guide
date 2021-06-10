#!/usr/bin/env python3

import os

import pandas as pd
import plotly.express as px

from evaluate_clusters import write_images_from_dict

_LABEL_MAT = {"Tax.dist": "Taxonomic distance",
              "MCC": "Matthews correlation coefficient"}


def load_mcc_tables(analysis_dir: str) -> pd.DataFrame:
    tbl_name = "prokaryotic_e5.proteomes_MCC_table.tsv"
    output_dirs = {"MCC_graftm_0.13.0_prok": "GraftM",
                   "MCC_treesapp_0.11.1_aelw_prok": "TreeSAPP-aELW",
                   "MCC_treesapp_0.11.1_max-lwr_prok": "TreeSAPP-Max(LWR)"}
    data_tables = []
    for dir_name, spl in output_dirs.items():
        tbl = pd.read_csv(os.path.join(analysis_dir, dir_name, tbl_name),
                          sep="\t", header=0)
        tbl["Method"] = spl
        data_tables.append(tbl)
    return pd.concat(data_tables)


def mcc_line_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "MCC_line")
    line_plt = px.line(df,
                       x="Tax.dist", y="MCC",
                       color="Method", line_group="Method",
                       color_discrete_sequence=px.colors.qualitative.Set3,
                       labels=_LABEL_MAT,
                       range_y=[0.1, 0.95],
                       title="TreeSAPP and GraftM classification performances<br>using EggNOG prokaryotic sequences")
    line_plt.update_traces(mode='markers+lines')
    line_plt.update_xaxes(autorange=True)
    line_plt.update_layout(legend=dict(itemsizing="constant"))
    line_plt.update_traces(line=dict(width=12),
                           marker=dict(size=12,
                                       line=dict(width=2,
                                                 color='DarkSlateGrey')))
    write_images_from_dict({plot_prefix: line_plt})
    return


def evaluate_taxonomic_summary_methods(analysis_dir: str, figures_dir: str, tables_dir: str) -> None:

    data_df = load_mcc_tables(analysis_dir)

    mcc_line_plot(data_df, figures_dir)

    data_df.to_csv(path_or_buf=os.path.join(tables_dir, "MCC_data.tsv"), sep="\t", index=False, header=True)

    return


if __name__ == '__main__':
    evaluate_taxonomic_summary_methods(analysis_dir="../tax_summary/",
                                       figures_dir="../manuscript/figures/",
                                       tables_dir="../manuscript/tables/")
