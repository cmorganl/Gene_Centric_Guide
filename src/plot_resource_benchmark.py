#!/usr/bin/env python3

import os

import pandas as pd
from scipy.signal import savgol_filter
import plotly.express as px
import plotly.graph_objs as go

from evaluate_clusters import write_images_from_dict

_LABEL_MAT = {"Fasta.Length": "Query size (aa/bp)",
              "Memory.Max (kbytes)": "Maximum memory (kbytes)"}
_CATEGORIES = {"Software": ["GraftM", "DIAMOND", "TreeSAPP", "TreeSAPP-aELW", "TreeSAPP-Max(LWR)"],
               "Molecule": ["nuc", "aa"]}


def update_figure_aesthetics(fig: go.Figure) -> None:
    fig.update_xaxes(matches=None)
    fig.update_xaxes(autorange=True)
    fig.update_traces(line=dict(width=12),
                      marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)
    fig.update_layout(legend=dict(itemsizing="constant"))
    return


def ram_plot(df: pd.DataFrame, output_dir: str) -> None:
    group = ["Software", "Fasta.Length", "Molecule"]
    plot_prefix = os.path.join(output_dir, "RAM_line")
    mem_df = df.groupby(group).max(numeric_only=True)["Memory.Max (kbytes)"].reset_index()
    mem_df["Memory.Max (kbytes)"] = savgol_filter(mem_df["Memory.Max (kbytes)"],
                                                  window_length=55, polyorder=5)
    fig = px.line(mem_df,
                  x="Fasta.Length", y="Memory.Max (kbytes)",
                  color="Software", line_group="Software",
                  color_discrete_sequence=px.colors.qualitative.Set3,
                  facet_col="Molecule",
                  line_shape="spline",
                  category_orders=_CATEGORIES,
                  labels=_LABEL_MAT,
                  render_mode="svg")
    update_figure_aesthetics(fig)
    write_images_from_dict({plot_prefix: fig})
    return


def time_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "compute_time")
    fig = px.line(df,
                  x="Fasta.Length", y="Time (m)",
                  color="Software", line_group="Software",
                  facet_col="Molecule",
                  facet_row="Threads",
                  color_discrete_sequence=px.colors.qualitative.Set3,
                  category_orders=_CATEGORIES,
                  labels=_LABEL_MAT,
                  render_mode="svg")

    update_figure_aesthetics(fig)
    fig.update_traces(mode='markers+lines')
    write_images_from_dict({plot_prefix: fig})
    return


def convert_hhmmss_to_float(vec: list) -> list:
    float_vec = []
    for t in vec:
        seconds = 0.0
        scalar = 1
        groups = [float(x) for x in t.split(':')]
        while groups:
            seconds += scalar * groups.pop()
            scalar = scalar * 60
        float_vec.append(seconds)
    return float_vec


def merge_graftm_resources(data_df: pd.DataFrame) -> pd.DataFrame:
    """GraftM can only analyze one marker at a time so its resources need to be merged"""
    gm_df = data_df[data_df["Software"] == "GraftM"]
    other_df = data_df[data_df["Software"] != "GraftM"]

    sum_gm_df = gm_df.groupby(["Software", "Sample.Name", "Fasta.Length", "Molecule", "Threads"]).sum(numeric_only=True)
    max_gm_df = gm_df.groupby(["Software", "Sample.Name", "Fasta.Length", "Molecule", "Threads"]).max(numeric_only=True)
    max_gm_df["Time (s)"] = sum_gm_df["Time (s)"]

    data_df = pd.concat([gm_df, other_df])
    return data_df


def main(time_table: str, figures_dir: str, tables_dir: str) -> None:
    data_df = pd.read_csv(time_table, sep=",", header=0)
    data_df["Time (s)"] = convert_hhmmss_to_float(data_df["Time (mm:ss)"])
    data_df = merge_graftm_resources(data_df)
    data_df["Time (m)"] = round(data_df["Time (s)"] / 60, 2)

    sum_df = data_df.groupby(["Software", "Threads", "Molecule"]).sum(numeric_only=True).reset_index()
    sum_df["chars/s"] = round(sum_df["Fasta.Length"] / sum_df["Time (s)"], 0)
    sum_df.drop(labels=["Memory.Max (kbytes)", "Time (m)"], axis=1, inplace=True)
    sum_df.to_csv(os.path.join(tables_dir, "sum_runtimes.tsv"),
                  sep="\t", header=True, index=False, float_format='%.2e')

    # Remove DIAMOND samples
    data_df = data_df[data_df["Software"] != "DIAMOND"]

    ram_plot(data_df, figures_dir)

    time_plot(data_df, figures_dir)

    return


if __name__ == '__main__':
    main(time_table="../manuscript/tables/runtime_log.csv",
         figures_dir="../manuscript/figures/",
         tables_dir="../manuscript/tables/")
