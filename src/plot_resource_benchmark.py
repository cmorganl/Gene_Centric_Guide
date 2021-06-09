#!/usr/bin/env python3

import os

import pandas as pd
import plotly.express as px

from evaluate_clusters import write_images_from_dict

_LABEL_MAT = {"Fasta.Length": "Query size (aa/bp)",
              "Memory.Max (kbytes)": "Maximum memory (kbytes)"}
_CATEGORIES = {"Software": ["GraftM", "DIAMOND", "TreeSAPP", "TreeSAPP-aELW", "TreeSAPP-Max(LWR)"]}


def ram_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "RAM_line")
    fig = px.line(df,
                  x="Fasta.Length", y="Memory.Max (kbytes)",
                  color="Software", line_group="Software",
                  color_discrete_sequence=px.colors.qualitative.Set3,
                  facet_col="Molecule",
                  category_orders=_CATEGORIES,
                  labels=_LABEL_MAT)
    fig.update_traces(mode='markers+lines')
    fig.update_xaxes(autorange=True)
    fig.update_traces(line=dict(width=4),
                      marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)
    write_images_from_dict({plot_prefix: fig})
    return


def time_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "compute_time")
    fig = px.line(df,
                  x="Fasta.Length", y="Time (s)",
                  color="Software", line_group="Software",
                  facet_col="Molecule",
                  color_discrete_sequence=px.colors.qualitative.Set3,
                  category_orders=_CATEGORIES,
                  labels=_LABEL_MAT)
    fig.update_traces(mode='markers+lines')
    fig.update_xaxes(autorange=True)
    fig.update_traces(line=dict(width=4),
                      marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)
    write_images_from_dict({plot_prefix: fig})
    return


def parallel_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "parallelization")
    fig = px.line(df,
                  x="Threads", y="Time (s)",
                  color="Software", line_group="Software",
                  facet_col="Molecule",
                  color_discrete_sequence=px.colors.qualitative.Set3,
                  category_orders=_CATEGORIES,
                  labels=_LABEL_MAT)
    fig.update_traces(mode='markers+lines')
    fig.update_xaxes(autorange=True)
    fig.update_traces(line=dict(width=4),
                      marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)
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


def main(time_table: str, figures_dir: str) -> None:
    data_df = pd.read_csv(time_table, sep=",", header=0)
    data_df["Time (s)"] = convert_hhmmss_to_float(data_df["Time (mm:ss)"])
    data_df = merge_graftm_resources(data_df)

    # Remove DIAMOND samples
    data_df = data_df[data_df["Software"] != "DIAMOND"]

    ram_plot(data_df, figures_dir)

    time_plot(data_df, figures_dir)

    parallel_plot(data_df, figures_dir)

    return


if __name__ == '__main__':
    main(time_table="../manuscript/tables/runtime_log.csv",
         figures_dir="../manuscript/figures/")
