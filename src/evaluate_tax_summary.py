#!/usr/bin/env python3

import os
import glob

import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
from sklearn import linear_model as lm
from sklearn.metrics import mean_squared_error

from treesapp import jplace_utils
from evaluate_clusters import write_images_from_dict
from category_maps import _METHODS_PALETTE_MAP, _LABEL_MAT, _CATEGORIES

pio.templates.default = "plotly_white"
pd.set_option('max_columns', 10)
pd.set_option('max_rows', 100)
pd.set_option('display.width', 120)

_OUTPUT_EXP_MAP = {"MCC_graftm_0.13.0_prok": "GraftM",
                   "MCC_treesapp_0.11.2_aelw_prok": "TreeSAPP-aELW",
                   "MCC_treesapp_0.11.2_LCA_prok": "TreeSAPP-LCA",
                   "MCC_treesapp_0.11.2_maxLWR_prok": "TreeSAPP-LM"}


def load_table_files_to_dataframe(sample_table_paths: dict, value_index="Method") -> pd.DataFrame:
    data_tables = []
    for tbl_path, spl in sample_table_paths.items():
        tbl = pd.read_csv(tbl_path, sep="\t", header=0)
        tbl[value_index] = spl
        data_tables.append(tbl)
    return pd.concat(data_tables)


def load_mcc_tables(analysis_dir: str) -> pd.DataFrame:
    tbl_name = "prokaryotic_e5.proteomes_MCC_table.tsv"
    output_dirs = list(_OUTPUT_EXP_MAP.keys())
    table_paths = {}
    for dir_name in output_dirs:
        table_paths[os.path.join(analysis_dir, dir_name, tbl_name)] = _OUTPUT_EXP_MAP[dir_name]
    return load_table_files_to_dataframe(table_paths)


def summarise_classification_performance(df: pd.DataFrame, output_dir: str) -> None:
    sig_digs = 4
    df["FPR"] = round(df["False.Pos"] / (df["False.Pos"] + df["True.Neg"]), sig_digs)
    df["FDR"] = round(df["False.Pos"] / (df["False.Pos"] + df["True.Pos"]), sig_digs)
    df["PPV"] = round(df["True.Pos"] / (df["True.Pos"] + df["False.Pos"]), sig_digs)
    # df["NPV"] = df["True.Neg"] / (df["True.Neg"] + df["False.Neg"])
    df.to_csv(path_or_buf=os.path.join(output_dir, "MCC_data.tsv"), sep="\t", index=False, header=True)
    return


def load_classification_tables(analysis_dir: str, output_dirs: list) -> pd.DataFrame:
    tbl_name = "prokaryotic_e5.proteomes_classifications.tsv"
    table_paths = {}
    for dir_name in output_dirs:
        table_paths[os.path.join(analysis_dir, dir_name, tbl_name)] = _OUTPUT_EXP_MAP[dir_name]
    return load_table_files_to_dataframe(table_paths)


def mcc_line_plot(df: pd.DataFrame, output_dir: str) -> None:
    plot_prefix = os.path.join(output_dir, "MCC_line")
    line_plt = px.line(df,
                       x="Tax.dist", y="MCC",
                       color="Method", line_group="Method",
                       color_discrete_map=_METHODS_PALETTE_MAP,
                       category_orders=_CATEGORIES,
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


def fit_line(x_train: np.array, y_train: np.array, num=100):
    model = lm.LinearRegression()
    model.fit(X=x_train.reshape(-1, 1),
              y=y_train)
    x_range = np.linspace(start=x_train.min(), stop=x_train.max(), num=num)
    y_range = model.predict(x_range.reshape(-1, 1))
    y_pred = model.predict(x_train.reshape(-1, 1))
    residuals = y_pred - y_train
    return x_range, y_range, y_pred, residuals


def generate_rmse_dataframe(df: pd.DataFrame, group_var: str, rep_var: str, value_var: str) -> pd.DataFrame:
    rmse_dat = {"Method": [],
                "RMSE": []}
    for group in pd.unique(df[group_var]):
        group_df = df[df[group_var] == group]

        for rep_name in pd.unique(group_df[rep_var]):
            rep_df = group_df[group_df[rep_var] == rep_name]
            rmse_dat["Method"].append(group)
            preds = rep_df[value_var]
            rmse_dat["RMSE"].append(mean_squared_error(y_true=[0]*len(preds), y_pred=preds, squared=False))

    return pd.DataFrame(rmse_dat)


def rmse_bar_plot(df: pd.DataFrame, output_dir: str) -> None:
    rmse_df = generate_rmse_dataframe(df, group_var="Method", rep_var="TrueLineage", value_var="TaxDist")
    fig = px.bar(rmse_df.groupby("Method").agg({"RMSE": ['mean', 'std']})["RMSE"].reset_index(),
                 x="Method", y="mean", error_y="std",
                 color="Method",
                 labels={"mean": "RMSE"},
                 category_orders=_CATEGORIES,
                 color_discrete_map=_METHODS_PALETTE_MAP)
    bar_prefix = os.path.join(output_dir, "tax_dist_rmse")
    write_images_from_dict({bar_prefix: fig}, fig_scale=2)
    return


def evo_dist_plot(df: pd.DataFrame, output_dir: str) -> None:
    # How well do the two methods handle large evolutionary distances?
    # Is taxonomic distance driven by evolutionary distance?
    fig = go.Figure(layout_xaxis_range=[0, 8])
    residual_dict = {"Method": [],
                     "prediction": [],
                     "residual": [],
                     "x": []}
    i = 1
    for group in df["Method"].unique():
        sub_df = df[df["Method"] == group]
        x, y, preds, res = fit_line(x_train=np.array(sub_df.EvoDist), y_train=sub_df.TaxDist)
        residual_dict["prediction"] += list(preds)
        residual_dict["residual"] += list(res)
        residual_dict["x"] += list(np.array(sub_df.EvoDist))
        residual_dict["Method"] += [group] * len(sub_df)
        if i == 1:
            side = 'positive'
        else:
            side = 'negative'
        fig.add_traces([go.Scatter(x=x, y=y, name=group, showlegend=False,
                                   line=dict(color=_METHODS_PALETTE_MAP[group], width=4)),
                        go.Violin(x=sub_df['EvoDist'],
                                  y=sub_df['TaxDist'],
                                  legendgroup=group, scalegroup=group, name=group,
                                  fillcolor=_METHODS_PALETTE_MAP[group],
                                  orientation='h', side=side, points=False, meanline={"visible": True})
                        ])
        i += 1
    fig.update_traces(line_color="gray", selector=dict(type='violin'))
    fig.update_traces(marker_line_color="gray", selector=dict(type='scatter'))
    fig.update_layout(violingap=0, violinmode='overlay',
                      xaxis_title="Evolutionary distance",
                      yaxis_title="Taxonomic distance")

    res_viol = px.scatter(pd.DataFrame(residual_dict),
                          x='prediction', y='residual',
                          marginal_y='violin',
                          color_discrete_map=_METHODS_PALETTE_MAP,
                          category_orders=_CATEGORIES,
                          labels=_LABEL_MAT,
                          color='Method', trendline='ols')

    boxes = px.box(df,
                   x="RefPkg", y="TaxDist",
                   color="Method", category_orders=_CATEGORIES,
                   labels=_LABEL_MAT, color_discrete_map=_METHODS_PALETTE_MAP,
                   points='outliers', boxmode="group")
    boxes.update_xaxes(categoryorder='category ascending',
                       tickangle=45,
                       title_font=None,
                       tickfont_size=6,
                       title_standoff=10)
    boxes.update_traces(marker_line_color='rgb(105,105,105)',
                        marker_line_width=0.5,
                        marker_size=3, selector=dict(type='box'))

    reg_prefix = os.path.join(output_dir, "evo_tax_corr")
    res_prefix = os.path.join(output_dir, "residuals")
    box_prefix = os.path.join(output_dir, "tax_summary_dist_boxes")
    write_images_from_dict({reg_prefix: fig, res_prefix: res_viol, box_prefix: boxes}, fig_scale=2)
    return


def summarise_jplace_placements(analysis_dir: str, tables_dir: str):
    output_dirs = ["MCC_treesapp_0.11.2_aelw_prok"]
    pplace_dat = {"RefPkg": [],
                  "Placements": [],
                  "LWR": []}

    for dir_name in output_dirs:
        pplace_dir = os.path.join(analysis_dir, dir_name, "TreeSAPP_output", "intermediates", "place")
        for jplace in glob.glob(pplace_dir + os.sep + "*jplace"):
            place_dat = jplace_utils.jplace_parser(jplace)
            pqueries = jplace_utils.demultiplex_pqueries(place_dat)
            for pquery in pqueries:
                for pplace in pquery.placements:
                    pplace_dat["RefPkg"].append(pquery.ref_name)
                    pplace_dat["Placements"].append(len(pquery.placements))
                    pplace_dat["LWR"].append(-1*float(pplace.like_weight_ratio))
    df = pd.DataFrame(pplace_dat)
    df_perc = df.groupby("RefPkg").quantile([.25, .5, .75]).reset_index()
    df_perc["LWR"] = round(df_perc["LWR"]*-1, 2)
    df_perc["Placements"] = df_perc["Placements"].astype(int)
    df_perc.rename(columns={'level_1': "Percentile"})
    df_perc.to_csv(os.path.join(tables_dir,
                                "placements_summary.csv"),
                   index=False,
                   mode='w')
    return


def evaluate_taxonomic_summary_methods(analysis_dir: str, figures_dir: str, tables_dir: str) -> None:
    summarise_jplace_placements(analysis_dir, tables_dir)

    mcc_df = load_mcc_tables(analysis_dir)
    summarise_classification_performance(mcc_df, tables_dir)

    mcc_line_plot(mcc_df, figures_dir)

    pquery_df = load_classification_tables(analysis_dir, output_dirs=["MCC_treesapp_0.11.2_aelw_prok",
                                                                      "MCC_treesapp_0.11.2_maxLWR_prok",
                                                                      "MCC_treesapp_0.11.2_LCA_prok"]).dropna(axis=0)
    rmse_bar_plot(pquery_df, figures_dir)

    pquery_df = load_classification_tables(analysis_dir, output_dirs=["MCC_treesapp_0.11.2_aelw_prok",
                                                                      "MCC_treesapp_0.11.2_maxLWR_prok"]).dropna(axis=0)
    rmse_bar_plot(pquery_df, figures_dir)
    evo_dist_plot(pquery_df, figures_dir)

    return


if __name__ == '__main__':
    evaluate_taxonomic_summary_methods(analysis_dir="../tax_summary/",
                                       figures_dir="../manuscript/figures/",
                                       tables_dir="../manuscript/tables/")
