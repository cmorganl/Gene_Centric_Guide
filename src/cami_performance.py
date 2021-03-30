#!/usr/bin/env python

import os
import re
import glob

import plotly.express as px
import plotly.graph_objects as go
import ptitprince as pt
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import upsetplot

from treesapp import file_parsers
from treesapp import classy
from treesapp import taxonomic_hierarchy as ts_th
from treesapp import entrez_utils as ts_entrez
from treesapp import phylo_seq as ts_phyloseq
from treesapp import refpkg as ts_refpkg
from treesapp import lca_calculations as ts_lca
from evaluate_clusters import ClusterExperiment

_LABEL_MAT = {"Length": "Query sequence length",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method",
              "TaxDist": "Taxonomic distance",
              "Size": "Number of queries"}


def generate_entrez_queries(taxids) -> list:
    entrez_queries = []
    unique_taxids = set(taxids)
    for taxid in unique_taxids:
        e_record = ts_entrez.EntrezRecord('NA', 'NA')
        e_record.ncbi_tax = taxid
        e_record.tracking_stamp()
        entrez_queries.append(e_record)
    return entrez_queries


def match_taxid_to_seqname(cami_pqueries: dict, taxon_mapping_file: str) -> set:
    # Match the CAMI sequence names to NCBI taxids using the mapping file
    assigned_seqnames = set()
    taxids = set()
    for pqueries in cami_pqueries.values():
        orf_names = [pq.seq_name for pq in pqueries]
        for orf in orf_names:
            # Convert ORF and read names to contig/mate-pair names
            assigned_seqnames.add(re.sub(r'/1$', '', '_'.join(orf.split('_')[:-1])))

    seqname_taxid_map = {}
    with open(taxon_mapping_file, 'r') as taxa_map:
        for line in taxa_map:
            line = line.strip()
            if not line or line[0] == '@':
                continue
            if len(assigned_seqnames) == 0:
                break
            try:
                seq_name, _, taxid, _ = line.split("\t")
            except ValueError:
                if len(line.split("\t")) == 3:
                    seq_name, _, taxid = line.split("\t")
            if seq_name in assigned_seqnames:
                seqname_taxid_map[seq_name] = taxid
                assigned_seqnames.remove(seq_name)

    if len(assigned_seqnames) > 0:
        print("{} sequence names were not found in the mapping file '{}'."
              "".format(len(assigned_seqnames), taxon_mapping_file))

    # Stash the NCBI taxid in PQuery
    for pqueries in cami_pqueries.values():
        for pquery in pqueries:  # type: ts_phyloseq.PQuery
            try:
                pquery.ncbi_tax = seqname_taxid_map[re.sub(r'/1$', '', '_'.join(pquery.seq_name.split('_')[:-1]))]
                taxids.add(pquery.ncbi_tax)
            except KeyError:
                pquery.ncbi_tax = None

    return taxids


def retrieve_lineages(cami_pqueries: dict, sample_taxon_map: dict) -> dict:
    """
    Determines the format of the query sequence header to extract the accession and/or NCBI taxid then
    queries the Entrez database using these information to collect the taxonomic lineage for each unique NCBI taxid
    NCBI taxid and lineage information are stored in self.tax_lineage_map

    :return: A dictionary mapping NCBI taxids to taxonomic lineages
    """
    # Gather the unique taxonomy IDs and store in EntrezRecord instances
    t_hierarchy = ts_th.TaxonomicHierarchy()
    unknowns = {}
    tax_lineage_map = {}

    # Map the classified sequence names to NCBI taxids using CAMI mapping file, create EntrezRecords for querying
    taxids = set()
    for sample, taxon_map in sample_taxon_map.items():
        taxids.update(match_taxid_to_seqname(cami_pqueries[sample], taxon_map))
    entrez_records = generate_entrez_queries(taxids)

    # Query the Entrez database for these unique taxonomy IDs
    ts_entrez.get_multiple_lineages(entrez_records, t_hierarchy, "prot")

    for e_record in entrez_records:  # type: ts_entrez.EntrezRecord
        if not e_record.lineage:
            try:
                unknowns[e_record.ncbi_tax] += 1
            except KeyError:
                unknowns[e_record.ncbi_tax] = 1
            tax_lineage_map[e_record.ncbi_tax] = ''
        else:
            tax_lineage_map[e_record.ncbi_tax] = "r__Root; " + e_record.lineage

    if unknowns:
        warn_str = ""
        for taxid in sorted(unknowns, key=int):
            warn_str += "\t{} query sequences with unknown taxid '{}'\n".format(unknowns[taxid], taxid)
        print("Lineage information unavailable for taxonomy IDs:\n" + warn_str)
    return tax_lineage_map


def calculate_taxonomic_distances(sample_assigned_pqueries, tax_lineage_map, refpkg_dict) -> pd.DataFrame:
    samples_ar = []
    ref_pkgs_ar = []
    dists_ar = []
    query_size_ar = []

    for sample, refpkg_pqueries in sample_assigned_pqueries.items():
        for ref_pkg_name, pqueries in refpkg_pqueries.items():
            try:
                ref_pkg = refpkg_dict[ref_pkg_name]  # type: ts_refpkg.ReferencePackage
            except KeyError:
                continue
            ref_pkg.taxa_trie.build_multifurcating_trie()
            n_queries = len(pqueries)
            for pq in pqueries:  # type: ts_phyloseq.PQuery
                # Find the optimal taxonomic assignment
                if not pq.ncbi_tax:
                    continue
                true_lineage = ref_pkg.taxa_trie.clean_lineage_string(tax_lineage_map[pq.ncbi_tax])
                optimal_taxon = ts_lca.optimal_taxonomic_assignment(ref_pkg.taxa_trie.trie, true_lineage)
                if not optimal_taxon:
                    print("Optimal taxonomic assignment '{}' for {}"
                          " not found in reference hierarchy.\n".format(true_lineage,
                                                                        pq.place_name))
                    continue
                tax_dist, status = ts_lca.compute_taxonomic_distance(pq.recommended_lineage, optimal_taxon)
                if status > 0:
                    print("Lineages didn't converge between:\n"
                          "'{}' and '{}' (taxid: {})\n".format(pq.recommended_lineage, optimal_taxon, pq.ncbi_tax))
                samples_ar.append(sample)
                ref_pkgs_ar.append(ref_pkg_name)
                dists_ar.append(tax_dist)
                query_size_ar.append(n_queries)

    return pd.DataFrame(dict(RefPkg=ref_pkgs_ar, Sample=samples_ar, TaxDist=dists_ar, Size=query_size_ar))


def generate_plotly_bubble_size_legend(data: list) -> go.Figure:
    # TODO: Finish implementing this unless there's an easy way to create a bubble size legend in Plotly
    fig_legend = go.Figure(data=[go.Scatter(
        x=[1, 1, 1], y=[min(data), max(data) / 2, max(data)],
        mode='markers', fillcolor="DarkSlateGrey")])
    fig_legend.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
    return fig_legend


def plot_taxonomic_distance_bubbles(tax_dist_dat: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.Antique

    # Filter to remove RefPkgs with fewer than 100 observations
    plt_dat = tax_dist_dat[tax_dist_dat["Size"] >= 100].groupby(["RefPkg", "Sample"]).mean().reset_index()

    # Define the bubble plot size legend
    # bubble_legend = generate_plotly_bubble_size_legend(tax_dist_dat["Size"])
    # bubble_legend.show()

    bubble_plt = px.scatter(plt_dat,
                            x="RefPkg", y="TaxDist", color="Sample", size="Size",
                            size_max=20,
                            color_discrete_sequence=palette,
                            labels=_LABEL_MAT,
                            title="")
    bubble_plt.update_traces(marker=dict(line=dict(width=2,
                                                   color='DarkSlateGrey')),
                             selector=dict(mode='markers'))
    bubble_plt.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)'})
    bubble_plt.update_xaxes(showgrid=True, gridwidth=1, tickangle=45)
    bubble_plt.update_yaxes(showgrid=True, gridwidth=1, dtick=1)

    bubble_plt.write_image(os.path.join(output_dir, "tax_dist_bubbles.png"), engine="kaleido", scale=4.0)
    # bubble_plt.write_image(os.path.join(output_dir, "tax_dist_bubbles.svg"), engine="kaleido", scale=4.0)
    return


def summary_stats(tax_dist_dat: pd.DataFrame) -> None:
    alpha = 1E-5
    print("Mean:\n", tax_dist_dat.groupby("Sample")["TaxDist"].mean())
    print("Variance:\n", tax_dist_dat.groupby("Sample")["TaxDist"].var())
    k2, p = stats.normaltest(tax_dist_dat["TaxDist"])

    if p < alpha:  # null hypothesis: x comes from a normal distribution
        print("Normality test P-value = {:.6g}. These data come from a normal distribution".format(p))
    else:
        print("The null hypothesis cannot be rejected; these data are not normally distributed.")

    x = tax_dist_dat[tax_dist_dat["Sample"] == "RH_S001_merged"]["TaxDist"]
    y = tax_dist_dat[tax_dist_dat["Sample"] == "gold_standard_high_single"]["TaxDist"]
    s, p = stats.ttest_ind(x, y)
    if p < alpha:
        print("These distributions are significantly different (P-value = {:.6g})".format(p))
    else:
        print("Distributions are not significantly different.")

    return


def merge_otu_matrices(refpkg_otu_matrices: dict, reset_otu_count=True) -> pd.DataFrame:
    key_name = "#OTU_ID"
    merged_df = pd.DataFrame()
    for refpkg_name in refpkg_otu_matrices:
        # Combine the OTU matrices by columns to generate a single DataFrame
        result = refpkg_otu_matrices[refpkg_name].pop()  # type: pd.DataFrame
        for otu_mat in refpkg_otu_matrices[refpkg_name]:
            result = result.join(otu_mat.set_index(key_name), on=key_name)

        result["RefPkg"] = refpkg_name

        # Merge the combined data frame with pOTU counts for all samples
        if merged_df.empty:
            merged_df = result
        else:
            merged_df = merged_df.append(result)

    merged_df = merged_df.rename(columns={"#OTU_ID": "OTU"})
    if reset_otu_count:
        merged_df["RefPkg_OTU"] = [x for x in range(0, len(merged_df))]
    if not reset_otu_count:
        raise AssertionError()
        # merged_df.assign(RefPkg_OTU=lambda x: x.RefPkg + '-' + x.OTU)

    # Drop the irrelevant and/or redundant variables
    merged_df = merged_df.drop(["OTU", "RefPkg"], axis=1)
    return merged_df


def get_potu_data(phylotu_outputs: list) -> (pd.DataFrame, pd.DataFrame):
    samples_ar = []
    potu_count_ar = []
    refpkg_ar = []
    if len(phylotu_outputs) == 0:
        print("No treesapp phylotu output directories were found.")
        raise AssertionError

    refpkg_otu_matrices = {}

    for phylotu_dir in phylotu_outputs:
        phylotu_exp = ClusterExperiment(directory=phylotu_dir)
        if not phylotu_exp.test_files():
            continue
        if not phylotu_exp.load_cluster_assignments():
            continue
        phylotu_exp.merge_query_clusters()
        try:
            refpkg_otu_matrices[phylotu_exp.pkg_name].append(phylotu_exp.load_otu_matrix())
        except KeyError:
            refpkg_otu_matrices[phylotu_exp.pkg_name] = [phylotu_exp.load_otu_matrix()]
        potu_count_ar.append(phylotu_exp.get_unique_potus())
        samples_ar.append(phylotu_dir.split(os.sep)[-2])
        refpkg_ar.append(phylotu_exp.pkg_name)

    merged_otu_df = merge_otu_matrices(refpkg_otu_matrices)

    return pd.DataFrame(dict(Sample=samples_ar, OTUs=potu_count_ar, RefPkg=refpkg_ar)), merged_otu_df


def plot_phylotu_boxes(potu_df: pd.DataFrame, output_dir: str) -> None:
    violin_plt = px.violin(potu_df,
                           x="Sample", y="OTUs",
                           box=True, points="all", range_y=[0, 310],
                           labels=_LABEL_MAT,
                           title="")
    violin_plt.update_traces(marker=dict(color='DarkSlateGrey'))
    violin_plt.update_yaxes(showgrid=True, gridwidth=1, dtick=100)

    violin_plt.write_image(os.path.join(output_dir, "pOTU_count_violin.png"), engine="kaleido", scale=4.0)
    return


def generate_otu_contents(potu_df, index=None) -> dict:
    memberships = {}
    potu_df = potu_df.reset_index()
    for i in potu_df.columns:
        if i in [index, "index"]:
            continue
        memberships[i] = [potu_df[index][x] for x, y in potu_df[i].items() if y > 0]
    return memberships


def plot_phylotu_upset(potu_df: pd.DataFrame, output_dir: str) -> None:
    index_name = "RefPkg_OTU"
    tmp_dat = pd.DataFrame()
    tmp_dat[index_name] = potu_df[index_name]
    otu_memberships = generate_otu_contents(potu_df, index=index_name)
    # Subset tmp_dat

    upset_dat = upsetplot.from_contents(otu_memberships)
    upsetplot.plot(upset_dat, sort_by="cardinality", show_percentages=True)
    # plt.savefig(fname=os.path.join(output_dir, "UpSet_pOTUs.svg"))
    plt.savefig(fname=os.path.join(output_dir, "UpSet_pOTUs.png"))

    return


def plot_rainclouds(potu_df: pd.DataFrame, output_dir: str) -> None:
    dx = "Sample"
    dy = "OTUs"
    pal = "Set2"
    sigma = .2
    f, ax = plt.subplots(figsize=(7, 5))
    ax = pt.RainCloud(x=dx, y=dy, data=potu_df, palette=pal,
                      bw=sigma, width_viol=.6, move=.2, ax=ax, orient="h")
    plt.savefig(fname=os.path.join(output_dir, "raining_pOTUs.png"))
    # TODO: add lines between reference package points
    return


def main(root_dir):
    data_dir = os.path.join(root_dir, "CAMI_experiments") + os.sep
    fig_dir = os.path.join(root_dir, "manuscript", "figures") + os.sep
    assign_outputs = {"gold_standard_high_single": data_dir + "gsa_mapping_pool.binning",
                      "RH_S001_merged": data_dir + "gs_read_mapping_1.binning",
                      "RH_S001_forward": data_dir + "gs_read_mapping_1.binning"}
    refpkg_dir = "refpkgs"
    refpkg_dict = ts_refpkg.gather_ref_packages(refpkg_data_dir=data_dir + refpkg_dir)

    sample_assigned_pqueries = {}

    # Count the number of pOTUs for each RefPkg in each output
    potu_df, potu_mat = get_potu_data(glob.glob(data_dir + "*/phylotu_out_*"))
    plot_phylotu_upset(potu_mat, fig_dir)

    # Read the assignments from each of the assign_outputs directories
    for assign_out in assign_outputs:
        sample_assigned_pqueries.update({assign_out:
                                             file_parsers.load_classified_sequences_from_assign_output(assign_output_dir=data_dir + assign_out,
                                                                                                       assigner_cls=classy.TreeSAPP("cami"))})

    # Fetch the taxonomic lineages for the unique NCBI taxonomy IDs from the classified query sequences
    tax_lineage_map = retrieve_lineages(sample_assigned_pqueries, assign_outputs)

    # Calculate the taxonomic distances by RefPkg and output, return a pandas dataframe
    tax_dist_dat = calculate_taxonomic_distances(sample_assigned_pqueries, tax_lineage_map, refpkg_dict)

    # Plot the taxonomic distances
    plot_taxonomic_distance_bubbles(tax_dist_dat, fig_dir)

    # Report the mean and variance of each sample and test for significance
    summary_stats(tax_dist_dat)

    # Plot the pOTUs
    plot_phylotu_boxes(potu_df, fig_dir)

    plot_rainclouds(potu_df, fig_dir)

    return


if __name__ == "__main__":
    main("/media/connor/Rufus/Gene_Centric_Guide/")
