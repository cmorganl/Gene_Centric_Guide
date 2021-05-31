#!/usr/bin/env python3

import glob
import os
import sys
import re
import time
from collections import Counter
import itertools

import pandas as pd
import treesapp as ts
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from bokeh import palettes
from tqdm import tqdm
from sklearn import metrics
from scipy.signal import savgol_filter
from scipy.stats import f_oneway, normaltest, ttest_ind
import statsmodels.api as sm

pio.templates.default = "plotly_white"
_REFPKG_DIR = "clustering_refpkgs"
_LABEL_MAT = {"Length": "Protein length percentage",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method",
              "Spline": "Cluster accuracy",
              "WTD": "Taxonomic Distinctness"}
_RANKS = ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
_REFPKGS = ["RecA", "RpoB", "PF01655", "NifH", "SoxY", "McrA"]
_CATEGORIES = {"Clustering": ["de_novo-aln", "de_novo-psc", "ref_guided"],
               "RefPkg": _REFPKGS,
               "Rank": _RANKS,
               "Resolution": _RANKS}
_RANK_PAL = palettes.linear_palette(px.colors.diverging.PuOr, len(_RANKS))
_RANK_PALETTE_MAP = {_RANKS[i]: _RANK_PAL[i] for i in range(0, len(_RANKS))}
_REFPKG_PAL = px.colors.qualitative.Safe
_REFPKG_PALETTE_MAP = {_REFPKGS[i]: _REFPKG_PAL[i] for i in range(0, len(_REFPKGS))}


class ClusterExperiment:
    length_parse_re = re.compile(r"length_(\d+)$")
    qseq_slider_re = re.compile(r"_sliding:\d+-\d+$")

    def __init__(self, directory: str):
        """Create an instance representing a single treesapp phylotu output."""
        self.dir_path = directory
        self.pquery_assignments_file = os.path.join(self.dir_path, "final_outputs", "phylotu_pquery_assignments.tsv")
        self.matrix_file = os.path.join(self.dir_path, "final_outputs", "phylotu_matrix.tsv")
        self.seq_length = None
        self.cluster_resolution = None
        self.cluster_mode = None
        self.precluster_mode = None
        self.pkg_name = None
        self.ref_pkg = ts.refpkg.ReferencePackage()
        self.cluster_assignments = {}
        self.cluster_members = {}
        self.phylo_place_map = {}
        self.entrez_query_dict = {}
        self.query_taxon_map = {}
        self.renamed_taxa = {}

        return

    def __str__(self):
        return "'{}' ClusterExperiment at resolution '{}'".format(self.pkg_name, self.cluster_resolution)

    def test_files(self) -> bool:
        if not os.path.exists(self.pquery_assignments_file):
            # raise AssertionError("No phylotu_pquery_assignments.tsv file found for '{}'.".format(self.dir_path))
            return False
        return True

    def parse_seq_length(self):
        for dir_name in self.dir_path.split(os.sep):
            if self.length_parse_re.match(dir_name):
                self.seq_length = int(self.length_parse_re.match(dir_name).group(1))
                break
        if not self.seq_length:
            raise AssertionError("Unable to determine the sequence length used in '{}'.".format(self.dir_path))
        try:
            hmm_len = self.ref_pkg.hmm_length()
        except ValueError:
            raise RuntimeError("Unable to load reference package's HMM length.")
        if self.seq_length > hmm_len:
            self.seq_length = hmm_len
        return

    def set_precluster_mode(self, append=False):
        for word in os.path.basename(self.dir_path).split('_'):
            if word in ["aln", "psc"]:
                self.precluster_mode = word
                if append:
                    self.cluster_mode += "-" + self.precluster_mode
        if self.precluster_mode is None:
            self.precluster_mode = "NA"
        return

    @staticmethod
    def validate_unity(collection_of_things: set):
        if len(collection_of_things) > 1:
            raise AssertionError("Expected a single element in {}.".format(collection_of_things))
        elif len(collection_of_things) == 0:
            raise AssertionError("Empty set.")
        return collection_of_things.pop()

    def load_cluster_assignments(self, delim="\t") -> int:
        assignments_handler = open(self.pquery_assignments_file)
        cluster_res = set()
        cluster_modes = set()
        pkg_names = set()
        header = assignments_handler.readline().strip().split(delim)
        if header != ['PQuery', 'RefPkg', 'Resolution', 'Mode', 'pOTU', "d_Distal", "d_Pendant", "d_MeanTip"]:
            raise AssertionError("Unexpected table format detected by header differences.")

        line = assignments_handler.readline()
        if not line:
            return 0
        while line:
            qseq_name, ref_pkg, resolution, mode, cluster_id, distal, pendant, mean_tip = line.strip().split(delim)
            self.cluster_assignments[qseq_name] = cluster_id
            pplace = ts.phylo_seq.PhyloPlace()
            pplace.name = qseq_name
            pplace.set_attribute_types(pquery_dists=','.join([distal, pendant, mean_tip]))
            self.phylo_place_map[qseq_name] = pplace
            cluster_res.add(resolution)
            cluster_modes.add(mode)
            pkg_names.add(ref_pkg)
            line = assignments_handler.readline()
        assignments_handler.close()
        self.cluster_resolution = self.validate_unity(cluster_res)
        self.cluster_mode = self.validate_unity(cluster_modes)
        self.pkg_name = self.validate_unity(pkg_names)

        return 1

    def convert_length_to_percentage(self) -> None:
        hmm_perc = float(100 * int(self.seq_length) / self.ref_pkg.hmm_length())
        if hmm_perc > 100.0:
            self.seq_length = 100.0
        else:
            self.seq_length = hmm_perc
        return

    def load_otu_matrix(self, delim="\t") -> pd.DataFrame:
        return pd.read_csv(self.matrix_file, sep=delim)

    def load_ref_pkg(self, refpkgs_root_dir: str):
        self.ref_pkg.f__pkl = os.path.join(refpkgs_root_dir, self.pkg_name + "_build.pkl")
        self.ref_pkg.slurp()
        return

    def merge_query_clusters(self):
        merged_cluster_assignments = {}
        merged_pplaces = {}
        for qseq_name, cluster in self.cluster_assignments.items():
            pplace = self.phylo_place_map[qseq_name]
            if self.qseq_slider_re.search(qseq_name):
                qseq_name = self.qseq_slider_re.sub(repl='', string=qseq_name)
            try:
                merged_cluster_assignments[qseq_name].append(cluster)
                merged_pplaces[qseq_name].append(pplace)
            except KeyError:
                merged_cluster_assignments[qseq_name] = [cluster]
                merged_pplaces[qseq_name] = [pplace]
        self.cluster_assignments = merged_cluster_assignments
        self.phylo_place_map = merged_pplaces
        return

    def taxonomically_resolve_clusters(self) -> dict:
        taxon_cluster_ids = {}
        rank_depth = self.ref_pkg.taxa_trie.accepted_ranks_depths[self.cluster_resolution]
        parent_rank = get_key(self.ref_pkg.taxa_trie.accepted_ranks_depths, rank_depth - 1)
        for query in self.cluster_assignments:  # type: str
            taxon = self.query_taxon_map[query]  # type: ts.taxonomic_hierarchy.Taxon
            taxon_resolved = taxon.get_rank_in_lineage(
                self.cluster_resolution)  # type: ts.taxonomic_hierarchy.Taxon
            if not taxon_resolved:
                try:
                    taxon_resolved = taxon.get_rank_in_lineage(parent_rank).name + " " + self.cluster_resolution
                except KeyError:
                    print("Unable to find the {} rank of taxon '{}'.\n".format(self.cluster_resolution,
                                                                               taxon.name))
                    continue
                self.renamed_taxa[taxon.name] = taxon_resolved
            else:
                taxon_resolved = taxon_resolved.name

            if taxon_resolved not in taxon_cluster_ids:
                taxon_cluster_ids[taxon_resolved] = []
            taxon_cluster_ids[taxon_resolved] += self.cluster_assignments[query]

        return taxon_cluster_ids

    def match_query_to_taxon_cluster(self, seq_name: str) -> str:
        taxon = self.query_taxon_map[seq_name]
        taxon_resolved = taxon.get_rank_in_lineage(self.cluster_resolution)
        if not taxon_resolved or taxon.name in self.renamed_taxa:
            return self.renamed_taxa[taxon.name]
        else:
            return taxon_resolved.name

    @staticmethod
    def generate_true_cluster_labels(cluster_ids: list) -> list:
        cluster_counts = Counter(cluster_ids).most_common()
        consensus_cluster = cluster_counts[0][0]
        return [consensus_cluster] * len(cluster_ids)

    @staticmethod
    def report_query_completeness(cluster_assignments: list) -> float:
        true_labels = ClusterExperiment.generate_true_cluster_labels(cluster_assignments)
        return round(metrics.completeness_score(labels_true=true_labels, labels_pred=cluster_assignments), 3)

    @staticmethod
    def find_clustering_accuracy(taxon_cluster_ids: list) -> float:
        true_labels = ClusterExperiment.generate_true_cluster_labels(taxon_cluster_ids)
        h, c, v = metrics.homogeneity_completeness_v_measure(labels_true=true_labels, labels_pred=taxon_cluster_ids)
        return round(v, 3)

    def generate_entrez_queries(self) -> None:
        header_registry = ts.fasta.register_headers(list(self.cluster_assignments.keys()))
        entrez_record_dict = ts.classy.get_header_info(header_registry)
        for index, e_record in entrez_record_dict.items():  # type: ts.entrez_utils.EntrezRecord
            self.entrez_query_dict[e_record.description] = e_record
        return

    def get_unique_potus(self) -> int:
        potus = set()
        for query_potus in self.cluster_assignments.values():
            potus.update(set(query_potus))
        return len(potus)

    def load_cluster_members(self):
        for query, potus in self.cluster_assignments.items():  # type: (str, list)
            for cluster_num in set(potus):
                try:
                    self.cluster_members[cluster_num].add(query)
                except KeyError:
                    self.cluster_members[cluster_num] = set([query])
        return


def retrieve_lineages(cluster_experiments) -> dict:
    """
    Determines the format of the query sequence header to extract the accession and/or NCBI taxid then
    queries the Entrez database using these information to collect the taxonomic lineage for each unique NCBI taxid
    NCBI taxid and lineage information are stored in self.tax_lineage_map

    :return: None
    """
    # Gather the unique taxonomy IDs and store in EntrezRecord instances
    t_hierarchy = ts.taxonomic_hierarchy.TaxonomicHierarchy()
    entrez_records = []
    acc_taxon_map = dict()
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        phylotu_exp.generate_entrez_queries()
        entrez_records += [phylotu_exp.entrez_query_dict[index] for index in phylotu_exp.entrez_query_dict]
    # Query the Entrez database for these unique taxonomy IDs
    ts.entrez_utils.get_multiple_lineages(entrez_records, t_hierarchy, "prot")
    t_hierarchy.root_domains(t_hierarchy.find_root_taxon())

    for e_record in entrez_records:  # type: ts.entrez_utils.EntrezRecord
        acc_taxon_map[e_record.ncbi_tax] = t_hierarchy.get_taxon(e_record.lineage.split(t_hierarchy.lin_sep)[-1])

    return acc_taxon_map


def map_queries_to_taxa(cluster_experiments: list, accession_taxon_map: dict) -> None:
    unmapped = {}
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        for query_name in sorted(phylotu_exp.cluster_assignments, key=lambda x: int(x.split('.')[0])):
            query_acc = query_name.split('.')[0]
            if accession_taxon_map[query_acc] is None:
                try:
                    unmapped[query_acc] += 1
                except KeyError:
                    unmapped[query_acc] = 1
            else:
                phylotu_exp.query_taxon_map[query_name] = accession_taxon_map[query_acc]
    if len(unmapped) > 0:
        print("{} sequences from {} unique taxa were not downloaded".format(sum(unmapped.values()), len(unmapped)))
    return


def filter_incomplete_lineages(cluster_experiments: list) -> None:
    removed = set()
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        for query_name in sorted(phylotu_exp.cluster_assignments, key=lambda x: int(x.split('.')[0])):
            try:
                taxon = phylotu_exp.query_taxon_map[query_name]
            except KeyError:
                removed.add(query_name)
                continue
            if not taxon.get_rank_in_lineage(phylotu_exp.cluster_resolution):
                removed.add(query_name)

        for query_name in removed:
            try:
                phylotu_exp.cluster_assignments.pop(query_name)
            except KeyError:
                pass
    print("Removed {} query sequences with incomplete or missing lineages".format(len(removed)))
    return


def filter_by_seq_lengths(cluster_experiments: list, min_perc=20) -> None:
    i = 0
    rm = 0
    while i < len(cluster_experiments):
        phylotu_exp = cluster_experiments[i]  # type: ClusterExperiment
        phylotu_exp.convert_length_to_percentage()
        if phylotu_exp.seq_length < min_perc:
            cluster_experiments.pop(i)
            rm += 1
        else:
            i += 1
    if rm:
        print("Removed {} datasets where the percentage of full-length sequence was less than {}%".format(rm, min_perc))
    return


def get_key(a_dict: dict, val):
    for key, value in a_dict.items():
        if val == value:
            return key


def prepare_clustering_cohesion_dataframe(cluster_experiments: list) -> pd.DataFrame:
    ref_pkgs = []
    resolutions = []
    lengths = []
    cluster_modes = []
    cohesion = []
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing cluster completeness dataframe")
    for phylotu_exp in sorted(cluster_experiments, key=lambda x: int(x.seq_length)):  # type: ClusterExperiment
        for seq_name in phylotu_exp.cluster_assignments:
            cohesion.append(phylotu_exp.report_query_completeness(phylotu_exp.cluster_assignments[seq_name]))
            cluster_modes.append(phylotu_exp.cluster_mode)
            ref_pkgs.append(phylotu_exp.pkg_name)
            resolutions.append(phylotu_exp.cluster_resolution)
            lengths.append(phylotu_exp.seq_length)
        p_bar.update()
    p_bar.close()

    return pd.DataFrame(dict(RefPkg=ref_pkgs, Completeness=cohesion,
                             Clustering=cluster_modes, Resolution=resolutions, Length=lengths))


def prepare_clustering_accuracy_dataframe(cluster_experiments: list) -> pd.DataFrame:
    re_map = {}

    # data arrays
    refpkgs = []
    lengths = []
    clusters = []
    resos = []
    taxa = []
    accurs = []
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing cluster accuracy dataframe")
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        if phylotu_exp.cluster_resolution not in re_map:
            re_map[phylotu_exp.cluster_resolution] = {}
        taxon_cluster_ids = phylotu_exp.taxonomically_resolve_clusters()
        re_map[phylotu_exp.cluster_resolution].update(phylotu_exp.renamed_taxa)
        for taxon in taxon_cluster_ids:  # type: str
            refpkgs.append(phylotu_exp.pkg_name)
            lengths.append(phylotu_exp.seq_length)
            clusters.append(phylotu_exp.cluster_mode)
            resos.append(phylotu_exp.cluster_resolution)
            taxa.append(taxon)
            accurs.append(phylotu_exp.find_clustering_accuracy(taxon_cluster_ids[taxon]))
        p_bar.update()
    p_bar.close()

    for cluster_res in re_map:  # type: str
        if re_map[cluster_res]:  # type: dict
            for taxon_name, parent_name in re_map[cluster_res].items():
                print("Set missing taxonomic {} of {} to '{}'.".format(cluster_res, taxon_name, parent_name))

    return pd.DataFrame(dict(RefPkg=refpkgs, Length=lengths, Clustering=clusters,
                             Resolution=resos, Taxon=taxa, Accuracy=accurs))


def prepare_taxonomic_summary_dataframe(cluster_experiments: list) -> pd.DataFrame:
    pqueries_accounted = set()
    data_dict = {"RefPkg": [],
                 "Taxon": [],
                 "Related": [],
                 "Rank": []}
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing taxonomic summary dataframe")
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        if set(phylotu_exp.query_taxon_map).issubset(pqueries_accounted):
            p_bar.update()
            continue
        phylotu_exp.ref_pkg.taxa_trie.build_multifurcating_trie()
        taxonomic_tree = phylotu_exp.ref_pkg.all_possible_assignments()
        for query_name in phylotu_exp.cluster_assignments:
            if query_name in pqueries_accounted:
                continue
            else:
                pqueries_accounted.add(query_name)
            try:
                query_taxon = phylotu_exp.query_taxon_map[query_name]  # type: ts.taxonomic_hierarchy.Taxon
            except KeyError:
                continue
            data_dict["RefPkg"].append(phylotu_exp.pkg_name)
            data_dict["Taxon"].append(query_taxon.name)
            query_lineage = "; ".join([t.prefix_taxon() for t in query_taxon.lineage() if
                                       t.rank in phylotu_exp.ref_pkg.taxa_trie.accepted_ranks_depths])
            relative = ts.lca_calculations.optimal_taxonomic_assignment(query_taxon=query_lineage, trie=taxonomic_tree)
            closest_taxon = relative.split(phylotu_exp.ref_pkg.taxa_trie.lin_sep)[-1]
            data_dict["Related"].append(closest_taxon)
            rank = phylotu_exp.ref_pkg.taxa_trie.get_taxon(closest_taxon).rank
            data_dict["Rank"].append(rank)
        p_bar.update()
    p_bar.close()

    return pd.DataFrame(data_dict)


def prepare_relationships_dataframe(cluster_experiments: list) -> pd.DataFrame:
    refpkg_pqueries = {}
    pqueries_accounted = set()
    rank_depths = None
    data_dict = {"RefPkg": [],
                 "Rank": [],
                 "Count": []}

    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        if not rank_depths:
            rank_depths = phylotu_exp.ref_pkg.taxa_trie.accepted_ranks_depths
        if set(phylotu_exp.query_taxon_map).issubset(pqueries_accounted):
            continue
        if phylotu_exp.pkg_name not in refpkg_pqueries:
            refpkg_pqueries[phylotu_exp.pkg_name] = []
        for pquery in phylotu_exp.query_taxon_map:
            if pquery not in pqueries_accounted:
                refpkg_pqueries[phylotu_exp.pkg_name].append(phylotu_exp.query_taxon_map[pquery])
                pqueries_accounted.add(pquery)

    p_bar = tqdm(ncols=100)
    for ref_pkg, taxa in refpkg_pqueries.items():
        p_bar.reset(total=len(taxa))
        p_bar.set_description(desc="Counting taxonomic relationships for {}".format(ref_pkg))
        counts = {}
        r_taxon = taxa.pop()  # type: ts.taxonomic_hierarchy.Taxon
        while taxa:
            lcr = 'root'  # lowest common rank
            for q_taxon in taxa:
                lca = ts.taxonomic_hierarchy.Taxon.lca(r_taxon, q_taxon)
                while lca.rank not in rank_depths:
                    lca = lca.parent
                if rank_depths[lca.rank] > rank_depths[lcr]:
                    lcr = lca.rank
                if lcr == "species":
                    break
            try:
                counts[lcr] += 1
            except KeyError:
                counts[lcr] = 1
            r_taxon = taxa.pop()
            p_bar.update()

        for rank in counts:
            data_dict["RefPkg"].append(ref_pkg)
            data_dict["Rank"].append(rank)
            data_dict["Count"].append(counts[rank])
    p_bar.close()
    return pd.DataFrame(data_dict)


def write_images_from_dict(fig_path_map: dict, fig_scale=4.0) -> None:
    for prefix, fig in fig_path_map.items():
        fig.write_image(prefix + ".png", engine="kaleido", scale=fig_scale)
        fig.update_layout(title="").write_image(prefix + ".pdf", engine="kaleido", scale=fig_scale)
    return


def taxonomic_relationships_plot(hierarchy_df: pd.DataFrame, output_dir: str) -> None:
    fig = px.bar(hierarchy_df,
                 x="RefPkg", y="Count",
                 color="Rank",
                 barmode="group",
                 color_discrete_map=_RANK_PALETTE_MAP,
                 category_orders=_CATEGORIES,
                 labels=_LABEL_MAT,
                 title="Taxonomic relationships between EggNOG query sequences")
    fig.update_traces(marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)

    write_images_from_dict({os.path.join(output_dir, "relationship_bars"): fig})
    return


def prepare_taxonomic_distinctness_dataframe(cluster_experiments: list) -> pd.DataFrame:
    data_dict = {"RefPkg": [],
                 "Resolution": [],
                 "Clustering": [],
                 "Length": [],
                 "Cluster": [],
                 "WTD": []}
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing taxonomic distinctness dataframe")
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        phylotu_exp.load_cluster_members()
        phylotu_exp.ref_pkg.taxa_trie.validate_rank_prefixes()
        for potu, members in phylotu_exp.cluster_members.items():  # type: (str, set)
            taxa = {m: phylotu_exp.query_taxon_map[m] for m in members}
            data_dict["RefPkg"].append(phylotu_exp.pkg_name)
            data_dict["Resolution"].append(phylotu_exp.cluster_resolution)
            data_dict["Clustering"].append(phylotu_exp.cluster_mode)
            data_dict["Length"].append(phylotu_exp.seq_length)
            data_dict["Cluster"].append(potu)
            data_dict["WTD"].append(ts.lca_calculations.taxonomic_distinctness(taxa,
                                                                               rank=phylotu_exp.cluster_resolution,
                                                                               rank_depths=phylotu_exp.ref_pkg.taxa_trie.accepted_ranks_depths))
        p_bar.update()
    p_bar.close()
    return pd.DataFrame(data_dict)


def prepare_evodist_accuracy_dataframe(cluster_experiments: list) -> pd.DataFrame:
    data_dict = {"RefPkg": [],
                 "Clustering": [],
                 "Length": [],
                 "Query": [],
                 "Proper": [],
                 "Pendant": [],
                 "Distal": [],
                 "Total": []}
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing evolutionary distance dataframe")
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        taxon_cluster_ids = phylotu_exp.taxonomically_resolve_clusters()
        true_clusters = {taxon: phylotu_exp.generate_true_cluster_labels(taxon_cluster_ids[taxon]) for
                         taxon in taxon_cluster_ids}
        for query_name, fragments in phylotu_exp.cluster_assignments.items():  # type: (str, list)
            repr_taxon = phylotu_exp.match_query_to_taxon_cluster(query_name)
            if not repr_taxon:
                continue
            data_dict["RefPkg"] += [phylotu_exp.pkg_name] * len(fragments)
            data_dict["Clustering"] += [phylotu_exp.cluster_mode] * len(fragments)
            data_dict["Length"] += [phylotu_exp.seq_length] * len(fragments)
            data_dict["Query"] += [query_name] * len(fragments)
            for potu in fragments:
                if potu in true_clusters[repr_taxon]:
                    data_dict["Proper"].append(True)
                else:
                    data_dict["Proper"].append(False)
            for pplace in phylotu_exp.phylo_place_map[query_name]:  # type: ts.phylo_seq.PhyloPlace
                data_dict["Pendant"].append(pplace.pendant_length)
                data_dict["Distal"].append(pplace.distal_length)
                data_dict["Total"].append(pplace.total_distance())
        p_bar.update()
    p_bar.close()
    evo_df = pd.DataFrame(data_dict)
    return evo_df.groupby(["Clustering", "RefPkg", "Query", "Proper"]).mean(numeric_only=True).reset_index()


def sequence_cohesion_plots(frag_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.T10
    line_path = os.path.join(output_dir, "completeness_line")
    bar_path = os.path.join(output_dir, "completeness_bar")
    violin_path = os.path.join(output_dir, "completeness_violin")

    line_plt = px.line(frag_df.groupby(["Clustering", "Length"]).mean().reset_index(),
                       x="Length", y="Completeness",
                       color="Clustering", line_group="Clustering",
                       color_discrete_sequence=palette,
                       labels=_LABEL_MAT,
                       title="Split-sequence cluster completeness as a function of query sequence length")
    line_plt.update_traces(line=dict(width=4))
    # line_plt.show()

    bar_plt = px.bar(frag_df.groupby(["Clustering", "RefPkg"]).mean().reset_index(),
                     x="RefPkg", y="Completeness",
                     barmode="group", color="Clustering",
                     color_discrete_sequence=palette,
                     category_orders=_CATEGORIES,
                     labels=_LABEL_MAT,
                     title="Mean cluster completeness across the different reference packages")
    # bar_plt.show()

    violin_plt = px.violin(frag_df.groupby(["RefPkg", "Clustering", "Resolution", "Length"]).mean().reset_index(),
                           x="Clustering", y="Completeness", color="Clustering",
                           color_discrete_sequence=palette,
                           box=True, range_y=[0, 1.01],
                           labels=_LABEL_MAT,
                           title="Comparing cluster completeness between reference-guided and de novo methods")
    violin_plt.update_layout(showlegend=False)
    # violin_plt.show()

    write_images_from_dict({line_path: line_plt, bar_path: bar_plt, violin_path: violin_plt})

    return


def smooth_sav_golay(df: pd.DataFrame, group_vars: list, sort_var: str, num_var: str,
                     ret_var="Spline", w=5, p=3) -> pd.DataFrame:
    ret_df = pd.DataFrame({})
    var_groups = {}
    group_vals_map = {k: [] for k in group_vars}
    for var in group_vars:
        var_groups.update({v: var for v in set(df[var])})
        group_vals_map[var] = set(df[var])
    combos = itertools.product(list(group_vals_map.values()).pop(0),
                               list(group_vals_map.values()).pop(1))
    for c in combos:
        group_df = df.copy()
        for val in c:
            group_df = group_df[group_df[var_groups[val]] == val].sort_values(by=sort_var)
        if len(group_df) == 0:
            continue
        group_df[ret_var] = savgol_filter(group_df[num_var], window_length=w, polyorder=p)
        ret_df = pd.concat([ret_df, group_df])

    return ret_df


def acc_summary_stats(accuracy_df: pd.DataFrame, kwargs: dict) -> pd.DataFrame:
    summary_dat = {"Mode": [],
                   "Resolution": [],
                   "AUC": [],
                   "Accuracy": []}
    rows = 0
    n_sig = 3
    # Summarise the accuracies across the different clustering resolutions and modes
    for cluster_mode in set(accuracy_df["Clustering"]):  # type: str
        mode_df = accuracy_df[accuracy_df["Clustering"] == cluster_mode]
        for res in set(mode_df["Resolution"]):
            res_df = mode_df[mode_df["Resolution"] == res]
            res_df = res_df.groupby(["Length", "Clustering"]).mean().reset_index().sort_values(by="Length")
            if len(res_df) <= 1:
                continue
            summary_dat["Mode"].append(cluster_mode)
            summary_dat["Resolution"].append(res)
            summary_dat["AUC"].append(round(metrics.auc(x=res_df["Length"], y=res_df["Accuracy"]), n_sig))
            summary_dat["Accuracy"].append(round(sum(res_df["Accuracy"]) / len(res_df["Accuracy"]), n_sig))
            rows += 1
    for keyword, val in kwargs.items():
        summary_dat[keyword] = [val] * rows
    summary_df = pd.DataFrame(summary_dat).sort_values(by=["Resolution", "Mode"])

    # Test whether the distributions of accuracy are significantly different
    anova_long_df = accuracy_df.groupby(["Clustering", "Length", "Resolution"]).mean().reset_index()
    print("Normality:", normaltest(anova_long_df["Accuracy"]))

    anova_df = anova_long_df.pivot(columns="Clustering", values="Accuracy")
    mode_acc_vals = []
    for mode in anova_df.columns:
        mode_acc_vals.append(anova_df[[mode]].dropna())
    print("One-way ANOVA:", f_oneway(*mode_acc_vals))
    res = sm.stats.multicomp.pairwise_tukeyhsd(endog=anova_df.stack().reset_index()[0],
                                               groups=anova_df.stack().reset_index()["Clustering"])
    print(res.summary())
    return summary_df


def evo_summary_stats(evo_dist_df: pd.DataFrame, tables_dir) -> None:
    """Write the mean evolutionary distances between the correct and incorrect clusters across different methods."""
    mean_dist_df = evo_dist_df.groupby(["Clustering", "Proper", "Query", "RefPkg"]).mean(numeric_only=True)[
        "Total"].reset_index()
    mean_dist_df = mean_dist_df.astype({"Proper": "str"})
    summary_df = mean_dist_df.groupby(["Clustering", "Proper"]).describe()
    summary_df.to_csv(os.path.join(tables_dir,
                                   "dist_summary_" + time.strftime("%d-%m-%y_%H%M%S", time.localtime()) + ".csv"),
                      index=False,
                      mode='w')
    print(summary_df)

    mode_evo_vals = []
    for mode in set(mean_dist_df["Clustering"]):
        mode_df = mean_dist_df[mean_dist_df["Clustering"] == mode]
        anova_df = mode_df.pivot(values="Total", columns="Proper")
        for outcome in anova_df.columns:
            mode_evo_vals.append(anova_df[[outcome]].dropna())
        print("T-test between cluster outcomes for '{}':".format(mode),
              ttest_ind(*mode_evo_vals))
    return


def completeness_summary_stats(comp_df: pd.DataFrame) -> None:
    comp_df.index = pd.MultiIndex.from_frame(comp_df,
                                             names=["refpkg", "completeness", "clustering", "resolution", "length"])
    print(comp_df.groupby(by=["Resolution", "Clustering"],
                          observed=True)["Completeness"].mean().reset_index())

    return


def acc_line(clustering_df: pd.DataFrame, palette, x_lims=None) -> go.Figure:
    mean_clust_df = clustering_df.groupby(["Resolution", "Length", "Clustering"]).mean(numeric_only=True).reset_index()
    smooth_clust_df = smooth_sav_golay(mean_clust_df,
                                       group_vars=["Resolution", "Clustering"],
                                       num_var="Accuracy",
                                       sort_var="Length", w=25, p=3)
    acc_line_plt = px.line(smooth_clust_df,
                           x="Length", y="Spline", color="Clustering",
                           color_discrete_sequence=palette, line_group="Clustering",
                           facet_col_spacing=0.05, range_x=x_lims,
                           line_shape="spline",
                           category_orders=_CATEGORIES,
                           facet_col="Resolution",
                           labels=_LABEL_MAT,
                           title="Effect of sequence length on clustering accuracy at different resolutions")
    acc_line_plt.update_traces(line=dict(width=4))
    acc_line_plt.update_yaxes(autorange=True)
    acc_line_plt.update_xaxes(tickangle=45,
                              title_font={"size": 10},
                              title_standoff=10)
    acc_line_plt.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # acc_line_plt.show()
    return acc_line_plt


def len_bars(clustering_df: pd.DataFrame, x_lims=None) -> go.Figure:
    fig = px.histogram(clustering_df.groupby(["RefPkg", "Length", "Taxon"]).mean().reset_index(),
                       x="Length", color="RefPkg", range_x=x_lims,
                       category_orders=_CATEGORIES,
                       labels=_LABEL_MAT,
                       color_discrete_map=_REFPKG_PALETTE_MAP,
                       height=400,
                       title="Number of taxa evaluated at each percentage of a reference package's sequence length")
    fig.update_layout(yaxis_title="No. Taxa")
    fig.update_traces(marker_line_color='rgb(105,105,105)',
                      marker_line_width=1)
    return fig


def refpkg_traces_for_plot(df: pd.DataFrame) -> list:
    traces = []
    for ref_pkg in _REFPKGS:
        rp_df = df[df["RefPkg"] == ref_pkg]
        x = []
        y = []
        for grp in set(rp_df["Clustering"]):
            x += list(rp_df[rp_df["Clustering"] == grp]["Clustering"])
            y += list(rp_df[rp_df["Clustering"] == grp]["Accuracy"])

        trace = go.Box({'x': x,
                        'y': y,
                        'name': ref_pkg,
                        'marker': {'color': _REFPKG_PALETTE_MAP[ref_pkg],
                                   'size': 6,
                                   'line': dict(width=1, color='DarkSlateGrey')
                                   }})
        trace.update({'type': 'box',
                      'boxpoints': 'all',
                      'fillcolor': 'rgba(255,255,255,0)',
                      'hoveron': 'points',
                      'hovertemplate': 'Clustering=%{x}<br>Accuracy=%{y}<extra></extra>',
                      'line': {'color': 'rgba(255,255,255,0)'},
                      'pointpos': -2,
                      'showlegend': True})
        traces.append(trace)
    return traces


def taxonomic_distinctness_plots(clustering_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.T10
    rp_box_path = os.path.join(output_dir, "wtd_refpkg_box")
    cm_box_path = os.path.join(output_dir, "wtd_method_box")
    line_path = os.path.join(output_dir, "wtd_line")
    rp_box_plt = px.box(
        clustering_df.groupby(["RefPkg", "Clustering", "Resolution", "Length"]).mean(numeric_only=True).reset_index(),
        x="Clustering", y="WTD", color="RefPkg",
        color_discrete_map=_REFPKG_PALETTE_MAP,
        category_orders=_CATEGORIES,
        facet_col="Resolution",
        labels=_LABEL_MAT,
        title="Taxonomic distinctness is greater in reference-guided clusters than de novo")
    rp_box_plt.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    cm_box_plt = px.box(
        clustering_df.groupby(["RefPkg", "Clustering", "Resolution", "Length"]).mean(numeric_only=True).reset_index(),
        x="Resolution", y="WTD", color="Clustering",
        color_discrete_sequence=palette,
        category_orders=_CATEGORIES,
        labels=_LABEL_MAT,
        title="Taxonomic distinctness is greater in reference-guided clusters than de novo")
    cm_box_plt.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # cm_box_plt.show()

    species_wtd_df = clustering_df[clustering_df["Resolution"] == "species"]
    line_plt = px.line(
        species_wtd_df.groupby(["Clustering", "RefPkg", "Length"]).mean(numeric_only=True).reset_index(),
        x="Length", y="WTD", color="Clustering",
        color_discrete_sequence=palette,
        facet_col_spacing=0.05,
        line_shape="spline",
        category_orders=_CATEGORIES,
        facet_col="RefPkg",
        labels=_LABEL_MAT,
        title="Effect of sequence length on cluster taxonomic distinctness at species resolution")
    line_plt.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    line_plt.update_traces(line=dict(width=2))
    line_plt.update_yaxes(autorange=True)
    line_plt.update_xaxes(tickangle=45,
                          title_font=None,
                          title_standoff=10)
    # line_plt.show()
    write_images_from_dict({rp_box_path: rp_box_plt, cm_box_path: cm_box_plt, line_path: line_plt})
    return


def acc_box(clustering_df: pd.DataFrame, palette) -> go.Figure:
    refpkg_acc_df = clustering_df.groupby(["RefPkg", "Clustering", "Length", "Resolution"]).mean().reset_index()
    box_plt = px.box(refpkg_acc_df,
                     x="Clustering", y="Accuracy", color="Clustering",
                     color_discrete_sequence=palette,
                     category_orders=_CATEGORIES,
                     labels=_LABEL_MAT,
                     title="Comparing cluster accuracy between reference-guided and de novo methods")
    for trace in box_plt['data']:
        trace['showlegend'] = False
    traces = refpkg_traces_for_plot(refpkg_acc_df)
    box_plt = go.Figure(({"data": tuple(list(box_plt["data"]) + traces),
                          "layout": box_plt["layout"]}))
    box_plt.update_layout(legend=dict(title="Reference Package"))
    # box_plt.show()
    return box_plt


def taxonomic_accuracy_plots(clustering_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.T10
    x_range = [20, 101]
    box_path = os.path.join(output_dir, "accuracy_boxes")
    line_path = os.path.join(output_dir, "accuracy_lines")
    bar_path = os.path.join(output_dir, "length_bars")

    bar_plt = len_bars(clustering_df, x_lims=x_range)

    line_plt = acc_line(clustering_df, palette, x_lims=x_range)

    box_plt = acc_box(clustering_df, palette)

    write_images_from_dict({line_path: line_plt, box_path: box_plt, bar_path: bar_plt})
    return


def taxonomic_summary_plots(taxa_df: pd.DataFrame, output_dir: str) -> None:
    plt_path = os.path.join(output_dir, "taxa_stack")
    taxa_df.index = pd.MultiIndex.from_frame(taxa_df, names=["refpkg", "taxon", "related", "rank"])
    count_df = pd.merge(left=taxa_df.count(level="refpkg").get("Related").reset_index(name="sum"),
                        right=taxa_df.groupby(["RefPkg", "Rank"]).count().get("Related").reset_index(name="count"),
                        how="inner", left_on="refpkg", right_on="RefPkg")
    count_df["Proportion"] = count_df["count"] / count_df["sum"]

    stacked_ranks_plt = px.bar(count_df,
                               x="RefPkg", y="Proportion",
                               color="Rank",
                               color_discrete_map=_RANK_PALETTE_MAP,
                               category_orders=_CATEGORIES,
                               labels=_LABEL_MAT,
                               title="Taxonomic relationships between reference package and query sequences")
    stacked_ranks_plt.update_traces(marker_line_color='rgb(105,105,105)',
                                    marker_line_width=1)
    # stacked_ranks_plt.show()
    write_images_from_dict({plt_path: stacked_ranks_plt})
    return


def evolutionary_summary_plots(evo_df: pd.DataFrame, output_dir: str) -> None:
    plt_path = os.path.join(output_dir, "evo_dist_scatter")
    ps_plt = px.scatter(evo_df,
                        x="Pendant", y="Distal",
                        color="RefPkg",
                        facet_col="Proper",
                        facet_row="Clustering",
                        color_discrete_map=_REFPKG_PALETTE_MAP,
                        category_orders=_CATEGORIES,
                        labels=_LABEL_MAT,
                        render_mode="svg",
                        title="Distribution evolutionary distances between query and reference sequences")
    ps_plt.update_traces(marker=dict(size=4,
                                     line=dict(width=0.5,
                                               color='DarkSlateGrey')),
                         selector=dict(mode='markers'))

    # Remove the redundant y-axis labels
    for axis in ps_plt.layout:
        if type(ps_plt.layout[axis]) == go.layout.YAxis and ps_plt.layout[axis].anchor != 'x3':
            ps_plt.layout[axis].title.text = ''

    # Now fix the facet_row labels by removing the variable name and rotating
    for annotation in ps_plt['layout']['annotations']:
        annotation['textangle'] = 0
        if annotation["text"].startswith("Clustering"):
            annotation["text"] = annotation["text"].split("=")[-1]

    # Customize the colour legend
    ps_plt.update_layout(legend_title_text='')
    ps_plt.update_layout(legend=dict(
        orientation="h",
        itemsizing="constant",
        yanchor="bottom",
        y=-0.25,
        xanchor="center",
        x=0.5
    ))

    write_images_from_dict({plt_path: ps_plt})
    return


def evaluate_clusters(project_path: str, n_examples=0, reso_ranks=None, **kwargs):
    cluster_experiments = []
    refpkg_map = {}
    if reso_ranks is None:
        reso_ranks = {"family", "genus", "species"}
    data_dir = os.path.join(project_path, "clustering_experiments") + os.sep
    refpkg_dir = os.path.join(data_dir, _REFPKG_DIR)
    fig_dir = os.path.join(project_path, "manuscript", "figures") + os.sep
    tab_dir = os.path.join(project_path, "manuscript", "tables") + os.sep
    for dir_path in [fig_dir, tab_dir]:
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    # Process the PhylOTU outputs
    dirs = glob.glob(data_dir + "length_*/phylotu_outputs/*")
    p_bar = tqdm(total=len(dirs),
                 ncols=100,
                 desc="Loading phylotu outputs")
    for phylotu_dir in dirs:
        if n_examples and p_bar.n >= n_examples:
            break
        phylotu_exp = ClusterExperiment(directory=phylotu_dir)
        if not phylotu_exp.test_files():
            p_bar.update()
            continue
        if not phylotu_exp.load_cluster_assignments():
            p_bar.update()
            continue
        if phylotu_exp.cluster_resolution not in reso_ranks:
            p_bar.update()
            continue

        phylotu_exp.set_precluster_mode(append=True)
        # Load the ClusterExperiment's associated reference package
        if phylotu_exp.pkg_name not in refpkg_map:
            phylotu_exp.load_ref_pkg(refpkg_dir)
            refpkg_map[phylotu_exp.pkg_name] = phylotu_exp.ref_pkg
        else:
            phylotu_exp.ref_pkg = refpkg_map[phylotu_exp.pkg_name]
        phylotu_exp.parse_seq_length()
        phylotu_exp.merge_query_clusters()
        cluster_experiments.append(phylotu_exp)
        p_bar.update()
    p_bar.close()

    print("Downloading lineage information")
    acc_taxon_map = retrieve_lineages(cluster_experiments)
    print("Mapping query sequences to taxa")
    map_queries_to_taxa(cluster_experiments, acc_taxon_map)
    filter_incomplete_lineages(cluster_experiments)
    filter_by_seq_lengths(cluster_experiments, min_perc=20)

    evo_dist_df = prepare_evodist_accuracy_dataframe(cluster_experiments)
    evo_summary_stats(evo_dist_df, tab_dir)
    evolutionary_summary_plots(evo_dist_df, fig_dir)

    taxonomic_distinctness_plots(prepare_taxonomic_distinctness_dataframe(cluster_experiments), fig_dir)

    # Percentage of sequences that were clustered together correctly at each length and rank
    acc_df = prepare_clustering_accuracy_dataframe(cluster_experiments)
    summary_df = acc_summary_stats(acc_df, kwargs=kwargs)

    # Cohesiveness of clusters for each sliced sequence at each length
    comp_df = prepare_clustering_cohesion_dataframe(cluster_experiments)
    completeness_summary_stats(comp_df)

    summary_df.to_csv(os.path.join(tab_dir,
                                   "acc_summary_" + time.strftime("%d-%m-%y_%H%M%S", time.localtime()) + ".csv"),
                      index=False,
                      mode='w')
    taxonomic_accuracy_plots(acc_df, fig_dir)

    taxonomic_relationships_plot(prepare_relationships_dataframe(cluster_experiments), fig_dir)

    taxonomic_summary_plots(prepare_taxonomic_summary_dataframe(cluster_experiments), fig_dir)

    sequence_cohesion_plots(comp_df, fig_dir)

    return


if __name__ == "__main__":
    root_dir = "/media/connor/Rufus/Gene_Centric_Guide/"
    if len(sys.argv) == 2:
        test = sys.argv[1]
    else:
        test = 0
    evaluate_clusters(root_dir, test, confidence=0.99, dist="interval")
