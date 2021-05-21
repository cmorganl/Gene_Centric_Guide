#!/usr/bin/env python3

import glob
import os
import sys
import re
from collections import Counter

import pandas as pd
import treesapp as ts
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from bokeh import palettes
from tqdm import tqdm
from sklearn.metrics import homogeneity_completeness_v_measure, completeness_score

pio.templates.default = "plotly_white"
_REFPKG_DIR = "clustering_refpkgs"
_LABEL_MAT = {"Length": "Protein length percentage",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method"}
_RANKS = ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
_CATEGORIES = {"Clustering": ["de_novo-aln", "de_novo-psc", "ref_guided"],
               "RefPkg": ["RecA", "RpoB", "PF01655",
                          "NifH", "SoxY", "McrA"],
               "Rank": _RANKS}
_RANK_PAL = palettes.linear_palette(px.colors.diverging.PuOr, len(_RANKS))
_RANK_PALETTE_MAP = {_RANKS[i]: _RANK_PAL[i] for i in range(0, len(_RANKS))}


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
        return round(completeness_score(labels_true=true_labels, labels_pred=cluster_assignments), 3)

    @staticmethod
    def find_clustering_accuracy(taxon_cluster_ids: list) -> float:
        true_labels = ClusterExperiment.generate_true_cluster_labels(taxon_cluster_ids)
        h, c, v = homogeneity_completeness_v_measure(labels_true=true_labels, labels_pred=taxon_cluster_ids)
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
        hmm_perc = float(100*int(phylotu_exp.seq_length)/phylotu_exp.ref_pkg.hmm_length())
        for seq_name in phylotu_exp.cluster_assignments:
            cohesion.append(phylotu_exp.report_query_completeness(phylotu_exp.cluster_assignments[seq_name]))
            cluster_modes.append(phylotu_exp.cluster_mode)
            ref_pkgs.append(phylotu_exp.pkg_name)
            resolutions.append(phylotu_exp.cluster_resolution)
            lengths.append(hmm_perc)
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
        hmm_perc = float(100*int(phylotu_exp.seq_length)/phylotu_exp.ref_pkg.hmm_length())
        for taxon in taxon_cluster_ids:  # type: str
            refpkgs.append(phylotu_exp.pkg_name)
            lengths.append(hmm_perc)
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
                      marker_line_width=1.5)

    fig.write_image(os.path.join(output_dir, "relationship_bars.png"), engine="kaleido", scale=4.0)
    return


def prepare_evodist_accuracy_dataframe(cluster_experiments: list) -> pd.DataFrame:
    data_dict = {"RefPkg": [],
                 "Clustering": [],
                 "Length": [],
                 "Query": [],
                 "Proper": [],
                 "Pendant": [],
                 "Distal": []}
    p_bar = tqdm(total=len(cluster_experiments), ncols=100, desc="Preparing evolutionary distance dataframe")
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        taxon_cluster_ids = phylotu_exp.taxonomically_resolve_clusters()
        true_clusters = {taxon: phylotu_exp.generate_true_cluster_labels(taxon_cluster_ids[taxon]) for
                         taxon in taxon_cluster_ids}
        for query_name, fragments in phylotu_exp.cluster_assignments.items():  # type: (str, list)
            repr_taxon = phylotu_exp.match_query_to_taxon_cluster(query_name)
            if not repr_taxon:
                continue
            data_dict["RefPkg"] += [phylotu_exp.pkg_name]*len(fragments)
            data_dict["Clustering"] += [phylotu_exp.cluster_mode]*len(fragments)
            data_dict["Length"] += [phylotu_exp.seq_length]*len(fragments)
            data_dict["Query"] += [query_name]*len(fragments)
            for potu in fragments:
                if potu in true_clusters[repr_taxon]:
                    data_dict["Proper"].append(True)
                else:
                    data_dict["Proper"].append(False)
            for pplace in phylotu_exp.phylo_place_map[query_name]:  # type: ts.phylo_seq.PhyloPlace
                data_dict["Pendant"].append(pplace.pendant_length)
                data_dict["Distal"].append(pplace.distal_length)
        p_bar.update()
    p_bar.close()
    evo_df = pd.DataFrame(data_dict)
    return evo_df.groupby(["Clustering", "RefPkg", "Query", "Proper"]).mean(numeric_only=True).reset_index()


def sequence_cohesion_plots(frag_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.T10
    line_plt = px.line(frag_df.groupby(["Clustering", "Length"]).mean().reset_index(),
                       x="Length", y="Completeness",
                       color="Clustering", line_group="Clustering",
                       color_discrete_sequence=palette,
                       labels=_LABEL_MAT,
                       title="Cluster completeness as a function of query sequence length")
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

    line_plt.write_image(os.path.join(output_dir, "completeness_line.png"), engine="kaleido", scale=4.0)
    bar_plt.write_image(os.path.join(output_dir, "completeness_bar.png"), engine="kaleido", scale=4.0)
    violin_plt.write_image(os.path.join(output_dir, "completeness_violin.png"), engine="kaleido", scale=4.0)

    return


def taxonomic_accuracy_plots(clustering_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.T10
    acc_line_plt = px.line(clustering_df.groupby(["Resolution", "Length", "Clustering"]).mean().reset_index(),
                           x="Length", y="Accuracy", color="Clustering",
                           color_discrete_sequence=palette, line_group="Clustering",
                           facet_col_spacing=0.05,
                           range_x=[10, 105],
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

    violin_plt = px.violin(clustering_df.groupby(["RefPkg", "Clustering", "Length", "Resolution"]).mean().reset_index(),
                           x="Clustering", y="Accuracy", color="Clustering",
                           color_discrete_sequence=palette,
                           box=True, points="all", range_y=[0, 1.01],
                           labels=_LABEL_MAT,
                           title="Comparing cluster accuracy between reference-guided and de novo methods")
    violin_plt.update_layout(showlegend=False)
    # violin_plt.show()

    acc_line_plt.write_image(os.path.join(output_dir, "accuracy_lines.png"), engine="kaleido", scale=4.0)
    # acc_line_plt.write_image(os.path.join(output_dir, "accuracy_lines.svg"), engine="kaleido", scale=4.0)
    violin_plt.write_image(os.path.join(output_dir, "accuracy_violin.png"), engine="kaleido", scale=4.0)
    return


def taxonomic_summary_plots(taxa_df: pd.DataFrame, output_dir: str) -> None:
    taxa_df.index = pd.MultiIndex.from_frame(taxa_df, names=["refpkg", "taxon", "related", "rank"])
    count_df = pd.merge(left=taxa_df.count(level="refpkg").get("Related").reset_index(name="sum"),
                        right=taxa_df.groupby(["RefPkg", "Rank"]).count().get("Related").reset_index(name="count"),
                        how="inner", left_on="refpkg", right_on="RefPkg")
    count_df["Proportion"] = count_df["count"]/count_df["sum"]

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
    stacked_ranks_plt.write_image(os.path.join(output_dir, "taxa_stack.png"), engine="kaleido", scale=4.0)
    return


def evolutionary_summary_plots(evo_df: pd.DataFrame, output_dir: str) -> None:
    palette = px.colors.qualitative.Safe
    ps_plt = px.scatter(evo_df,
                        x="Pendant", y="Distal",
                        color="RefPkg",
                        facet_col="Proper",
                        facet_row="Clustering",
                        color_discrete_sequence=palette,
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

    ps_plt.write_image(os.path.join(output_dir, "evo_dist_scatter.png"), engine="kaleido", scale=4.0)
    return


def evaluate_clusters(project_path: str, n_examples=0):
    cluster_experiments = []
    refpkg_map = {}
    data_dir = os.path.join(project_path, "clustering_experiments") + os.sep
    refpkg_dir = os.path.join(data_dir, _REFPKG_DIR)
    fig_dir = os.path.join(project_path, "manuscript", "figures") + os.sep
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
        # Load the ClusterExperiment's associated reference package
        phylotu_exp.set_precluster_mode(append=True)
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

    taxonomic_relationships_plot(prepare_relationships_dataframe(cluster_experiments), fig_dir)

    # Taxonomic summary plots
    taxonomic_summary_plots(prepare_taxonomic_summary_dataframe(cluster_experiments), fig_dir)

    evolutionary_summary_plots(prepare_evodist_accuracy_dataframe(cluster_experiments), fig_dir)

    # Cohesiveness of clusters for each sliced sequence at each length
    sequence_cohesion_plots(prepare_clustering_cohesion_dataframe(cluster_experiments), fig_dir)

    # Percentage of sequences that were clustered together correctly at each length and rank
    taxonomic_accuracy_plots(prepare_clustering_accuracy_dataframe(cluster_experiments), fig_dir)
    return


if __name__ == "__main__":
    root_dir = "/media/connor/Rufus/Gene_Centric_Guide/"
    if len(sys.argv) == 2:
        test = sys.argv[1]
    else:
        test = 0
    evaluate_clusters(root_dir, test)

