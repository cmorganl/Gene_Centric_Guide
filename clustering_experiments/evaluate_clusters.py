#!/usr/bin/env python3

import re
import os
import glob

import plotly.io as pio
import plotly.express as px
import pandas as pd
from sklearn.metrics import homogeneity_completeness_v_measure, completeness_score
from collections import Counter

from treesapp import refpkg
from treesapp import fasta as ts_fasta
from treesapp import classy as ts_classes
from treesapp import entrez_utils as ts_entrez
from treesapp import taxonomic_hierarchy as ts_tax


pio.templates.default = "plotly_white"
_REFPKG_DIR = "clustering_refpkgs"
_LABEL_MAT = {"Length": "Query sequence length",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method"}


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
        self.pkg_name = None
        self.ref_pkg = refpkg.ReferencePackage()
        self.cluster_assignments = {}
        self.entrez_query_dict = {}
        self.query_taxon_map = {}

        return

    def test_files(self) -> bool:
        if not os.path.exists(self.pquery_assignments_file):
            # raise AssertionError("No phylotu_pquery_assignments.tsv file found for '{}'.".format(self.dir_path))
            return False
        return True

    def parse_seq_length(self):
        for dir_name in self.dir_path.split(os.sep):
            if self.length_parse_re.match(dir_name):
                self.seq_length = self.length_parse_re.match(dir_name).group(1)
                break
            elif dir_name == "length_full":
                self.seq_length = "Full-length"
        if not self.seq_length:
            raise AssertionError("Unable to determine the sequence length used in '{}'.".format(self.dir_path))
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
        if header != ['PQuery', 'RefPkg', 'Resolution', 'Mode', 'OTU_ID']:
            raise AssertionError("Unexpected table format detected by header differences.")

        line = assignments_handler.readline()
        if not line:
            return 0
        while line:
            qseq_name, ref_pkg, resolution, mode, cluster_id = line.strip().split(delim)
            self.cluster_assignments[qseq_name] = cluster_id
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

    def load_ref_pkg(self):
        self.ref_pkg.f__json = os.path.join(_REFPKG_DIR), self.pkg_name + "_build.pkl"
        self.ref_pkg.slurp()
        return

    def merge_query_clusters(self):
        merged_cluster_assignments = {}
        for qseq_name, cluster in self.cluster_assignments.items():
            if self.qseq_slider_re.search(qseq_name):
                qseq_name = self.qseq_slider_re.sub(repl='', string=qseq_name)
            try:
                merged_cluster_assignments[qseq_name].append(cluster)
            except KeyError:
                merged_cluster_assignments[qseq_name] = [cluster]
        self.cluster_assignments = merged_cluster_assignments
        return

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
        header_registry = ts_fasta.register_headers(list(self.cluster_assignments.keys()))
        entrez_record_dict = ts_classes.get_header_info(header_registry)
        for index, e_record in entrez_record_dict.items():  # type: ts_entrez.EntrezRecord
            self.entrez_query_dict[e_record.description] = e_record
        return


def retrieve_lineages(cluster_experiments) -> dict:
    """
    Determines the format of the query sequence header to extract the accession and/or NCBI taxid then
    queries the Entrez database using these information to collect the taxonomic lineage for each unique NCBI taxid
    NCBI taxid and lineage information are stored in self.tax_lineage_map

    :return: None
    """
    # Gather the unique taxonomy IDs and store in EntrezRecord instances
    t_hierarchy = ts_tax.TaxonomicHierarchy()
    entrez_records = []
    acc_taxon_map = dict()
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        phylotu_exp.generate_entrez_queries()
        entrez_records += [phylotu_exp.entrez_query_dict[index] for index in phylotu_exp.entrez_query_dict]
    # Query the Entrez database for these unique taxonomy IDs
    ts_entrez.get_multiple_lineages(entrez_records, t_hierarchy, "prot")
    t_hierarchy.root_domains(t_hierarchy.find_root_taxon())

    for e_record in entrez_records:  # type: ts_entrez.EntrezRecord
        acc_taxon_map[e_record.ncbi_tax] = t_hierarchy.get_taxon(e_record.lineage.split(t_hierarchy.lin_sep)[-1])

    return acc_taxon_map


def map_queries_to_taxa(cluster_experiments: list, accession_taxon_map: dict) -> None:
    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        for query_name in phylotu_exp.cluster_assignments:
            query_acc = query_name.split('.')[0]
            phylotu_exp.query_taxon_map[query_name] = accession_taxon_map[query_acc]
    return


def get_key(a_dict: dict, val):
    for key, value in a_dict.items():
        if val == value:
            return key


def taxonomically_resolve_clusters(phylotu_exp: ClusterExperiment) -> (dict, dict):
    taxon_cluster_ids = {}
    re_map = {}
    for query in phylotu_exp.cluster_assignments:  # type: str
        taxon = phylotu_exp.query_taxon_map[query]  # type: ts_tax.Taxon
        taxon_resolved = taxon.get_rank_in_lineage(phylotu_exp.cluster_resolution)  # type: ts_tax.Taxon
        if not taxon_resolved:
            rank_depth = phylotu_exp.ref_pkg.taxa_trie.accepted_ranks_depths[phylotu_exp.cluster_resolution]
            parent_rank = get_key(phylotu_exp.ref_pkg.taxa_trie.accepted_ranks_depths, rank_depth-1)
            try:
                taxon_resolved = taxon.get_rank_in_lineage(parent_rank).name + " " + phylotu_exp.cluster_resolution
            except KeyError:
                print("Unable to find the {} rank of taxon '{}'.\n".format(phylotu_exp.cluster_resolution, taxon.name))
                continue
            re_map[taxon.name] = taxon_resolved
        else:
            taxon_resolved = taxon_resolved.name

        try:
            taxon_cluster_ids[taxon_resolved] += phylotu_exp.cluster_assignments[query]
        except KeyError:
            taxon_cluster_ids[taxon_resolved] = phylotu_exp.cluster_assignments[query]

    return taxon_cluster_ids, re_map


def prepare_clustering_cohesion_dataframe(cluster_experiments: list) -> pd.DataFrame:
    ref_pkgs = []
    resolutions = []
    lengths = []
    cluster_modes = []
    cohesion = []
    for phylotu_exp in sorted(cluster_experiments, key=lambda x: int(x.seq_length)):  # type: ClusterExperiment
        for seq_name in phylotu_exp.cluster_assignments:
            cohesion.append(phylotu_exp.report_query_completeness(phylotu_exp.cluster_assignments[seq_name]))
            cluster_modes.append(phylotu_exp.cluster_mode)
            ref_pkgs.append(phylotu_exp.pkg_name)
            resolutions.append(phylotu_exp.cluster_resolution)
            lengths.append(int(phylotu_exp.seq_length))

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

    for phylotu_exp in cluster_experiments:  # type: ClusterExperiment
        if phylotu_exp.cluster_resolution not in re_map:
            re_map[phylotu_exp.cluster_resolution] = {}
        taxon_cluster_ids, renamed_taxa = taxonomically_resolve_clusters(phylotu_exp)
        re_map[phylotu_exp.cluster_resolution].update(renamed_taxa)
        for taxon in taxon_cluster_ids:  # type: str
            refpkgs.append(phylotu_exp.pkg_name)
            lengths.append(int(phylotu_exp.seq_length))
            clusters.append(phylotu_exp.cluster_mode)
            resos.append(phylotu_exp.cluster_resolution)
            taxa.append(taxon)
            accurs.append(phylotu_exp.find_clustering_accuracy(taxon_cluster_ids[taxon]))

    for cluster_res in re_map:  # type: str
        if re_map[cluster_res]:  # type: dict
            for taxon_name, parent_name in re_map[cluster_res].items():
                print("Set missing taxonomic {} of {} to '{}'.".format(cluster_res, taxon_name, parent_name))

    return pd.DataFrame(dict(RefPkg=refpkgs, Length=lengths, Clustering=clusters,
                             Resolution=resos, Taxon=taxa, Accuracy=accurs))


def sequence_cohesion_plots(frag_df: pd.DataFrame) -> None:
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
                     labels=_LABEL_MAT,
                     title="Mean cluster completeness across the different reference packages")
    # bar_plt.show()

    violin_plt = px.violin(frag_df.groupby(["RefPkg", "Clustering", "Resolution", "Length"]).mean().reset_index(),
                           x="Clustering", y="Completeness", color="Clustering",
                           color_discrete_sequence=palette,
                           box=True, points="all", range_y=[0, 1.01],
                           labels=_LABEL_MAT,
                           title="Comparing cluster completeness between reference-guided and de novo methods")
    # violin_plt.show()

    line_plt.write_image("completeness_line.png", engine="kaleido", scale=4.0)
    bar_plt.write_image("completeness_bar.png", engine="kaleido", scale=4.0)
    violin_plt.write_image("completeness_violin.png", engine="kaleido", scale=4.0)

    return


def taxonomic_accuracy_plots(clustering_df: pd.DataFrame) -> None:
    palette = px.colors.qualitative.T10
    bar_plt = px.bar(clustering_df.groupby(["Resolution", "Length", "Clustering"]).mean().reset_index(),
                     x="Length", y="Accuracy", color="Clustering",
                     color_discrete_sequence=palette, barmode="group",
                     facet_col="Resolution",
                     labels=_LABEL_MAT,
                     title="Comparing clustering accuracy between reference-guided and de novo methods")
    # bar_plt.show()

    violin_plt = px.violin(clustering_df.groupby(["RefPkg", "Clustering", "Length", "Resolution"]).mean().reset_index(),
                           x="Clustering", y="Accuracy", color="Clustering",
                           color_discrete_sequence=palette,
                           box=True, points="all", range_y=[0, 1.01],
                           labels=_LABEL_MAT,
                           title="Comparing cluster accuracy between reference-guided and de novo methods")
    # violin_plt.show()

    bar_plt.write_image("accuracy_bars.png", engine="kaleido", scale=4.0)
    violin_plt.write_image("accuracy_violin.png", engine="kaleido", scale=4.0)
    return


def evaluate_clusters():
    cluster_experiments = []
    # Process the PhylOTU outputs
    for phylotu_dir in glob.glob("length_*/phylotu_outputs/*"):
        phylotu_exp = ClusterExperiment(directory=phylotu_dir)
        if not phylotu_exp.test_files():
            continue
        phylotu_exp.parse_seq_length()
        if not phylotu_exp.load_cluster_assignments():
            continue
        phylotu_exp.merge_query_clusters()
        cluster_experiments.append(phylotu_exp)

    acc_taxon_map = retrieve_lineages(cluster_experiments)
    map_queries_to_taxa(cluster_experiments, acc_taxon_map)

    # Cohesiveness of clusters for each sliced sequence at each length
    sequence_cohesion_plots(prepare_clustering_cohesion_dataframe(cluster_experiments))

    # Percentage of sequences that were clustered together correctly at each length and rank
    taxonomic_accuracy_plots(prepare_clustering_accuracy_dataframe(cluster_experiments))
    return


if __name__ == "__main__":
    evaluate_clusters()
