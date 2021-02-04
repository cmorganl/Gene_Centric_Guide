#!/usr/bin/env python3

import re
import os
import glob
import unittest

import plotly.express as px
import pandas as pd

from treesapp import refpkg


_REFPKG_DIR = "clustering_refpkgs"
# _MODE_NAMES = {"rg": "Reference-guided", "dn": "de novo"}
# _RANK_PREFIX_MAP = {'s': "Species", 'f': "Family", 'c': "Class"}
_ACC_TAXON_MAP = {}


class Tester(unittest.TestCase):
    def test_report_query_cohesiveness(self):
        mock_ce = ClusterExperiment("./")
        mock_ce.report_query_cohesiveness({"seq1": ['1', '2', '1', '2']})
        return


class ClusterExperiment:
    length_parse_re = re.compile(r"length_(\d+)$")
    qseq_slider_re = re.compile(r"_sliding:\d+-\d+$")

    def __init__(self, directory: str):
        """Create an instance representing a single treesapp phylotu output."""
        self.dir_path = directory
        self.pquery_assignments_file = os.path.join(self.dir_path, "final_outputs", "phylotu_pquery_assignments.tsv")
        self.seq_length = None
        self.cluster_resolution = None
        self.cluster_mode = None
        self.pkg_name = None
        self.ref_pkg = refpkg.ReferencePackage()
        self.cluster_assignments = {}

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
    def report_query_cohesiveness(cluster_assignments: dict) -> list:
        fragment_ratio = []
        for clusters in cluster_assignments.values():  # type: list
            fragment_ratio.append((len(clusters) / (len(clusters) * len(set(clusters)))))

        return fragment_ratio


def cohesiveness_plots(cluster_experiments: list) -> None:
    ref_pkgs = []
    resolutions = []
    lengths = []
    cluster_modes = []
    cohesivity = []
    for phylotu_exp in sorted(cluster_experiments, key=lambda x: int(x.seq_length)):  # type: ClusterExperiment
        exp_cohesivity = phylotu_exp.report_query_cohesiveness(phylotu_exp.cluster_assignments)
        cluster_modes += [phylotu_exp.cluster_mode] * len(exp_cohesivity)
        ref_pkgs += [phylotu_exp.pkg_name] * len(exp_cohesivity)
        resolutions += [phylotu_exp.cluster_resolution] * len(exp_cohesivity)
        lengths += [int(phylotu_exp.seq_length)] * len(exp_cohesivity)
        cohesivity += exp_cohesivity

    frag_df = pd.DataFrame(dict(RefPkg=ref_pkgs,
                                Cohesivity=cohesivity,
                                Clustering=cluster_modes,
                                Resolution=resolutions,
                                Length=lengths))

    line_plt = px.line(frag_df.groupby(["Clustering", "Length"]).mean().reset_index(),
                       x="Length", y="Cohesivity",
                       color="Clustering", line_group="Clustering",
                       title="Cluster cohesiveness as a function of query sequence length")
    line_plt.show()

    bar_plt = px.bar(frag_df.groupby(["Clustering", "RefPkg"]).mean().reset_index(),
                     x="RefPkg", y="Cohesivity",
                     barmode="group", color="Clustering",
                     title="Cluster cohesiveness as a function of query sequence length")
    bar_plt.show()

    violin_plt = px.violin(frag_df.groupby(["RefPkg", "Clustering", "Resolution", "Length"]).mean().reset_index(),
                           x="Clustering", y="Cohesivity", color="Clustering",
                           color_discrete_sequence=px.colors.qualitative.T10,
                           box=True, points="all", range_y=[0, 1.01],
                           labels={"Clustering": "Clustering method",
                                   "Cohesivity": "Mean Cluster Cohesivity Index"},
                           title="Cluster cohesiveness as a function of query sequence length")
    violin_plt.show()

    line_plt.write_image("cohesiveness_line.png", engine="kaleido")
    bar_plt.write_image("cohesiveness_bar.png", engine="kaleido")
    violin_plt.write_image("cohesiveness_violin.png", engine="kaleido")

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

    # Cohesiveness of clusters for each sliced sequence at each length
    cohesiveness_plots(cluster_experiments)

    # TODO: Percentage of sequences that were clustered together correctly at each length and rank
    return


if __name__ == "__main__":
    evaluate_clusters()