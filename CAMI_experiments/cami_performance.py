#!/usr/bin/env python

import os
import re

import plotly.io as pio
import plotly.express as px
import pandas as pd

from treesapp import file_parsers
from treesapp import classy
from treesapp import taxonomic_hierarchy as ts_th
from treesapp import entrez_utils as ts_entrez
from treesapp import phylo_seq as ts_phyloseq
from treesapp import refpkg as ts_refpkg
from treesapp import lca_calculations as ts_lca


_LABEL_MAT = {"Length": "Query sequence length",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method",
              "TaxDist": "Taxonomic distance"}


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
                pquery.ncbi_tax = seqname_taxid_map['_'.join(pquery.seq_name.split('_')[:-1])]
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

    for sample, refpkg_pqueries in sample_assigned_pqueries.items():
        for ref_pkg_name, pqueries in refpkg_pqueries.items():
            ref_pkg = refpkg_dict[ref_pkg_name]  # type: ts_refpkg.ReferencePackage
            ref_pkg.taxa_trie.build_multifurcating_trie()
            for pq in pqueries:  # type: ts_phyloseq.PQuery
                # Find the optimal taxonomic assignment
                true_lineage = ref_pkg.taxa_trie.clean_lineage_string(tax_lineage_map[pq.ncbi_tax])
                optimal_taxon = ts_lca.optimal_taxonomic_assignment(ref_pkg.taxa_trie.trie, true_lineage)
                if not optimal_taxon:
                    print("Optimal taxonomic assignment '{}' for {}"
                          " not found in reference hierarchy.\n".format(true_lineage,
                                                                        pq.place_name))
                    continue
                tax_dist, status = ts_lca.compute_taxonomic_distance(pq.lineage, optimal_taxon)
                if status > 0:
                    print("Lineages didn't converge between:\n"
                          "'{}' and '{}' (taxid: {})\n".format(pq.lineage, optimal_taxon, pq.ncbi_tax))
                samples_ar.append(sample)
                ref_pkgs_ar.append(ref_pkg_name)
                dists_ar.append(tax_dist)

    return pd.DataFrame(dict(RefPkg=ref_pkgs_ar, Sample=samples_ar, TaxDist=dists_ar))


def plot_taxonomic_distance_bubbles(tax_dist_dat: pd.DataFrame) -> None:
    palette = px.colors.qualitative.T10
    bubble_plt = px.scatter(tax_dist_dat.groupby(["RefPkg", "Sample"]).mean().reset_index(),
                            x="RefPkg", y="TaxDist", color="Sample",
                            color_discrete_sequence=palette,
                            labels=_LABEL_MAT,
                            title="")
    bubble_plt.show()
    bubble_plt.write_image("tax_dist_violin.png", engine="kaleido", scale=4.0)
    return


def main():
    assign_outputs = {"gold_standard_high_single": "gsa_mapping_pool.binning",
                      "RH_S001_merged": "gs_read_mapping_1.binning"}
    refpkg_dir = "refpkgs"
    refpkg_dict = ts_refpkg.gather_ref_packages(refpkg_data_dir=refpkg_dir)

    sample_assigned_pqueries = {}

    # Read the assignments from each of the assign_outputs directories
    for assign_out in assign_outputs:
        sample_assigned_pqueries.update({assign_out:
                                             file_parsers.load_classified_sequences_from_assign_output(assign_output_dir=assign_out,
                                                                                                       assigner_cls=classy.TreeSAPP("cami"))})

    # Fetch the taxonomic lineages for the unique NCBI taxonomy IDs from the classified query sequences
    tax_lineage_map = retrieve_lineages(sample_assigned_pqueries, assign_outputs)

    # Calculate the taxonomic distances by RefPkg and output, return a pandas dataframe
    tax_dist_dat = calculate_taxonomic_distances(sample_assigned_pqueries, tax_lineage_map, refpkg_dict)

    # Plot the taxonomic distances
    plot_taxonomic_distance_bubbles(tax_dist_dat)

    # TODO: Count the number of pOTUs for each RefPkg in each output

    # TODO: Plot the pOTUs

    return


if __name__ == "__main__":
    main()
