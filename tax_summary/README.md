# Taxonomic summary benchmark

The purpose of this is to determine the best taxonomic summary method for TreeSAPP,
and potentially other taxonomic classification tools using phylogenetic placements.

Here, we tested accumulated Likelihood Weight Ratio (aELW) and a linear model of placement distances with the placement with greatest LWR.

The scripts used to facilitate and compare these methods were 
[tax_summary.commands.sh](https://github.com/cmorganl/Gene_Centric_Guide/blob/main/src/tax_summary.commands.sh) and
[evaluate_tax_summary.py](https://github.com/cmorganl/Gene_Centric_Guide/blob/main/src/evaluate_tax_summary.py).

The data used for the evaluation were the Prokaryotic sequences in EggNOG v5.0.
