import plotly.express as px
from bokeh import palettes


_LABEL_MAT = {"Tax.dist": "Taxonomic distance",
              "TaxDist": "Taxonomic distance",
              "EvoDist": "Evolutionary distance",
              "MCC": "Matthews correlation coefficient",
              "Fasta.Length": "Query size (aa/bp)",
              "Mb": "Query size (M aa/bp)",
              "Memory.Max (kbytes)": "Maximum memory (kbytes)",
              "Memory.Max (Mb)": "Maximum memory (Mb)",
              "Length": "Protein length percentage",
              "Completeness": "Cluster completeness score",
              "Accuracy": "Cluster accuracy",
              "RefPkg": "Reference package",
              "Resolution": "Cluster resolution",
              "Clustering": "Clustering method",
              "Spline": "Cluster accuracy",
              "WTD": "Taxonomic Distinctness",
              "Size": "Number of queries"}

_RANKS = ['root', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
_REFPKGS = ["RecA", "RpoB", "PF01655", "NifH", "SoxY", "McrA"]
_METHODS = ["GraftM", "TreeSAPP-LCA", "TreeSAPP-aELW", "TreeSAPP-LM", "TreeSAPP", "TreeSAPP v0.6.8", "DIAMOND"]
_CATEGORIES = {"Clustering": ["local", "de_novo-aln", "de_novo-psc", "ref_guided"],
               "RefPkg": _REFPKGS,
               "Rank": _RANKS,
               "Resolution": _RANKS,
               "Software": _METHODS,
               "Method": _METHODS,
               "Molecule": ["nuc", "aa"]}

_RANK_PAL = palettes.linear_palette(px.colors.diverging.PuOr, len(_RANKS))
_REFPKG_PAL = px.colors.qualitative.Safe
_METHODS_PAL = px.colors.qualitative.Bold
_ASM_PAL = px.colors.qualitative.Antique

_RANK_PALETTE_MAP = {_RANKS[i]: _RANK_PAL[i] for i in range(0, len(_RANKS))}
_REFPKG_PALETTE_MAP = {_REFPKGS[i]: _REFPKG_PAL[i] for i in range(0, len(_REFPKGS))}
_METHODS_PALETTE_MAP = {_METHODS[i]: _METHODS_PAL[i] for i in range(0, len(_METHODS))}
