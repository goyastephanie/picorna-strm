#!/usr/bin/env python3
"""
Rhinovirus amplicon contamination-network analysis.

Inputs
------
1) <run_prefix>_summary.tsv
   One row per sample/reference hit with abundance, coverage, depth and
   consensus-quality metrics. The run prefix can be any user-assigned name
   (for example: run, 20260806, QiaCap_49-104).

2) <run_prefix>_contamination.tsv
   Pairwise consensus comparisons with flags:
      - possible_contamination
      - possible_same_genome
      - ok

Core idea
---------
1. Collapse "possible_same_genome" reference hits ONLY WITHIN EACH SAMPLE.
   These are treated as alternative assemblies/reference assignments of the
   same biological genome, not as coinfections.

2. Build a graph whose nodes are collapsed sample-genome observations and
   whose edges are "possible_contamination" relationships.

3. For each contamination component, infer candidate source(s) from:
      - within-sample abundance
      - source:recipient abundance ratio
      - consensus identity
      - consensus quality
   The pipeline is deliberately conservative: ambiguous high-abundance nodes
   are retained rather than automatically removed.

4. Produce a cleanup recommendation:
      KEEP_SOURCE
      REMOVE_CONTAMINANT
      KEEP_UNRESOLVED
      KEEP_NO_CONTAMINATION_LINK

The pipeline does not modify FASTA files. The REMOVE_CONTAMINANT rows are a
decision table that can be used by a downstream cleanup step.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import math
import pandas as pd
import numpy as np


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def safe_num(x):
    return pd.to_numeric(x, errors="coerce")


def choose_metric_columns(summary: pd.DataFrame, metric: str):
    """
    Return abundance/depth/coverage columns.

    'c1' is useful when pct_reads_mapped_c1 represents the final/competitive
    mapping result. 'ref' uses the per-reference mapping columns.
    """
    if metric == "c1":
        cols = ("pct_reads_mapped_c1", "mean_depth_c1", "coverage_c1")
    else:
        cols = ("pct_reads_mapped_ref", "mean_depth_ref", "coverage_ref")
    missing = [x for x in cols if x not in summary.columns]
    if missing:
        raise ValueError(f"Missing required columns for metric={metric}: {missing}")
    return cols


def collapse_same_genome(summary, comparisons, metric="c1"):
    """
    Collapse possible_same_genome hits WITHIN each sample only.

    Returns
    -------
    nodes : one row per collapsed sample-genome observation
    original_to_node : dict[(sample, tag)] -> node_id
    """
    abundance_col, depth_col, coverage_col = choose_metric_columns(summary, metric)

    # Every sample/tag represented in the summary is a node before collapsing.
    keys = [(r.sample_name, r.ref_tag) for r in summary.itertuples()]
    uf = UnionFind(keys)

    same = comparisons[comparisons["flag"] == "possible_same_genome"].copy()
    for r in same.itertuples():
        # Critical safeguard: collapse only if both hits belong to the same sample.
        if r.sample_a != r.sample_b:
            continue
        a = (r.sample_a, r.tag_a)
        b = (r.sample_b, r.tag_b)
        if a in uf.parent and b in uf.parent:
            uf.union(a, b)

    groups = defaultdict(list)
    for k in keys:
        groups[uf.find(k)].append(k)

    rows = []
    original_to_node = {}

    # Keep one representative row: the hit with the largest abundance metric.
    summary_idx = summary.set_index(["sample_name", "ref_tag"], drop=False)

    for i, members in enumerate(sorted(groups.values(), key=lambda z: (z[0][0], str(z)))):
        sample = members[0][0]
        member_rows = []
        for key in members:
            x = summary_idx.loc[key]
            if isinstance(x, pd.DataFrame):
                member_rows.extend([row for _, row in x.iterrows()])
            else:
                member_rows.append(x)
        g = pd.DataFrame(member_rows)

        g["_abundance"] = safe_num(g[abundance_col])
        g["_depth"] = safe_num(g[depth_col])
        g["_coverage"] = safe_num(g[coverage_col])
        g["_pct_ns"] = safe_num(g["pct_ns"]) if "pct_ns" in g else np.nan

        # Select representative by abundance, then depth, then lower N%.
        g = g.sort_values(
            ["_abundance", "_depth", "_pct_ns"],
            ascending=[False, False, True],
            na_position="last",
        )
        rep = g.iloc[0]

        node_id = f"N{i+1:04d}"
        tags = sorted(set(str(x[1]) for x in members))
        refs = sorted(set(str(x) for x in g["ref_acc"].dropna())) if "ref_acc" in g else []

        row = {
            "node_id": node_id,
            "sample_name": sample,
            "representative_tag": rep["ref_tag"],
            "collapsed_tags": ",".join(tags),
            "n_collapsed_hits": len(tags),
            "abundance_pct": rep["_abundance"],
            "mean_depth": rep["_depth"],
            "coverage_pct": rep["_coverage"],
            "pct_ns": rep["_pct_ns"],
            "consensus_length": safe_num(rep.get("consensus_length", np.nan)),
            "num_ambiguous": safe_num(rep.get("num_ambiguous", np.nan)),
            "representative_ref_acc": rep.get("ref_acc", ""),
            "all_ref_acc": ",".join(refs),
        }
        rows.append(row)
        for key in members:
            original_to_node[key] = node_id

    nodes = pd.DataFrame(rows)

    # Within-sample abundance rank/context.
    nodes["sample_max_abundance_pct"] = nodes.groupby("sample_name")["abundance_pct"].transform("max")
    nodes["fraction_of_sample_max"] = nodes["abundance_pct"] / nodes["sample_max_abundance_pct"].replace(0, np.nan)
    nodes["abundance_rank_in_sample"] = (
        nodes.groupby("sample_name")["abundance_pct"]
        .rank(method="min", ascending=False)
        .astype("Int64")
    )
    return nodes, original_to_node


def build_contamination_edges(comparisons, original_to_node):
    cont = comparisons[comparisons["flag"] == "possible_contamination"].copy()
    rows = []

    for r in cont.itertuples():
        a = original_to_node.get((r.sample_a, r.tag_a))
        b = original_to_node.get((r.sample_b, r.tag_b))
        if a is None or b is None or a == b:
            continue
        # undirected edge: source direction is inferred later
        x, y = sorted([a, b])
        rows.append({
            "node_a": x,
            "node_b": y,
            "sample_a_original": r.sample_a,
            "tag_a_original": r.tag_a,
            "sample_b_original": r.sample_b,
            "tag_b_original": r.tag_b,
            "comparable_sites": safe_num(r.comparable_sites),
            "differences": safe_num(r.differences),
            "pct_identity": safe_num(r.pct_identity),
        })

    if not rows:
        return pd.DataFrame(columns=[
            "node_a","node_b","max_pct_identity","min_differences",
            "max_comparable_sites","n_supporting_comparisons"
        ])

    x = pd.DataFrame(rows)

    # Multiple raw comparisons can map to the same collapsed node pair because
    # of alternative reference labels. Summarize them into a single edge.
    edges = (
        x.groupby(["node_a", "node_b"], as_index=False)
        .agg(
            max_pct_identity=("pct_identity", "max"),
            min_differences=("differences", "min"),
            max_comparable_sites=("comparable_sites", "max"),
            n_supporting_comparisons=("pct_identity", "size"),
            raw_tag_pairs=("tag_a_original",
                           lambda _: "")  # filled below
        )
    )

    pair_labels = (
        x.assign(pair=x["tag_a_original"].astype(str) + "~" + x["tag_b_original"].astype(str))
         .groupby(["node_a","node_b"])["pair"]
         .apply(lambda z: ",".join(sorted(set(z))))
         .reset_index(name="raw_tag_pairs")
    )
    edges = edges.drop(columns=["raw_tag_pairs"]).merge(pair_labels, on=["node_a","node_b"], how="left")
    return edges


def connected_components(node_ids, edges):
    adj = {n: set() for n in node_ids}
    for r in edges.itertuples():
        adj[r.node_a].add(r.node_b)
        adj[r.node_b].add(r.node_a)

    seen = set()
    components = []
    for n in node_ids:
        if n in seen:
            continue
        stack = [n]
        comp = []
        seen.add(n)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        components.append(sorted(comp))
    return components, adj



def shortest_path(adj, start, target):
    """Return an unweighted shortest path between two nodes, or None."""
    if start == target:
        return [start]
    queue = [(start, [start])]
    seen = {start}
    while queue:
        u, path = queue.pop(0)
        for v in adj[u]:
            if v in seen:
                continue
            new_path = path + [v]
            if v == target:
                return new_path
            seen.add(v)
            queue.append((v, new_path))
    return None


def infer_network(
    nodes,
    edges,
    source_min_pct=20.0,
    contaminant_max_pct=10.0,
    source_ratio=5.0,
    high_conf_ratio=10.0,
    high_conf_identity=99.95,
):
    """
    Conservative source/recipient inference.

    A node can be removed only if:
      1) it has a direct possible_contamination edge to a plausible source;
      2) its abundance <= contaminant_max_pct;
      3) the candidate source abundance >= source_min_pct;
      4) source abundance / recipient abundance >= source_ratio.

    High-abundance observations are never removed automatically merely because
    another sample has a nearly identical genome.
    """
    nodes = nodes.copy()
    node_lookup = nodes.set_index("node_id").to_dict("index")
    components, adj = connected_components(nodes["node_id"].tolist(), edges)

    edge_lookup = {}
    for r in edges.itertuples():
        edge_lookup[(r.node_a, r.node_b)] = r
        edge_lookup[(r.node_b, r.node_a)] = r

    component_id = {}
    for i, comp in enumerate(components, start=1):
        for n in comp:
            component_id[n] = f"C{i:03d}"

    decisions = []
    directed_edges = []

    for n in nodes["node_id"]:
        info = node_lookup[n]
        neighbors = sorted(adj[n])

        if not neighbors:
            decisions.append({
                "node_id": n,
                "component_id": component_id[n],
                "classification": "no_contamination_link",
                "cleanup_action": "KEEP_NO_CONTAMINATION_LINK",
                "candidate_source_node": "",
                "candidate_source_sample": "",
                "source_abundance_pct": np.nan,
                "recipient_abundance_pct": info["abundance_pct"],
                "source_to_recipient_ratio": np.nan,
                "edge_pct_identity": np.nan,
                "confidence": "NA",
                "inference_mode": "none",
                "path_length_to_source": np.nan,
                "reason": "No possible_contamination edge after same-genome collapsing.",
            })
            continue

        recipient_ab = info["abundance_pct"]
        candidates = []
        for src in neighbors:
            sinfo = node_lookup[src]
            src_ab = sinfo["abundance_pct"]
            ratio = np.inf if (pd.notna(recipient_ab) and recipient_ab == 0) else (
                src_ab / recipient_ab if pd.notna(src_ab) and pd.notna(recipient_ab) and recipient_ab > 0 else np.nan
            )
            e = edge_lookup[(n, src)]
            if (
                pd.notna(recipient_ab)
                and recipient_ab <= contaminant_max_pct
                and pd.notna(src_ab)
                and src_ab >= source_min_pct
                and pd.notna(ratio)
                and ratio >= source_ratio
            ):
                candidates.append((src, src_ab, ratio, e.max_pct_identity, e.min_differences))

        if candidates:
            # Best source: ratio, abundance, identity.
            candidates.sort(key=lambda z: (
                -np.inf if pd.isna(z[2]) else z[2],
                -np.inf if pd.isna(z[1]) else z[1],
                -np.inf if pd.isna(z[3]) else z[3],
            ), reverse=True)
            src, src_ab, ratio, ident, diffs = candidates[0]
            sinfo = node_lookup[src]

            if (
                ratio >= high_conf_ratio
                and pd.notna(ident)
                and ident >= high_conf_identity
            ):
                conf = "HIGH"
            else:
                conf = "MODERATE"

            decisions.append({
                "node_id": n,
                "component_id": component_id[n],
                "classification": "likely_contaminant",
                "cleanup_action": "REMOVE_CONTAMINANT",
                "candidate_source_node": src,
                "candidate_source_sample": sinfo["sample_name"],
                "source_abundance_pct": src_ab,
                "recipient_abundance_pct": recipient_ab,
                "source_to_recipient_ratio": ratio,
                "edge_pct_identity": ident,
                "confidence": conf,
                "inference_mode": "direct",
                "path_length_to_source": 1,
                "reason": (
                    f"Low-abundance node has a directly linked high-abundance candidate source "
                    f"(ratio={ratio:.2f}x; identity={ident:.4f}%)."
                ),
            })
            directed_edges.append({
                "source_node": src,
                "recipient_node": n,
                "source_sample": sinfo["sample_name"],
                "recipient_sample": info["sample_name"],
                "source_tag": sinfo["representative_tag"],
                "recipient_tag": info["representative_tag"],
                "source_abundance_pct": src_ab,
                "recipient_abundance_pct": recipient_ab,
                "abundance_ratio": ratio,
                "pct_identity": ident,
                "differences": diffs,
                "confidence": conf,
                "inference_mode": "direct",
                "path_length": 1,
                "network_path": f"{src}>{n}",
            })
        else:
            # Keep high/intermediate abundance or otherwise unresolved nodes.
            plausible_recipient_neighbors = []
            for other in neighbors:
                oinfo = node_lookup[other]
                if pd.notna(oinfo["abundance_pct"]) and oinfo["abundance_pct"] <= contaminant_max_pct:
                    plausible_recipient_neighbors.append(other)

            if pd.notna(recipient_ab) and recipient_ab >= source_min_pct and plausible_recipient_neighbors:
                cls = "likely_source"
                action = "KEEP_SOURCE"
                conf = "HIGH" if recipient_ab >= 50 else "MODERATE"
                reason = (
                    "High-abundance observation linked to one or more low-abundance observations; "
                    "retain as a candidate native/source genome."
                )
            else:
                cls = "unresolved"
                action = "KEEP_UNRESOLVED"
                conf = "REVIEW"
                reason = (
                    "Contamination relationship exists, but abundance evidence is insufficient "
                    "to infer a safe source->recipient direction. Retain pending review."
                )

            decisions.append({
                "node_id": n,
                "component_id": component_id[n],
                "classification": cls,
                "cleanup_action": action,
                "candidate_source_node": "",
                "candidate_source_sample": "",
                "source_abundance_pct": np.nan,
                "recipient_abundance_pct": recipient_ab,
                "source_to_recipient_ratio": np.nan,
                "edge_pct_identity": np.nan,
                "confidence": conf,
                "inference_mode": "unresolved",
                "path_length_to_source": np.nan,
                "reason": reason,
            })

    decisions = pd.DataFrame(decisions)

    # ------------------------------------------------------------------
    # Second-pass root-source inference.
    #
    # A low-abundance unresolved node can be assigned to a component-level
    # source even without a direct source edge when the contamination network
    # has a UNIQUE high-abundance root. This captures paths such as:
    #
    #     true source -> contaminated intermediate -> contaminated sample
    #
    # Indirect calls are kept separate from direct calls and receive lower
    # confidence. If a component has multiple plausible high-abundance roots,
    # no indirect direction is imposed.
    # ------------------------------------------------------------------
    node_to_component = {n: component_id[n] for n in nodes["node_id"]}
    for cid in sorted(set(node_to_component.values())):
        comp_nodes = [n for n, c in node_to_component.items() if c == cid]
        high_nodes = [
            n for n in comp_nodes
            if pd.notna(node_lookup[n]["abundance_pct"])
            and node_lookup[n]["abundance_pct"] >= source_min_pct
        ]

        # Only propagate a root when the component has a unique high-abundance
        # candidate. Multiple high-abundance genomes are preserved as ambiguous.
        if len(high_nodes) != 1:
            continue

        root = high_nodes[0]
        root_info = node_lookup[root]
        root_ab = root_info["abundance_pct"]

        for idx, d in decisions.iterrows():
            n = d["node_id"]
            if node_to_component[n] != cid:
                continue
            if d["cleanup_action"] != "KEEP_UNRESOLVED":
                continue

            rec_ab = node_lookup[n]["abundance_pct"]
            if pd.isna(rec_ab) or rec_ab > contaminant_max_pct:
                continue
            if rec_ab <= 0:
                ratio = np.inf
            else:
                ratio = root_ab / rec_ab
            if ratio < source_ratio:
                continue

            path = shortest_path(adj, root, n)
            if not path or len(path) < 2:
                continue

            # Conservative indirect call: network connectivity plus a unique
            # dominant root and a strong abundance asymmetry.
            decisions.at[idx, "classification"] = "likely_contaminant"
            decisions.at[idx, "cleanup_action"] = "REMOVE_CONTAMINANT"
            decisions.at[idx, "candidate_source_node"] = root
            decisions.at[idx, "candidate_source_sample"] = root_info["sample_name"]
            decisions.at[idx, "source_abundance_pct"] = root_ab
            decisions.at[idx, "recipient_abundance_pct"] = rec_ab
            decisions.at[idx, "source_to_recipient_ratio"] = ratio
            decisions.at[idx, "edge_pct_identity"] = np.nan
            decisions.at[idx, "confidence"] = "MODERATE"
            decisions.at[idx, "inference_mode"] = "indirect_component_root"
            decisions.at[idx, "path_length_to_source"] = len(path) - 1
            decisions.at[idx, "reason"] = (
                f"Low-abundance node belongs to a contamination component with a unique "
                f"high-abundance root source (ratio={ratio:.2f}x); source is reached through "
                f"{len(path)-1} contamination edge(s)."
            )

            directed_edges.append({
                "source_node": root,
                "recipient_node": n,
                "source_sample": root_info["sample_name"],
                "recipient_sample": node_lookup[n]["sample_name"],
                "source_tag": root_info["representative_tag"],
                "recipient_tag": node_lookup[n]["representative_tag"],
                "source_abundance_pct": root_ab,
                "recipient_abundance_pct": rec_ab,
                "abundance_ratio": ratio,
                "pct_identity": np.nan,
                "differences": np.nan,
                "confidence": "MODERATE",
                "inference_mode": "indirect_component_root",
                "path_length": len(path) - 1,
                "network_path": ">".join(path),
            })

    # A node identified as source for any recipient should always be explicitly retained.
    source_nodes = set(decisions.loc[decisions["cleanup_action"] == "REMOVE_CONTAMINANT", "candidate_source_node"])
    mask = decisions["node_id"].isin(source_nodes) & (decisions["cleanup_action"] != "REMOVE_CONTAMINANT")
    decisions.loc[mask, "classification"] = "likely_source"
    decisions.loc[mask, "cleanup_action"] = "KEEP_SOURCE"
    decisions.loc[mask, "inference_mode"] = "source"
    decisions.loc[mask, "path_length_to_source"] = 0
    decisions.loc[mask, "reason"] = (
        "Selected as the high-abundance source for at least one contaminant."
    )

    directed = pd.DataFrame(directed_edges)
    return decisions, directed, component_id


def component_summary(nodes, decisions, edges):
    x = nodes.merge(
        decisions[["node_id","component_id","classification","cleanup_action"]],
        on="node_id", how="left"
    )
    rows = []
    for cid, g in x.groupby("component_id"):
        ids = set(g["node_id"])
        ge = edges[edges["node_a"].isin(ids) & edges["node_b"].isin(ids)]
        rows.append({
            "component_id": cid,
            "n_nodes": len(g),
            "n_samples": g["sample_name"].nunique(),
            "samples": ",".join(sorted(g["sample_name"].unique())),
            "representative_tags": ",".join(sorted(set(g["representative_tag"].astype(str)))),
            "max_abundance_pct": g["abundance_pct"].max(),
            "n_keep_source": (g["cleanup_action"] == "KEEP_SOURCE").sum(),
            "n_remove_contaminant": (g["cleanup_action"] == "REMOVE_CONTAMINANT").sum(),
            "n_keep_unresolved": (g["cleanup_action"] == "KEEP_UNRESOLVED").sum(),
            "n_edges": len(ge),
            "max_edge_identity": ge["max_pct_identity"].max() if len(ge) else np.nan,
        })
    return pd.DataFrame(rows).sort_values(
        ["n_remove_contaminant","n_nodes"], ascending=[False,False]
    )


def make_cleanup_table(nodes, decisions):
    x = nodes.merge(decisions, on="node_id", how="left")
    cols = [
        "sample_name","node_id","representative_tag","collapsed_tags",
        "n_collapsed_hits","abundance_pct","mean_depth","coverage_pct","pct_ns",
        "abundance_rank_in_sample","fraction_of_sample_max",
        "component_id","classification","cleanup_action","confidence","inference_mode","path_length_to_source",
        "candidate_source_sample","candidate_source_node",
        "source_abundance_pct","source_to_recipient_ratio","edge_pct_identity",
        "reason"
    ]
    return x[cols].sort_values(
        ["sample_name","cleanup_action","abundance_pct"],
        ascending=[True,True,False]
    )



def resolve_input_files(summary_arg=None, comparisons_arg=None, input_dir=".", run_prefix=None):
    """Resolve summary/comparison input paths while supporting arbitrary run prefixes.

    Backward-compatible modes
    -------------------------
    1) Explicit paths:
         --summary FILE --comparisons FILE

    2) Prefix-based discovery:
         --input-dir DIR --run-prefix PREFIX
       resolves:
         DIR/PREFIX_summary.tsv
         DIR/PREFIX_contamination.tsv
         DIR/PREFIX_vp1_genotypes.tsv  (optional; reported only)

    3) Automatic discovery:
         --input-dir DIR
       finds a unique pair of files sharing the same prefix:
         *_summary.tsv + *_contamination.tsv

    The VP1 file is not required by this contamination-network script, but if
    <prefix>_vp1_genotypes.tsv exists it is returned for logging/provenance.
    """
    input_dir = Path(input_dir)

    if (summary_arg is None) ^ (comparisons_arg is None):
        raise ValueError(
            "Provide both --summary and --comparisons, or neither. "
            "If neither is provided, use --run-prefix and/or --input-dir for discovery."
        )

    if summary_arg is not None and comparisons_arg is not None:
        summary = Path(summary_arg)
        comparisons = Path(comparisons_arg)
        if not summary.exists():
            raise FileNotFoundError(f"Summary file not found: {summary}")
        if not comparisons.exists():
            raise FileNotFoundError(f"Contamination file not found: {comparisons}")

        prefix = summary.name[:-len("_summary.tsv")] if summary.name.endswith("_summary.tsv") else summary.stem
        vp1 = summary.parent / f"{prefix}_vp1_genotypes.tsv"
        return summary, comparisons, (vp1 if vp1.exists() else None), prefix

    if not input_dir.exists() or not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found or not a directory: {input_dir}")

    if run_prefix:
        prefix = str(run_prefix)
        summary = input_dir / f"{prefix}_summary.tsv"
        comparisons = input_dir / f"{prefix}_contamination.tsv"
        vp1 = input_dir / f"{prefix}_vp1_genotypes.tsv"

        missing = [str(x) for x in (summary, comparisons) if not x.exists()]
        if missing:
            raise FileNotFoundError(
                "Could not resolve required files for run prefix "
                f"'{prefix}'. Missing: {missing}"
            )
        return summary, comparisons, (vp1 if vp1.exists() else None), prefix

    # Auto-discovery: only accept summary files that have a matching
    # <same-prefix>_contamination.tsv in the same directory.
    candidates = []
    for summary in sorted(input_dir.glob("*_summary.tsv")):
        prefix = summary.name[:-len("_summary.tsv")]
        comparisons = input_dir / f"{prefix}_contamination.tsv"
        if comparisons.exists():
            vp1 = input_dir / f"{prefix}_vp1_genotypes.tsv"
            candidates.append((summary, comparisons, vp1 if vp1.exists() else None, prefix))

    if not candidates:
        raise FileNotFoundError(
            f"No matching '*_summary.tsv' + '*_contamination.tsv' pair found in {input_dir}. "
            "Use --run-prefix PREFIX or provide --summary/--comparisons explicitly."
        )

    if len(candidates) > 1:
        prefixes = ", ".join(x[3] for x in candidates)
        raise ValueError(
            f"Multiple run prefixes were detected in {input_dir}: {prefixes}. "
            "Specify the desired one with --run-prefix PREFIX."
        )

    return candidates[0]

def main():
    p = argparse.ArgumentParser(
        description="Infer rhinovirus amplicon contamination networks and candidate source samples."
    )
    p.add_argument(
        "--summary", default=None,
        help="Explicit summary TSV path. Optional when using --run-prefix/--input-dir discovery."
    )
    p.add_argument(
        "--comparisons", default=None,
        help="Explicit contamination TSV path. Optional when using --run-prefix/--input-dir discovery."
    )
    p.add_argument(
        "--input-dir", default=".",
        help="Directory containing <run_prefix>_summary.tsv and <run_prefix>_contamination.tsv (default: current directory)."
    )
    p.add_argument(
        "--run-prefix", default=None,
        help="User-assigned run/file prefix, e.g. 20260806 or QiaCap_49-104. If omitted, a unique matching prefix is auto-detected."
    )
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument(
        "--metric", choices=["c1","ref"], default="c1",
        help="Abundance metric: c1 (default) or ref."
    )
    p.add_argument("--source-min-pct", type=float, default=20.0,
                   help="Minimum abundance %% for a candidate source (default: 20).")
    p.add_argument("--contaminant-max-pct", type=float, default=10.0,
                   help="Maximum abundance %% for automatic contaminant removal (default: 10).")
    p.add_argument("--source-ratio", type=float, default=5.0,
                   help="Minimum source:recipient abundance ratio (default: 5x).")
    p.add_argument("--high-conf-ratio", type=float, default=10.0,
                   help="Ratio used for HIGH confidence (default: 10x).")
    p.add_argument("--high-conf-identity", type=float, default=99.95,
                   help="Identity %% used for HIGH confidence (default: 99.95).")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_path, comparisons_path, vp1_path, detected_prefix = resolve_input_files(
        summary_arg=args.summary,
        comparisons_arg=args.comparisons,
        input_dir=args.input_dir,
        run_prefix=args.run_prefix,
    )

    print(f"Detected run prefix: {detected_prefix}")
    print(f"Using summary file: {summary_path}")
    print(f"Using contamination file: {comparisons_path}")
    if vp1_path is not None:
        print(f"Detected VP1 genotype file (not used by this script): {vp1_path}")
    else:
        print("VP1 genotype file: not found (optional for this script)")

    summary = pd.read_csv(summary_path, sep="\t")
    comparisons = pd.read_csv(comparisons_path, sep="\t")

    required_summary = {"sample_name","ref_tag","pct_ns"}
    required_comparisons = {
        "sample_a","tag_a","sample_b","tag_b",
        "comparable_sites","differences","pct_identity","flag"
    }
    ms = required_summary - set(summary.columns)
    mc = required_comparisons - set(comparisons.columns)
    if ms:
        raise ValueError(f"{summary_path.name} is missing columns: {sorted(ms)}")
    if mc:
        raise ValueError(f"{comparisons_path.name} is missing columns: {sorted(mc)}")

    nodes, mapping = collapse_same_genome(summary, comparisons, args.metric)
    edges = build_contamination_edges(comparisons, mapping)
    decisions, directed, component_id = infer_network(
        nodes, edges,
        source_min_pct=args.source_min_pct,
        contaminant_max_pct=args.contaminant_max_pct,
        source_ratio=args.source_ratio,
        high_conf_ratio=args.high_conf_ratio,
        high_conf_identity=args.high_conf_identity,
    )
    components = component_summary(nodes, decisions, edges)
    cleanup = make_cleanup_table(nodes, decisions)

    nodes.to_csv(outdir / "collapsed_genome_nodes.tsv", sep="\t", index=False)
    edges.to_csv(outdir / "contamination_edges.tsv", sep="\t", index=False)
    directed.to_csv(outdir / "inferred_source_recipient_edges.tsv", sep="\t", index=False)
    components.to_csv(outdir / "network_components.tsv", sep="\t", index=False)
    cleanup.to_csv(outdir / "cleanup_recommendations.tsv", sep="\t", index=False)

    # Expanded removal manifest at the original sample/reference-hit level.
    removal_rows = []
    for r in cleanup[cleanup["cleanup_action"] == "REMOVE_CONTAMINANT"].itertuples():
        for tag in str(r.collapsed_tags).split(","):
            removal_rows.append({
                "sample_name": r.sample_name,
                "ref_tag": tag,
                "node_id": r.node_id,
                "representative_tag": r.representative_tag,
                "candidate_source_sample": r.candidate_source_sample,
                "candidate_source_node": r.candidate_source_node,
                "confidence": r.confidence,
                "inference_mode": r.inference_mode,
                "reason": r.reason,
            })
    removal_manifest = pd.DataFrame(removal_rows)
    if len(removal_manifest):
        removal_manifest = removal_manifest.drop_duplicates(["sample_name","ref_tag"])
    else:
        removal_manifest = pd.DataFrame(columns=[
            "sample_name","ref_tag","node_id","representative_tag",
            "candidate_source_sample","candidate_source_node","confidence",
            "inference_mode","reason"
        ])
    removal_manifest.to_csv(outdir / "removal_manifest.tsv", sep="\t", index=False)

    # A cleaned copy of run_summary.tsv. This removes only original sample/tag
    # rows represented by REMOVE_CONTAMINANT nodes. Source and unresolved rows
    # are retained.
    remove_keys = set(zip(removal_manifest["sample_name"], removal_manifest["ref_tag"]))
    keep_mask = [
        (row.sample_name, row.ref_tag) not in remove_keys
        for row in summary.itertuples()
    ]
    cleaned_summary = summary.loc[keep_mask].copy()
    cleaned_summary.to_csv(outdir / "cleaned_run_summary.tsv", sep="\t", index=False)

    # Original sample/tag -> collapsed biological observation.
    map_rows = []
    for (sample, tag), node in mapping.items():
        map_rows.append({"sample_name": sample, "ref_tag": tag, "node_id": node})
    pd.DataFrame(map_rows).sort_values(["sample_name","ref_tag"]).to_csv(
        outdir / "same_genome_collapse_map.tsv", sep="\t", index=False
    )

    # Small text report.
    n_remove = int((cleanup["cleanup_action"] == "REMOVE_CONTAMINANT").sum())
    n_source = int((cleanup["cleanup_action"] == "KEEP_SOURCE").sum())
    n_review = int((cleanup["cleanup_action"] == "KEEP_UNRESOLVED").sum())

    report = [
        "Rhinovirus contamination-network analysis",
        "==========================================",
        f"Detected run prefix: {detected_prefix}",
        f"Summary input: {summary_path}",
        f"Contamination input: {comparisons_path}",
        f"VP1 input: {vp1_path if vp1_path is not None else 'not found (optional)'}",
        f"Metric: {args.metric}",
        f"Candidate source minimum abundance: {args.source_min_pct}%",
        f"Automatic contaminant maximum abundance: {args.contaminant_max_pct}%",
        f"Minimum source:recipient abundance ratio: {args.source_ratio}x",
        f"High-confidence ratio: {args.high_conf_ratio}x",
        f"High-confidence identity: {args.high_conf_identity}%",
        "",
        f"Collapsed biological observations: {len(nodes)}",
        f"Contamination edges after collapsing: {len(edges)}",
        f"Inferred source->recipient edges: {len(directed)}",
        f"KEEP_SOURCE: {n_source}",
        f"REMOVE_CONTAMINANT: {n_remove}",
        f"KEEP_UNRESOLVED: {n_review}",
        f"Original summary rows removed: {len(removal_manifest)}",
        f"Original summary rows retained: {len(cleaned_summary)}",
        "",
        "IMPORTANT:",
        "REMOVE_CONTAMINANT is intentionally conservative and should be used as a",
        "decision table. High-abundance or directionally ambiguous observations are",
        "retained rather than deleted automatically.",
    ]
    (outdir / "analysis_report.txt").write_text("\n".join(report) + "\n")

    print("\n".join(report))
    print(f"\nOutputs written to: {outdir}")


if __name__ == "__main__":
    main()
