#!/usr/bin/env python3

import sys
from pathlib import Path
import pickle
from collections import namedtuple, defaultdict
from typing import Dict, List, Tuple
from itertools import chain
import logging
from functools import partial

import networkx as nx
import pandas as pd


def get_graph_description(G: nx.Graph) -> Dict[str, int]:
    """Return a dictionary with the number of nodes and edges in the graph."""
    return {
        "number_of_nodes": G.number_of_nodes(),
        "number_of_edges": G.number_of_edges(),
    }


def read_mmseqs_search_output(
    fn: str, usecols: str = "query,target,pident,qcov,tcov"
) -> pd.DataFrame:
    header = "query,target,pident,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov".split(
        ","
    )
    return pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=header,
        usecols=usecols.split(","),
    )


def read_network_from_pickle(fn: str) -> nx.Graph:
    with open(fn, "rb") as f:
        return pickle.load(f)


def get_subgraph(G: nx.Graph, start_nodes: list, n_edges: int = 5) -> nx.Graph:
    """Returns a subgraph of G with n_edges of distance from start_nodes"""
    connected_nodes = set()
    for node in start_nodes:
        connected_nodes.update(nx.ego_graph(G, node, radius=n_edges).nodes())

    return G.subgraph(connected_nodes)


def query_network(G: nx.Graph, query_node: str, radius: int) -> set:
    try:
        return set(nx.ego_graph(G, query_node, radius=radius).nodes())
    except nx.exception.NodeNotFound:
        return set()


def get_node_info(G: nx.Graph, node: str) -> List[Tuple[int, str, str]]:
    return G.nodes(data=True)[node]["info"]


def get_edge_info(G: nx.Graph, edge: str) -> List[Tuple[int, str, str]]:
    return G.edges(data=True)[edge]["info"]


def extract_edge_info(G: nx.Graph) -> Tuple[str, str, List[Tuple[str, str, int]]]:
    data = defaultdict(list)

    for edge in G.edges(data=True):
        source, target, infos = edge
        for info in infos["info"]:
            oncontig_pos, assembly, contig = info
            data[(assembly, contig)].append((source, target, oncontig_pos))

    for (assembly, contig), infos in data.items():
        yield assembly, contig, infos


def format_edge_info(
    edge_info: Tuple[str, str, List[Tuple[str, str, int]]]
) -> pd.DataFrame:
    return (
        pd.DataFrame(
            data=edge_info[2],
            columns=["gsource", "gtarget", "on_contig_position"],
        )
        .sort_values("on_contig_position")
        .assign(
            windows=lambda df_: df_.on_contig_position.diff().ne(1).cumsum(),
            assembly=edge_info[0],
            contig=edge_info[1],
        )
        .groupby(["assembly", "contig", "windows"])
        .apply(lambda df_: set(chain.from_iterable(df_[["gsource", "gtarget"]].values)))
        .rename("nodes")
        .reset_index()
    )


def split_col_to_multiindex(
    df: pd.DataFrame, sep: str = "-", names: None | list = None
) -> pd.DataFrame:
    df.columns = pd.MultiIndex.from_tuples(
        df.columns.to_series().str.split(sep).apply(tuple), names=names
    )
    return df


def get_simple_freq(df: pd.DataFrame, total_assemblies: int) -> pd.DataFrame:
    return (
        df.droplevel(["contig", "windows"], axis=0)
        .groupby(["assembly"])
        .sum()
        .applymap(bool)
        .applymap(int)
        .sum()
        .apply(lambda x: x / total_assemblies)
        .rename("frequency")
        .to_frame()
    )


def get_colloc_freq(df: pd.DataFrame, total_assemblies: int) -> pd.DataFrame:
    return (
        df.groupby(["system"], axis=1)
        .all()
        .droplevel(["contig", "windows"], axis=0)
        .groupby(["assembly"])
        .any()
        .sum()
        .apply(lambda x: x / total_assemblies)
        .rename("frequency")
        .to_frame()
    )


EdgeInfo = namedtuple("EdgeInfo", ["on_contig_position", "assembly", "contig"])


def log_df(
    df: pd.DataFrame,
    handler: logging.Logger,
    *message: str,
    level: str = "info",
) -> pd.DataFrame:
    getattr(handler, level)(*message)
    return df


def main(
    search_table, ge_pident, ge_cov, G_pickle, output_dir, query_radius: int = 5
) -> None:
    ge_pident = float(ge_pident)
    ge_cov = float(ge_cov)
    output_dir = Path(output_dir)

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s | %(name)s | %(message)s",
        handlers=[logging.StreamHandler()],
    )

    logger = logging.getLogger(Path(__file__).stem)

    logger.info(f"Reading network from pickle file: {G_pickle} ...")
    G = read_network_from_pickle(G_pickle)
    total_assemblies = len(
        set(
            map(
                lambda x: x.assembly,
                chain.from_iterable(
                    map(lambda x: x[-1].get("info"), G.edges(data=True))
                ),
            )
        )
    )

    logger.info(f"Reading mmseqs search table: {search_table} ...")

    query_network_preloaded = partial(query_network, G, radius=query_radius)

    system_matches = (
        read_mmseqs_search_output(search_table)
        .loc[
            lambda df_: (
                df_.pident.ge(ge_pident) & df_.qcov.ge(ge_cov) & df_.tcov.ge(ge_cov)
            ),
            ["query", "target"],
        ]
        .pipe(
            lambda df_: pd.DataFrame(
                data=df_.target.values,
                index=pd.MultiIndex.from_tuples(
                    df_["query"].apply(lambda x: x.split("-")),
                    names=["system", "component"],
                ),
                columns=["target"],
            )
        )
        .pipe(log_df, logger, "Quering network ...", level="info")
        .assign(
            close_nodes=lambda df_: df_.target.apply(query_network_preloaded),
        )
        .groupby(["system", "component"])
        .agg(
            {
                "target": set,
                "close_nodes": lambda s_: set(chain.from_iterable(s_)),
            }
        )
        .loc[lambda df_: df_.close_nodes.apply(bool)]
    )

    all_simple_freq = []
    all_colloc_freq = []

    for i, system_df in system_matches.groupby("system"):
        logger.info(f"Processing system: {i} ...")

        system_G = G.subgraph(set(chain.from_iterable(system_df.close_nodes)))

        positive_matches = set(chain.from_iterable(system_df.target))

        edge_infos = extract_edge_info(system_G)

        remove_contigs_without_matches = filter(
            lambda tup: (
                set(chain.from_iterable((i[:2] for i in tup[2]))) & positive_matches
            ),
            edge_infos,
        )

        as_dfs = map(format_edge_info, remove_contigs_without_matches)

        all_positive_windows = (
            pd.concat(
                as_dfs,
                ignore_index=True,
            )
            .loc[lambda df_: df_.nodes.apply(lambda x: bool(x & positive_matches))]
            .reset_index(drop=True)
        )

        logger.info(f"Getting frequency for system: {i} ...")
        prep_freq = (
            system_df.assign(
                merged_index=lambda df_: df_.index.to_series().apply(
                    lambda x: "-".join(x)
                ),
            )
            .set_index("merged_index")
            .pipe(
                lambda df_: pd.DataFrame(
                    data=zip(
                        *map(
                            lambda x: all_positive_windows.nodes.apply(
                                lambda y: bool(x & y)
                            ),
                            df_.target,
                        )
                    ),
                    columns=df_.index.values,
                )
            )
            .pipe(lambda df_: pd.concat([all_positive_windows, df_], axis=1))
            .drop("nodes", axis=1)
            .set_index(["assembly", "contig", "windows"])
            .pipe(split_col_to_multiindex, names=["system", "component"])
        )

        all_simple_freq.append(get_simple_freq(prep_freq, total_assemblies))
        all_colloc_freq.append(get_colloc_freq(prep_freq, total_assemblies))

    logger.info("Saving results ...")
    pd.concat(all_simple_freq).to_csv(output_dir / "simple_freq.tsv", sep="\t")
    pd.concat(all_colloc_freq).to_csv(output_dir / "colloc_freq.tsv", sep="\t")


if __name__ == "__main__":
    main(*sys.argv[1:])
