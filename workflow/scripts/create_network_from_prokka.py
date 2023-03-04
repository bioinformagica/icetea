#!/usr/bin/env python3

# Native
import sys
from pathlib import Path
import pickle
from collections import defaultdict
from typing import Generator
import logging

# External
import networkx as nx


def pickle_to_dict(fn: str) -> dict:
    """ Loads a pickle file and returns a dictionary """
    with open(fn, 'rb') as f:
        return pickle.load(f)

def get_protein_id_from_gff_info(info: str) -> str:
    """ Returns the protein_id from a gff info field """
    return dict(map(lambda x: x.split('='), info.split(';'))).get('locus_tag')

def parse_merged_gff(fn: str) -> Generator[tuple[str], None, None]:
    """ Parses a merged gff file and returns a generator of tuples """

    with open(fn) as f:
        lines = map(str.strip, f)
        remove_headers = filter(lambda x: not x.startswith('#'), lines)

        while (line := next(remove_headers, None)):
            contig, *_, info = line.split('\t')
            protein_id = get_protein_id_from_gff_info(info)
            assembly, contig = contig.split('_')
            yield (assembly, contig, protein_id)

def mmseqs_to_dict(fn: str) -> dict:
    """ Parses a mmseqs cluster file and returns a dictionary """
    with open(fn, 'r', encoding='utf-8') as f:
        lines = map(str.strip, f)
        return dict(map(lambda x: x.strip().split('\t'), lines))

def get_edges(locus_tags: tuple) -> set:
    """ Returns a set of edges from a tuple of locus tags """
    edges = zip(locus_tags, locus_tags[1:])
    remove_fragments = filter(lambda x: x[0] != x[1], edges)
    remove_none = filter(all, remove_fragments)
    add_on_contig_position = map(lambda x: (*x[-1], x[0]), enumerate(remove_none, start=1))
    return tuple(add_on_contig_position)

def dict_from_pickle(fn: str) -> dict:
    """ Loads a pickle file and returns a dictionary """
    with open(fn, 'rb') as f:
        return pickle.load(f)

def init_logger(default_level: str) -> logging.Logger:
    """ Initializes a logger object """

    corredted_level = default_level.upper()
    allowed_levels = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}
    assert corredted_level in allowed_levels, f'Logging level {corredted_level} not allowed. Allowed levels are {allowed_levels}'

    logging.basicConfig(
        level=getattr(logging, default_level.upper()),
        format="%(asctime)s | %(name)s | %(message)s",
        handlers=[logging.StreamHandler()],
    )

    logger = logging.getLogger(Path(__file__).stem)

    return logger



def main(clusters_mmseqs: str, merged_gff: str, output: str, logging_level: str = 'info', report_every: int = 100):

    logger = init_logger(logging_level)

    logger.info('Parsing mmseqs cluster file')
    cluster_assignments = mmseqs_to_dict(clusters_mmseqs)
    
    logger.info('Parsing merged gff file')
    all_cds = parse_merged_gff(merged_gff)

    add_cluster_assignments = map(lambda x: (*x[:-1], cluster_assignments.get(x[-1])), all_cds)

    get_contigs = defaultdict(list)

    counter_genomes = 0
    all_genomes = set()

    for i in add_cluster_assignments:
        get_contigs[i[:-1]].append(i[-1])
        if not i[0] in all_genomes:
            all_genomes.add(i[0])
            counter_genomes += 1
            if counter_genomes % report_every == 0:
                logger.info(f'Parsed {counter_genomes} genomes')

    logger.info(f'Parsed {counter_genomes} genomes in total')

    get_contigs = {k: get_edges(v) for k, v in get_contigs.items()}

    logger.info(f'Creating network and writing to {output} ...')
    edges =  defaultdict(list)

    for k, v in get_contigs.items():
        for edge in v:
            edges[edge[:-1]].append((*k, edge[-1]))

    G = nx.Graph()

    G.add_edges_from(map(lambda x: (*x[0], {'info' : tuple(x[1])}), edges.items()))

    with open(output, 'wb') as f:
        pickle.dump(G, f)

    logger.info('Done !')
if __name__ == '__main__':
    main(*sys.argv[1:])