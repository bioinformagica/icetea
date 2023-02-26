#!/usr/bin/env python3

import sys
from pathlib import Path
import pickle
from collections import defaultdict
import networkx as nx
import logging


def read_mmseqs_cluster_tsv(filename):
    """
    Reads a mmseqs cluster tsv file and returns a dictionary
    with the member_id as key and the cluster_id as value
    """
    with open(filename, 'r', encoding='utf-8') as f:
        return dict(map(lambda x: x.strip().split('\t')[::-1], f))

def dict_to_pickle(d, filename):
    """
    Saves a dictionary to a pickle file
    """
    with open(filename, 'wb') as f:
        pickle.dump(d, f)

def pickle_to_dict(filename):
    """
    Reads a dictionary from a pickle file
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)


def get_protein_id_from_gff_info(info):
    """ Returns the protein_id from a gff info field """
    return dict(map(lambda x: x.split('='), info.split(';'))).get('protein_id')


def parse_merged_gff(fn: str):
    """ Parses a merged gff file and returns a generator of tuples """

    with open(fn) as f:
        lines = map(str.strip, f)

        while (line := next(lines, None)):
            if line.startswith('#!genome-build-accession NCBI_Assembly:'):
                assembly = line.split(':')[-1]
                contigs = defaultdict(list)
                continue

            if line.startswith('###'):
                contigs = dict(map(lambda x: (x[0], tuple(x[1])), contigs.items()))
                yield (assembly, contigs)
                continue

            if line.startswith('#'):
                continue

            contig, *_, info = line.split('\t')

            if (protein_id := get_protein_id_from_gff_info(info)): #### <-- Walrus operator
                contigs[contig].append(protein_id)


def main(merged_gff, cluster_table, output_pickle_network):

    logger = logging.getLogger(Path(__file__).stem)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s | %(name)s | %(message)s')

    # print log to stderr
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)



    # Create a dictionary of the cluster table
    pickled_cluster_table = Path(cluster_table).with_suffix('.pickle')

    if pickled_cluster_table.exists():

        # print(f'Loading cluster table from {pickled_cluster_table}')
        logger.info(f'Loading cluster table from {pickled_cluster_table}')
        clusters_members = pickle_to_dict(pickled_cluster_table)
        logger.info(f'Done loading cluster table from {pickled_cluster_table}')

    else:
        logger.info(f'Creating cluster table from {cluster_table}')
        clusters_members = read_mmseqs_cluster_tsv(cluster_table)
        logger.info(f'Done creating cluster table from {cluster_table}')

        logger.info(f'Saving cluster table to {pickled_cluster_table}')
        dict_to_pickle(clusters_members, pickled_cluster_table)
        logger.info(f'Done saving cluster table to {pickled_cluster_table}')


    logger.info(f'Creating network from {merged_gff}')
    number_of_processed_asssemblies = 0

    edges = defaultdict(set)
    nodes = set()

    for i in parse_merged_gff(merged_gff):
        assembly = i[0]
        for contig, proteins in i[1].items():

            clusters = tuple(map(clusters_members.get, proteins))
            
            nodes.update((contig_nodes := set(filter(None, clusters))))

            if len(contig_nodes) <= 1:
                continue
            
            #  zip: get edges: (1, 2, 3, 4) -> (1, 2), (2, 3), (3, 4)
            #  filter: some times is None so we filter out those
            #  filter: because of fragmentation and duplication,
            #          some times the same cluster is repeated
            #          so we filter out those too. The information
            #          not lost because: (1, 2, 2, 3) -> (1, 2), (2, 3)
            #  map: we sort the tuple so direction doesn't matter
            #  map: we convert the tuple to a tuple so it can be a set element
            #  set: we convert the list to a set so we can filter out duplicates
            contig_edges = set(map(tuple, map(sorted, filter(lambda x: x[0] != x[1], filter(all, zip(clusters, clusters[1:]))))))

            if len(contig_edges) < 1:
                continue

            for edge in contig_edges:
                edges[edge].add((assembly, contig))
                
        number_of_processed_asssemblies += 1
        # print counter every 100 assemblies
        if number_of_processed_asssemblies % 100 == 0:
            logger.info(f'processed {number_of_processed_asssemblies} assemblies')


    logger.info(f'Done creating network from {merged_gff}')
    logger.info(f'processed {number_of_processed_asssemblies} assemblies')

    logger.info('Creating networkx graph')
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(( *k, {'contigs' : tuple(v)}) for k, v in edges.items())
    logger.info('Done creating networkx graph')

    logger.info(f'Saving network to {output_pickle_network}')
    with open(output_pickle_network, 'wb') as f:
        pickle.dump(G, f)
    logger.info(f'Done saving network to {output_pickle_network}')



if __name__ == '__main__':
    main(*sys.argv[1:])
