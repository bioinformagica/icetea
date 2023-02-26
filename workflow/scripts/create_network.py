#!/usr/bin/env python3

import sys
from pathlib import Path
import pickle
from collections import defaultdict
import networkx as nx
import logging
from collections import namedtuple


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


EdgeInfo = namedtuple('EdgeInfo', ['on_contig_position', 'assembly', 'contig'])

def main(merged_gff, cluster_table, output_pickle_network, log_each=100):

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

            removing_duplicates = tuple({cluster : None for cluster in clusters if cluster })
            
            nodes.update(removing_duplicates)

            if len(removing_duplicates) <= 1:
                continue

            get_edges = zip(removing_duplicates, removing_duplicates[1:])

            sort_each_tuple = map(tuple, map(sorted, get_edges))

            add_edge_position = map(lambda x: (*x[1], x[0]), enumerate(sort_each_tuple, start=1))

            as_edge_namedtuple = map(lambda x: (x[:2], EdgeInfo(x[2], assembly, contig)), add_edge_position)

            for edge_pos, edge_info in as_edge_namedtuple:
                edges[edge_pos].add(edge_info)
                
        number_of_processed_asssemblies += 1
        if number_of_processed_asssemblies % log_each == 0:
            logger.info(f'processed {number_of_processed_asssemblies} assemblies')


    logger.info(f'Done creating network from {merged_gff}')
    logger.info(f'processed {number_of_processed_asssemblies} assemblies')

    logger.info('Creating networkx graph')
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(( *k, {'info' : tuple(v)}) for k, v in edges.items())
    logger.info('Done creating networkx graph')

    logger.info(f'Saving network to {output_pickle_network}')
    with open(output_pickle_network, 'wb') as f:
        pickle.dump(G, f)
    logger.info(f'Done saving network to {output_pickle_network}')



if __name__ == '__main__':
    main(*sys.argv[1:])
