#!/usr/bin/env python3

import sys
from collections import defaultdict
import pandas as pd
from pathlib import Path
import json
import logging


def read_file(fn):
    with open(fn, 'r') as f:
        yield from map(str.strip, f)

def read_search_table(fn):
    yield from map(lambda x: x.split('\t')[:2], read_file(fn))
    
    
def get_systems(fn):
    systems = defaultdict(lambda: defaultdict(set))
    for subject, target in read_search_table(fn):
        system, component = subject.split('-')
        systems[system][component].add(target)
    return systems

def get_protein_id_from_gff_info(info):
    """ Returns the protein_id from a gff info field """
    return dict(map(lambda x: x.split('='), info.split(';'))).get('protein_id', 'NA')

def parse_merged_gff(fn: str):
    """ Parses a merged gff file and returns a generator of tuples """

    lines = read_file(fn)

    while (line := next(lines, None)):
        if line.startswith('#!genome-build-accession NCBI_Assembly:'):
            assembly = line.split(':')[-1]
            continue

        if line.startswith('#'):
            continue

        contig, *_, info = line.split('\t')
        yield assembly, contig, get_protein_id_from_gff_info(info)

def mmseqs_to_dict(fn: str) -> dict:
    return dict(map(lambda x: x.strip().split('\t')[::-1], read_file(fn)))

def search_matches(merged_gff, systems, cluster_map):
    """ Search matches in the gff file and returns a generator of tuples """
    current_assembly_contig = None

    current_assembly = None
    assembly_count = 0

    for assembly, contig, protein_id in parse_merged_gff(merged_gff):

        if current_assembly != assembly:
            current_assembly = assembly
            assembly_count += 1
            if assembly_count % 100 == 0:
                print(f'Processed: {assembly_count}', end='\r', file=sys.stderr)


        if current_assembly_contig != (assembly, contig):
            current_assembly_contig = (assembly, contig)
            current_index = 0

        current_index += 1
        as_cluster = cluster_map.get(protein_id, 'NA')

        for system, components in systems.items():
            for component, targets in components.items():
                if as_cluster in targets:
                    yield (assembly, contig, current_index, system, component, as_cluster)
    

def get_simple_frequency(df, total_assemblies):
    return (
        df
        .loc[:, ['assembly', 'system', 'component']]
        .drop_duplicates()
        .groupby(['system', 'component'])
        .size()
        .apply(lambda x: x / total_assemblies)
        .rename('frequency')
        .to_frame()
    )



def get_windows(diff, window_size):
    """ Returns a dataframe with the windows for each system """
    result = []

    for i in diff:
        if not result:
            slide = 1
            result.append(slide)
            continue
            
        if i <= window_size:
            result.append(slide)
            continue
            
        slide += 1
        result.append(slide)

    return result

def get_colloc_freq(df, window_size, system_components, total_assemblies):
    return (
        df
        .groupby(['assembly', 'contig', 'system'],  group_keys=False)
        .apply(
            lambda df_: (
                df_.assign(
                    window=lambda df: get_windows(df.on_contig_index.diff().fillna(False), window_size)
                )
            )
        )
        .loc[:, ['assembly', 'contig', 'system', 'component', 'window']]
        .groupby(['assembly', 'contig', 'system', 'window'])
        .agg(set)
        .reset_index()
        .assign(
            is_complete = lambda df: df[['system','component']].apply(lambda x: x[1] == system_components[x[0]], axis=1)
        )
        .loc[:, ['assembly', 'system', 'is_complete']]
        .groupby(['assembly', 'system'])
        .any()
        .reset_index()
        .loc[lambda df: df.is_complete, ['assembly', 'system']]
        .groupby('system')
        .size()
        .rename('frequency')
        .apply(lambda x: x / total_assemblies)
        .to_frame()
    )

def setup_logging(name = Path(__file__).stem.replace('.py', ''), level='info'):
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s | %(name)s | %(levelname)s | %(message)s',
        datefmt='%Y-%m-%d %H:%M',
        handlers=[
            logging.StreamHandler(sys.stderr)
        ]
    )

    return logging.getLogger(name)



def main(merged_gff, search_table, cluster_map, total_assemblies, output_prefix, window_size = 5):
    args = locals()
    logger = setup_logging()


    # logging arguments
    logger.info(f'Arguments: {json.dumps(args, indent=4)}')

    logger.info('Reading cluster map ...')
    cluster_map = mmseqs_to_dict(cluster_map)

    total_assemblies = int(total_assemblies)
    output_prefix = Path(output_prefix)

    logger.info('Reading systems ...')
    systems = get_systems(search_table)
    system_components = dict(map(lambda x: (x[0], set(x[1].keys())), systems.items()))

    logger.info('Reading matches ...')
    matches = pd.DataFrame(
        search_matches(merged_gff, systems, cluster_map),
        columns=['assembly', 'contig', 'on_contig_index', 'system', 'component', 'cluster']
    )

    logger.info('Processing simple frequency ...')
    get_simple_frequency(matches, total_assemblies).to_csv(output_prefix / 'simple_frequency.csv')

    logger.info('Processing colloc frequency ...')
    get_colloc_freq(matches, window_size, system_components, total_assemblies).to_csv(output_prefix / 'colloc_frequency.csv')

    logger.info('Done !')

if __name__ == '__main__':
    main(*sys.argv[1:])