#!/usr/bin/env python3


import sys
import gzip
from collections import defaultdict
from pathlib import Path
from itertools import groupby
import logging
import pandas as pd



def read_file(fn, open_function=open, mode='r'):
    try:
        with open_function(fn, mode) as f:
            yield from map(str.strip, f)
    except UnicodeDecodeError:
        with gzip.open(fn, mode="rt") as f:
            yield from map(str.strip, f)

def parsed_cds_file(fn):
    lines = read_file(fn)
    yield from map(lambda x: x.split('\t'), lines)

def get_systems(fn):
    systems = defaultdict(lambda: defaultdict(set))
    for subject, target, *_ in map(lambda x: x.split('\t'), read_file(fn)):
        system, component = subject.split('-')
        systems[system][component].add(target)
    return systems

def get_only_positive_contigs(merged_gff, systems):

    def grouby_helper_get_contig(record):
        return record[0], record[1]

    current_assembly = None
    assembly_count = 0

    for (assembly, contig), records in groupby(parsed_cds_file(merged_gff), grouby_helper_get_contig):
        if current_assembly != assembly:
            current_assembly = assembly
            assembly_count += 1
            print(f'Procced: {assembly_count}', file=sys.stderr, end='\r')

        for record in records:            
            *_, current_index, cluster_match = record
            for system, components in systems.items():
                for component, targets in components.items():
                    if cluster_match in targets:
                        yield (assembly, contig, int(current_index), system, component, cluster_match)
    
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

def to_csv_with_comment(df, fn, comment):
    with open(fn, 'w') as f:
        print(comment, file=f)
        df.to_csv(f)

def main(
        search_table,
        cds_info_only_positive_matches,
        system_info,
        total_assemblies,
        window_size,
        simple_freq_fn,
        colloc_freq_fn
    ):

    args = locals()
    logger = setup_logging()

    logger.info(f'Arguments: {args}')    

    output_parent_dir = simple_freq_fn.parents[0]
    logger.info(f'Output directory: {output_parent_dir}')

    logger.info(f'Reading search table: {search_table}')
    systems = get_systems(search_table)
    system_components = dict(map(lambda x: (x[0], set(x[1].keys())), systems.items()))

    logger.info(f'Reading merged gff: {cds_info_only_positive_matches}')
    matches = pd.DataFrame(
        get_only_positive_contigs(cds_info_only_positive_matches, systems),
        columns=['assembly', 'contig', 'on_contig_index', 'system', 'component', 'cluster']
    )

    logger.info(f'Total assemblies: {total_assemblies}')

    logger.info(f'Calculating simple frequency ...')
    simple_freq = get_simple_frequency(matches, total_assemblies)
    logger.info(f'Simple frequency: {simple_freq}')

    logger.info(f'Calculating colloc frequency ...')
    colloc_freq = get_colloc_freq(matches, window_size, system_components, total_assemblies)
    logger.info(f'Colloc frequency: {colloc_freq}')

    logger.info(f'Writing output to: {output_parent_dir}')
    to_csv_with_comment(simple_freq, simple_freq_fn, f'# window_size: {window_size}\n# total_assemblies: {total_assemblies}')
    to_csv_with_comment(colloc_freq, colloc_freq_fn, f'# window_size: {window_size}\n# total_assemblies: {total_assemblies}')

    logger.info('Done !')


if __name__ == '__main__':
    args = sys.argv[1:]
    search_table, cds_info_only_positive_matches, system_info, *_, simple_freq_fn, colloc_freq_fn = map(Path, args)
    total_assemblies, window_size = map(int, args[3:5])
    main(
        search_table=search_table,
        cds_info_only_positive_matches=cds_info_only_positive_matches,
        system_info=system_info,
        total_assemblies=total_assemblies,
        window_size=window_size,
        simple_freq_fn=simple_freq_fn,
        colloc_freq_fn=colloc_freq_fn
    )