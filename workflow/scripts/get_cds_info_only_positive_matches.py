#!/usr/bin/env python3

import sys
import gzip
from pathlib import Path
from itertools import groupby
import logging

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


def read_file(fn, open_function=open, mode='r'):
    try:
        with open_function(fn, mode) as f:
            yield from map(str.strip, f)
    except UnicodeDecodeError:
        with gzip.open(fn, mode="rt") as f:
            yield from map(str.strip, f)



def main(clusters_assigment_file, merged_gff, output):
    args = locals()
    logger = setup_logging()
    logger.info(f"Starting {__file__} with args: {args} ...")

    logger.info(f"Reading {clusters_assigment_file} ...")
    clusters_assigment = dict(map(lambda x: x.split()[::-1], read_file(clusters_assigment_file)))

    logger.info(f"Reading {merged_gff} ...")
    gff_records = map(lambda x: tuple(x.split('\t')), read_file(merged_gff))
    contigs = groupby(gff_records, lambda x: tuple(x[:2]))
    cds_enumerated = map(lambda x: (x[0], enumerate(x[1], 1)), contigs)
    
    total_assemblies = 0
    current_assembly = None

    with open(output, 'w') as f:
        for (assembly, _), cdss in cds_enumerated:
            if assembly != current_assembly:
                total_assemblies += 1
                current_assembly = assembly
                print(f'Prosessed {total_assemblies} assemblies', end='\r', file=sys.stderr)
            for cds in cdss:
                if (cluster := clusters_assigment.get(cds[-1][-1])):
                    f.write('\t'.join((*cds[-1], str(cds[0]), cluster, '\n')))

    logger.info(f"Done ! Total assemblies: {total_assemblies}")


if __name__ == '__main__':
    try:
        main(*sys.argv[1:])
    except (KeyboardInterrupt, BrokenPipeError) as e:
        sys.exit(0)