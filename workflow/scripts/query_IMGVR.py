#!/usr/bin/env python3

import sys
from pathlib import Path
import gzip
import sqlite3
import logging

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

def usage(scrip_name):
    print(f'Usage: {scrip_name} <filtered_search_mmseqs2.tsv> <cds_db>', file=sys.stderr)

def read_file(fn, open_function=open, mode='r'):
    try:
        with open_function(fn, mode) as f:
            yield from map(str.strip, f)
    except UnicodeDecodeError:
        with gzip.open(fn, mode="rt") as f:
            yield from map(str.strip, f)

def main(filtered_search_mmseqs2, cds_db):
    positive_matches = ((i,) for i in set(map(lambda x: x.split('\t')[1], read_file(filtered_search_mmseqs2))))

    with sqlite3.connect(cds_db) as conn:
        c = conn.cursor()

        c.execute(
            """
            CREATE TEMPORARY TABLE clusters (
                cluster_id TEXT PRIMARY KEY
            )
            """
        )

        c.executemany(
            """
            INSERT INTO clusters (cluster_id) VALUES (?)
            """,
            positive_matches
        )

        c.execute(
            """
            SELECT * FROM cds_info_with_cluster_assignment
            WHERE cluster_id IN (
                SELECT cluster_id FROM clusters
            )
            """
        )

    for result in c:
        print(*result, sep='\t', file=sys.stdout)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 2:
        usage(sys.argv[0])
        sys.exit(1)

    filtered_search_mmseqs2, cds_db = args

    main(filtered_search_mmseqs2, cds_db)