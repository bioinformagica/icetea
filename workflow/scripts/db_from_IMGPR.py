#!/usr/bin/env python3

import sys
from pathlib import Path
import gzip
import sqlite3
from itertools import tee
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
    print(f'Usage: {scrip_name} <cds_info_with_cluster_assignment> <output_db>', file=sys.stderr)

def read_file(fn, open_function=open, mode='r'):
    try:
        with open_function(fn, mode) as f:
            yield from map(str.strip, f)
    except UnicodeDecodeError:
        with gzip.open(fn, mode="rt") as f:
            yield from map(str.strip, f)

def format_record(record):
    genome, _, locus_tag, cds_index, cluster_id = record.split('\t')
    return (locus_tag, cluster_id, genome, int(cds_index))

def main(cds_info_with_cluster_assignment, output_db):
    args = locals()
    logger = setup_logging()
    logger.info(f"Starting {__file__} with args: {args} ...")

    logger.info(f"Reading {cds_info_with_cluster_assignment} ...")
    records =  map(format_record, read_file(cds_info_with_cluster_assignment))

    logger.info(f"Writing {output_db} ...")
    with sqlite3.connect(output_db) as conn:
        c = conn.cursor()
        c.execute("DROP TABLE IF EXISTS cds_info_with_cluster_assignment")
        # Use the CREATE TABLE statement to create a new table with primary key locus_tag
        c.execute(
            """
            CREATE TABLE cds_info_with_cluster_assignment (
                locus_tag TEXT PRIMARY KEY,
                cluster_id TEXT,
                genome TEXT,
                cds_index INTEGER
            )
            """
        )

        # Create an index on the cluster_id column
        c.execute(
            """
            CREATE INDEX cluster_id_index
            ON cds_info_with_cluster_assignment (cluster_id)
            """
        )

        # apply the records to the database
        c.executemany(
            """
            INSERT INTO cds_info_with_cluster_assignment
            VALUES (?, ?, ?, ?)
            """,
            records
        )

        conn.commit()

    logger.info(f"Done!")

if __name__ == '__main__':
    args = sys.argv[1:]
    print(args)
    if len(args) != 2:
        usage(sys.argv[0])
        sys.exit(1)

    cds_info_with_cluster_assignment, output_db = args

    main(
        cds_info_with_cluster_assignment,
        output_db   
    )