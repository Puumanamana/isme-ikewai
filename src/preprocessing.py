#!/usr/bin/env python3

import argparse
from pathlib import Path
import logging

import pandas as pd
import numpy as np


logging.basicConfig(level=logging.INFO,
                    format='{asctime} {levelname}: {message}',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    style='{')

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('--abundance', type=str, help='Path to abundance data (mothur .shared format)')
    parser.add_argument('--taxonomy', type=str, help='Path to taxonomy data (mothur .taxonomy format)')
    parser.add_argument('--metadata', type=str, help='Path to metadata (.csv format)')
    parser.add_argument('--output', type=str, default='output', help='Folder path to store outputs')
    parser.add_argument('--discard-singletons', action='store_true', help='Discard OTUs that appear in only 1 sample')            
    args = parser.parse_args()

    return args

def load_abundance(path, discard_singletons=True):
    """
    Loads abundance table in path, sum replicates

    Args:
        path (str): path to .shared abundance table from mothur
        discard_singletons (bool): Discard OTUs that are in less than 2 samples
    Returns:
        pandas.DataFrame
    """
    header = next(open(path)).split('\t')

    data = pd.read_csv(
        path, sep='\t', index_col=0,
        dtype=dict((x, str) if x=='Group' else (x, int) for x in header),
        low_memory=False,
        usecols=lambda x: x not in {'label', 'numOtus'})

    # Group replicates
    data = data.groupby(data.index.str.replace('_\d$', '')).sum()

    # Filter singletons
    if discard_singletons:
        is_singleton = (data > 0).sum(axis=0) < 2
        data = data.loc[:, ~is_singleton]

    return data

def load_taxonomy(path):
    """
    Loads taxonomy table in path

    Args:
        path (str): path to .taxonomy table from mothur
    Returns:
        pandas.DataFrame
    """
    
    data = pd.read_csv(path, index_col=0, sep='\t').Taxonomy
    data = (
        data.str.strip(';')
        .str.split(';', expand=True)
        .set_axis(['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'], axis=1)
    )
    
    return data


def load_metadata(path, categorical_cols):
    """
    Loads metadata table in path, convert columns in `categorical_cols`. 
    Average replicate values for quantitative columns, take the first value for qualitative columns
    
    Args:
        path (str): path to .csv metadata file. Row names need to match abundance rows
        categorical_cols (list): list of numeric columns to be treated as factors
    Returns:
        pandas.DataFrame
    """
    
    data = pd.read_csv(path, index_col=0)
    data['eruption'] = [
        'pre-eruption' if year == 2017 or (year==2018 and season=='March')
        else 'post-eruption' for (year, season) in data[['year', 'season']].values
    ]

    categorical_cols = set(categorical_cols).union(
        data.select_dtypes(exclude='number').columns
    )

    # Group replicates
    data = (
        data.groupby(data.index.str.replace('_\d$', ''))
        .agg(dict((col, 'first') if col in categorical_cols else (col, 'mean')
                  for col in data.columns))
    )

    return data

def main():

    args = parse_args()

    logging.info('Loading abundance table')
    abundance = load_abundance(args.abundance, discard_singletons=args.discard_singletons)
    logging.info('Loading taxonomy table')    
    taxonomy = load_taxonomy(args.taxonomy)
    logging.info('Loading metadata table')    
    metadata = load_metadata(args.metadata, ['site_id', 'year'])

    logging.info('Coalesce data tables')
    common_samples = abundance.index.intersection(metadata.index)
    common_otus = abundance.columns.intersection(taxonomy.index)

    # Remove any OTUs that appear in no samples
    null_otus = abundance.loc[common_samples, common_otus].sum(axis=0) == 0
    common_otus = common_otus[~null_otus]

    logging.info(f'Write table to disk')
    # Write final tables to disk
    Path(args.output).mkdir(exist_ok=True, parents=True)
    abundance.loc[common_samples, common_otus].to_csv(f'{args.output}/abundance.csv')
    taxonomy.loc[common_otus].to_csv(f'{args.output}/taxonomy.csv')
    metadata.loc[common_samples].to_csv(f'{args.output}/metadata.csv')

if __name__ == '__main__':
    main()
