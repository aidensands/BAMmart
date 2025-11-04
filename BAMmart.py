import pandas as pd
from pybiomart import Dataset
import pysam
import os
import sys
import argparse


def find_bam_files(root_dir: str) -> list:
    """
    Recursively walks through root_dir and finds all files ending with .bam.
    """
    bam_files = []
    print(f"Checking for .bam files...")
    for dirpath, _, filenames in os.walk(root_dir):
        for f in filenames:
            if f.endswith('.bam'):
                full_path = os.path.join(dirpath, f)
                bam_files.append(full_path)
    print(f"Found {len(bam_files)} .bam files.")
    return bam_files


def term_help(search_term: str) -> None:
    """Takes a search term and returns a list of args containing the given search term"""
    database = Dataset(name='hsapiens_gene_ensembl',
                       host='http://www.ensembl.org')
    attributes = list()
    filters = list()
    for term in database.attributes:
        if search_term in term:
            attributes.append(term)

    for term in database.filters:
        if search_term in term:
            filters.append(term)

    if len(attributes) != 0:
        print("Attributes: ")
        print(attributes)
    else:
        print("No Attributes")

    if len(filters) != 0:
        print("Filters: ")
        print(filters)
    else:
        print("No Filters")


def parse_enst(bam) -> list:
    """Extract unique ENST IDs from BAM file's RNAME field using pysam."""
    ids = set()
    print(f"1. Extracting IDs from {bam}...")
    try:
        with pysam.AlignmentFile(bam, 'rb') as samfile:
            for read in samfile.fetch(until_eof=True):
                if not read.is_unmapped and read.reference_name:
                    ids.add(read.reference_name)

        print(f"   Found {len(ids)} unique ENST IDs.")
        return list(ids)
    except FileNotFoundError:
        print(f"Error: BAM file not found at {bam}")
        return []
    except Exception as e:
        print(f"An error occurred during pysam processing for {bam}: {e}")
        return []


def parse_ensg(bam):
    ids = set()
    try:
        with pysam.AlignmentFile(filename=bam, mode='rb') as samfile:
            for read in samfile.fetch(until_eof=True):
                if read.is_unmapped:
                    continue

                gene_id = None
                for tag in ['GX', 'GE', 'GN']:
                    try:
                        val = read.get_tag(tag)
                        if val.startswith('ENSG'):
                            gene_id = val
                            break

                    except KeyError:
                        continue

                if gene_id:
                    ids.add(gene_id)

    except FileNotFoundError:
        print(f"No file found at {bam}")
        return []


def biomart_query(ids, filter_name, attributes_list, batch_size) -> pd.DataFrame:
    """Query BioMart in batches for accessions for the given ENST IDs."""
    if not ids:
        return pd.DataFrame()

    print(
        f"\n2. Connecting to BioMart and querying in batches of {batch_size}...")
    all_results = []

    try:
        database = Dataset(name='hsapiens_gene_ensembl',
                           host='http://www.ensembl.org')

        attributes = attributes_list

        for i in range(0, len(ids), batch_size):
            batch = ids[i:i + batch_size]
            print(
                f"Querying batch {i // batch_size + 1} of {(len(ids) // batch_size) + 1} ({len(batch)} IDs)...")

            filters = {filter_name: batch}
            map_df = database.query(
                attributes=attributes,
                filters=filters,
                only_unique=True)

            print(
                f"   Returned columns: {list(map_df.columns)} | Rows: {len(map_df)}")
            all_results.append(map_df)

        final_df = pd.concat(all_results, ignore_index=True)

        final_df.columns = [c.strip().replace(" ", "_").lower()
                            for c in final_df.columns]

        print("\n--- Final Combined DataFrame (normalized) ---")
        print(final_df.head(100))

        print(f"Total rows: {len(final_df)}\n")

        return final_df

    except Exception as e:
        print(f'\n--- BioMart Query FAILED ---')
        print(f'An exception occurred during batch processing: {e}')
        return pd.DataFrame()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="This tool takes bam files to read transcript ids and query"
        "biomart to translate to other ids or information. If you don't know the name of "
        "any filters or attributes try ensembl_gene_id or biomart.py --search_terms for a list of filters and "
        "attributes"
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    search_parser = subparsers.add_parser(
        'query', help='Run the biomart search')

    search_parser.add_argument(
        '--root_dir',
        '-r',
        type=str,
        help='Root directory where biomart will search for bam files recursively',
        required=True
    )

    search_parser.add_argument(
        '--filter',
        '-f',
        type=str,
        default='link_ensembl_transcript_stable_id',
        help='The input id (use link_ensembl_transcript_stable_id for ensts)'
    )

    search_parser.add_argument(
        '--attributes',
        '-a',
        nargs='+',
        help='One or a list of args for the information you want biomart to fetch (eg. try transcript_biotype)',
        required=True
    )

    search_parser.add_argument(
        '--output',
        '-o',
        type=str,
        required=True,
        help='Output file name (only csv supported for now)'
    )

    search_parser.add_argument(
        '--batch_size',
        '-b',
        type=int,
        required=False,
        help='The number of ids to be sent in each query request. Note that Biomart suggests sizes less than 500'
    )

    help_parser = subparsers.add_parser(
        'helper', help='utilities on using biomart search')

    help_parser.add_argument(
        '--search_term',
        '-s',
        type=str,
        help='Takes a search term and returns the related attributes and filters'
    )

    args = parser.parse_args()

    if '--helper' in sys.argv:
        if args.search_term:
            term_help(args.search_term)
            sys.exit()

    all_bam_files = find_bam_files(args.root_dir)

    if not all_bam_files:
        print(f"No .bam files found in {args.root_dir}. Exiting.")
    else:
        all_enst_ids = set()
        for i, bam_file in enumerate(all_bam_files):
            print(f"\n--- Processing file {i+1}/{len(all_bam_files)} ---")
            enst_list_from_file = parse_enst(bam_file)
            all_enst_ids.update(enst_list_from_file)

        print(f"\n--- Aggregation Complete ---")
        print(
            f"Found {len(all_enst_ids)} total unique ENST IDs across {len(all_bam_files)} files.")

        final_enst_list = list(all_enst_ids)
        biomart_df = biomart_query(
            final_enst_list, args.filter, args.attributes, batch_size=args.batch_size)

        print("Biomart Queries Complete :)")
        print(f"{len(biomart_df)} entries before drop")
        biomart_df = biomart_df.dropna()
        print(f"{len(biomart_df)} entries after drop")
        biomart_df.to_csv(args.output)
