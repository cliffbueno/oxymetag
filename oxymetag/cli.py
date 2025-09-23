#!/usr/bin/env python3
"""
Command line interface for OxyMetaG
"""

import sys
import argparse
import logging

from . import __version__
from .core import extract_reads, profile_samples, predict_aerobes
from .utils import check_dependencies, run_kraken2_setup, OxyMetaGError

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('oxymetag')


def main():
    """Main CLI interface"""
    parser = argparse.ArgumentParser(
        description="OxyMetaG: Oxygen metabolism profiling from metagenomic data",
        prog="oxymetag"
    )
    parser.add_argument('--version', action='version', version=f'OxyMetaG {__version__}')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Setup command
    setup_parser = subparsers.add_parser('setup', help='Setup Kraken2 database')
    
    # Extract command
    extract_parser = subparsers.add_parser('extract', help='Extract bacterial reads')
    extract_parser.add_argument('-i', '--input', nargs='+', required=True,
                               help='Input fastq.gz files')
    extract_parser.add_argument('-o', '--output', default='BactReads',
                               help='Output directory (default: BactReads)')
    extract_parser.add_argument('-t', '--threads', type=int, default=4,
                               help='Number of threads (default: 4)')
    extract_parser.add_argument('--kraken-db', default='kraken2_db',
                               help='Kraken2 database path (default: kraken2_db)')
    
    # Profile command  
    profile_parser = subparsers.add_parser('profile', help='Profile samples with DIAMOND')
    profile_parser.add_argument('-i', '--input', default='BactReads',
                               help='Input directory (default: BactReads)')
    profile_parser.add_argument('-o', '--output', default='diamond_output',
                               help='Output directory (default: diamond_output)')
    profile_parser.add_argument('-t', '--threads', type=int, default=4,
                               help='Number of threads (default: 4)')
    profile_parser.add_argument('--diamond-db',
                               help='DIAMOND database path (default: package database)')
    
    # Predict command
    predict_parser = subparsers.add_parser('predict', help='Predict aerobe levels')
    predict_parser.add_argument('-i', '--input', default='diamond_output',
                               help='Input directory (default: diamond_output)')
    predict_parser.add_argument('-o', '--output', default='oxygen_predictions.tsv',
                               help='Output file (default: oxygen_predictions.tsv)')
    predict_parser.add_argument('-t', '--threads', type=int, default=4,
                               help='Number of threads (default: 4)')
    predict_parser.add_argument('-m', '--mode', choices=['modern', 'ancient', 'custom'],
                               default='modern', help='Filtering mode (default: modern)')
    predict_parser.add_argument('--idcut', type=float,
                               help='Custom identity cutoff (for custom mode)')
    predict_parser.add_argument('--bitcut', type=float,
                               help='Custom bitscore cutoff (for custom mode)')
    predict_parser.add_argument('--ecut', type=float,
                               help='Custom e-value cutoff (for custom mode)')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        check_dependencies()
        
        if args.command == 'setup':
            run_kraken2_setup()
            
        elif args.command == 'extract':
            extract_reads(args.input, args.output, args.threads, args.kraken_db)
            
        elif args.command == 'profile':
            profile_samples(args.input, args.output, args.threads, args.diamond_db)
            
        elif args.command == 'predict':
            predict_aerobes(args.input, args.output, args.mode,
                           args.idcut, args.bitcut, args.ecut, args.threads)
        
        logger.info("Command completed successfully")
        
    except OxyMetaGError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
