#!/usr/bin/env python3
"""
Core functions for OxyMetaG
"""

import os
import sys
import subprocess
import glob
from pathlib import Path
import pandas as pd
import logging
from typing import List, Optional

logger = logging.getLogger('oxymetag')

class OxyMetaGError(Exception):
    """Custom exception for OxyMetaG errors"""
    pass

def extract_reads(input_files: List[str], output_dir: str = "BactReads", 
                 threads: int = 48, kraken_db: str = "kraken2_db"):
    """
    Extract bacterial reads from metagenomic samples using Kraken2
    """
    logger.info(f"Starting bacterial read extraction with {threads} threads")
    
    if not Path(kraken_db).exists():
        raise OxyMetaGError(f"Kraken2 database not found: {kraken_db}")
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    for input_file in input_files:
        input_path = Path(input_file)
        if not input_path.exists():
            logger.warning(f"Input file not found: {input_file}")
            continue
            
        logger.info(f"Processing {input_file}")
        
        base_name = input_path.stem.replace('.fastq', '').replace('.gz', '')
        
        if '_R1' in base_name or '_1' in base_name:
            kraken_base = base_name.replace('_R1', '').replace('_1', '')
        else:
            kraken_base = base_name
        
        kraken_output = output_path / f"{kraken_base}_kraken.out"
        kraken_report = output_path / f"{kraken_base}_report.txt"
        
        if '_R1' in base_name or '_1' in base_name:
            r2_file = str(input_path).replace('_R1', '_R2').replace('_1', '_2')
            if Path(r2_file).exists():
                cmd = [
                    'kraken2', '--db', kraken_db, '--threads', str(threads),
                    '--output', str(kraken_output), '--report', str(kraken_report),
                    '--paired', str(input_path), r2_file
                ]
            else:
                logger.warning(f"R2 file not found for {input_file}, treating as single-end")
                cmd = [
                    'kraken2', '--db', kraken_db, '--threads', str(threads),
                    '--output', str(kraken_output), '--report', str(kraken_report),
                    str(input_path)
                ]
        else:
            cmd = [
                'kraken2', '--db', kraken_db, '--threads', str(threads),
                '--output', str(kraken_output), '--report', str(kraken_report),
                str(input_path)
            ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"Kraken2 classification completed for {input_file}")
            
            bacterial_reads = output_path / f"{base_name}_bacterial.fastq"
            
            cmd = [
                'extract_kraken_reads.py',
                '-k', str(kraken_output),
                '-s', str(input_path),
                '-o', str(bacterial_reads),
                '-r', str(kraken_report),
                '--taxid', '2',
                '--include-children'
            ]
            
            if ('_R1' in base_name or '_1' in base_name) and Path(r2_file).exists():
                if '_R1' in base_name:
                    r2_output = output_path / f"{base_name.replace('_R1', '_R2')}_bacterial.fastq"
                else:
                    r2_output = output_path / f"{base_name.replace('_1', '_2')}_bacterial.fastq"
                cmd.extend(['-s2', r2_file])
                cmd.extend(['-o2', str(r2_output)])
            
            subprocess.run(cmd, check=True)
            subprocess.run(['gzip', str(bacterial_reads)], check=True)
            
            if ('_R1' in base_name or '_1' in base_name) and Path(r2_file).exists():
                subprocess.run(['gzip', str(r2_output)], check=True)
            
            logger.info(f"Bacterial reads extracted for {input_file}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {input_file}: {e}")
            continue


def profile_samples(input_dir: str = "BactReads", output_dir: str = None, 
                   threads: int = 4, method: str = "diamond",
                   diamond_db: str = None, mmseqs_db: str = None):
    """
    Profile samples using DIAMOND or MMseqs2 against Pfam database
    """
    from .utils import get_package_data_path
    
    logger.info(f"Starting sample profiling with {method} using {threads} threads")
    
    if output_dir is None:
        output_dir = 'diamond_output' if method == 'diamond' else 'mmseqs_output'
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    input_path = Path(input_dir)
    input_files = []
    
    patterns = [
        '*_R1_bacterial.fastq.gz',
        '*_1_bacterial.fastq.gz',
        '*_bacterial_R1.fastq.gz',
        '*_bacterial_1.fastq.gz',
        '*_bacterial.fastq.gz'
    ]
    
    for pattern in patterns:
        found_files = list(input_path.glob(pattern))
        if found_files:
            input_files.extend(found_files)
            logger.info(f"Found {len(found_files)} files using pattern: {pattern}")
            break
    
    if not input_files:
        all_files = list(input_path.glob("*.fastq.gz"))
        logger.error(f"FASTQ files in {input_dir}: {[f.name for f in all_files[:5]]}")
        raise OxyMetaGError(f"No bacterial read files found in {input_dir}")
    
    if method == 'diamond':
        _profile_with_diamond(input_files, output_path, threads, diamond_db)
    elif method == 'mmseqs2':
        _profile_with_mmseqs(input_files, output_path, threads, mmseqs_db)
    else:
        raise OxyMetaGError(f"Unknown method: {method}")


def _profile_with_diamond(input_files: List[Path], output_path: Path, 
                         threads: int, diamond_db: str = None):
    """Profile samples using DIAMOND blastx"""
    from .utils import get_package_data_path
    
    if diamond_db is None:
        diamond_db = get_package_data_path("oxymetag_pfams.dmnd")
    
    if not Path(diamond_db).exists():
        raise OxyMetaGError(f"DIAMOND database not found: {diamond_db}")
    
    for input_file in input_files:
        base_name = input_file.stem.replace('.fastq', '').replace('.gz', '')
        base_name = (base_name.replace('_R1_bacterial', '')
                             .replace('_1_bacterial', '')
                             .replace('_bacterial_R1', '')
                             .replace('_bacterial_1', '')
                             .replace('_bacterial', ''))
        
        logger.info(f"Processing {input_file}")
        
        output_file = output_path / f"{base_name}_diamond.tsv"
        
        cmd = [
            'diamond', 'blastx',
            '-d', diamond_db,
            '-q', str(input_file),
            '-o', str(output_file),
            '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"DIAMOND profiling completed for {input_file}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {input_file}: {e}")
            continue


def _profile_with_mmseqs(input_files: List[Path], output_path: Path,
                        threads: int, mmseqs_db: str = None):
    """Profile samples using MMseqs2 easy-search"""
    from .utils import get_package_data_path
    
    if mmseqs_db is None:
        mmseqs_db = get_package_data_path("oxymetag_pfams_n117_db")
    
    if not Path(mmseqs_db).exists():
        raise OxyMetaGError(f"MMseqs2 database not found: {mmseqs_db}")
    
    data_dir = Path(get_package_data_path(""))
    vtml_matrix = data_dir / "VTML20.out"
    nucl_matrix = data_dir / "nucleotide.out"
    
    if not vtml_matrix.exists():
        raise OxyMetaGError(f"VTML20.out matrix not found: {vtml_matrix}")
    if not nucl_matrix.exists():
        raise OxyMetaGError(f"nucleotide.out matrix not found: {nucl_matrix}")
    
    for input_file in input_files:
        base_name = input_file.stem.replace('.fastq', '').replace('.gz', '')
        base_name = (base_name.replace('_R1_bacterial', '')
                             .replace('_1_bacterial', '')
                             .replace('_bacterial_R1', '')
                             .replace('_bacterial_1', '')
                             .replace('_bacterial', ''))
        
        logger.info(f"Processing {input_file} with MMseqs2")
        
        output_file = output_path / f"{base_name}_mmseqs.tsv"
        tmp_dir = output_path / f"{base_name}_tmp"
        tmp_dir.mkdir(exist_ok=True)
        
        cmd = [
            'mmseqs', 'easy-search',
            str(input_file),
            str(mmseqs_db),
            str(output_file),
            str(tmp_dir),
            '--min-length', '12',
            '-e', '10.0',
            '--min-seq-id', '0.86',
            '-c', '0.65',
            '--cov-mode', '2',
            '--format-mode', '0',
            '--format-output', 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,cigar,qaln,taln',
            '--comp-bias-corr', '0',
            '--mask', '0',
            '--exact-kmer-matching', '1',
            '--sub-mat', f'aa:{vtml_matrix},nucl:{nucl_matrix}',
            '--seed-sub-mat', f'aa:{vtml_matrix},nucl:{nucl_matrix}',
            '-s', '2',
            '-k', '6',
            '--spaced-kmer-pattern', '11011101',
            '--max-seqs', '10000',
            '--max-rejected', '10',
            '--threads', str(threads),
            '--remove-tmp-files', '0',
            '--use-all-table-starts', '1'
        ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"MMseqs2 profiling completed for {input_file}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to process {input_file}: {e}")
            continue


def predict_aerobes(input_dir: str = None, output_file: str = "per_aerobe_predictions.tsv",
                   mode: str = "modern", id_cut: float = None, bit_cut: float = None, 
                   e_cut: float = None, threads: int = 4):
    """
    Predict aerobe levels from DIAMOND or MMseqs2 results
    Mode determines method: modern=DIAMOND, ancient=MMseqs2
    """
    from .utils import get_package_data_path
    
    logger.info(f"Starting aerobe level prediction in {mode} mode")
    
    if input_dir is None:
        input_dir = 'diamond_output' if mode == 'modern' else 'mmseqs_output'
    
    if mode == "modern":
        identity_cutoff, bitscore_cutoff, evalue_cutoff = 60.0, 50.0, 0.001
    elif mode == "ancient":
        identity_cutoff, bitscore_cutoff, evalue_cutoff = 86.0, 50.0, 0.001
    elif mode == "custom":
        if any(x is None for x in [id_cut, bit_cut, e_cut]):
            raise OxyMetaGError("Custom mode requires id_cut, bit_cut, and e_cut parameters")
        identity_cutoff, bitscore_cutoff, evalue_cutoff = id_cut, bit_cut, e_cut
    else:
        raise OxyMetaGError("Mode must be 'modern', 'ancient', or 'custom'")
    
    package_data_dir = str(Path(get_package_data_path("")).parent / "data")
    package_base = Path(__file__).parent
    r_script_path = str(package_base / "scripts" / "predict_oxygen.R")
    
    if not Path(input_dir).exists():
        raise OxyMetaGError(f"Input directory not found: {input_dir}")
    
    cmd = [
        'Rscript', r_script_path,
        input_dir, output_file, package_data_dir, mode,
        str(identity_cutoff), str(evalue_cutoff), str(bitscore_cutoff)
    ]
    
    try:
        logger.info(f"Calling R script: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("R script completed successfully")
        if result.stdout:
            logger.info(f"R output: {result.stdout}")
            
    except subprocess.CalledProcessError as e:
        logger.error(f"R script failed: {e}")
        if e.stderr:
            logger.error(f"R stderr: {e.stderr}")
        raise OxyMetaGError(f"Aerobe prediction failed: {e}")
    
    if Path(output_file).exists():
        results_df = pd.read_csv(output_file, sep='\t')
        logger.info(f"Results saved to {output_file}")
        return results_df
    else:
        raise OxyMetaGError(f"Output file was not created: {output_file}")
