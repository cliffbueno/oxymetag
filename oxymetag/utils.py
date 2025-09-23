#!/usr/bin/env python3
"""
Utility functions for OxyMetaG
"""

import subprocess
import pkg_resources
from pathlib import Path
import logging

logger = logging.getLogger('oxymetag')

class OxyMetaGError(Exception):
    """Custom exception for OxyMetaG errors"""
    pass

def check_dependencies():
    """Check if required external tools are available"""
    required_tools = ['kraken2', 'diamond', 'Rscript']
    missing_tools = []
    
    for tool in required_tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        raise OxyMetaGError(f"Missing required tools: {', '.join(missing_tools)}")

def get_package_data_path(filename: str) -> str:
    """Get path to package data files"""
    try:
        return pkg_resources.resource_filename('oxymetag', f'data/{filename}')
    except:
        package_dir = Path(__file__).parent
        return str(package_dir / 'data' / filename)

def run_kraken2_setup():
    """Download and set up standard Kraken2 database"""
    logger.info("Setting up Kraken2 standard database...")
    
    db_path = Path.cwd() / "kraken2_db"
    db_path.mkdir(exist_ok=True)
    
    try:
        cmd = ['kraken2-build', '--download-taxonomy', '--db', str(db_path)]
        subprocess.run(cmd, check=True)
        logger.info("Taxonomy downloaded successfully")
        
        libraries = ['bacteria', 'archaea', 'viral', 'fungi']
        for lib in libraries:
            cmd = ['kraken2-build', '--download-library', lib, '--db', str(db_path)]
            logger.info(f"Downloading {lib} library...")
            subprocess.run(cmd, check=True)
        
        cmd = ['kraken2-build', '--build', '--db', str(db_path), '--threads', '4']
        logger.info("Building Kraken2 database...")
        subprocess.run(cmd, check=True)
        
        logger.info(f"Kraken2 database setup complete: {db_path}")
        
    except subprocess.CalledProcessError as e:
        raise OxyMetaGError(f"Failed to setup Kraken2 database: {e}")
