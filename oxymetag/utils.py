#!/usr/bin/env python3
"""
Utility functions for OxyMetaG
"""

import subprocess
from pathlib import Path
import logging

# Use importlib.resources instead of deprecated pkg_resources
try:
    from importlib.resources import files
    use_importlib = True
except ImportError:
    use_importlib = False

logger = logging.getLogger('oxymetag')

class OxyMetaGError(Exception):
    """Custom exception for OxyMetaG errors"""
    pass

def get_package_data_path(filename: str) -> str:
    """Get path to package data files"""
    if use_importlib:
        try:
            package_files = files('oxymetag')
            if filename.startswith('../'):
                parts = filename.split('/')
                for part in parts:
                    if part == '..':
                        package_files = package_files.parent
                    elif part and part != '.':
                        package_files = package_files / part
                return str(package_files)
            else:
                return str(package_files / 'data' / filename)
        except:
            pass
    
    package_dir = Path(__file__).parent
    if filename.startswith('../'):
        return str(package_dir / filename)
    else:
        return str(package_dir / 'data' / filename)

def check_dependencies():
    """Check if required external tools are available"""
    required_tools = ['kraken2', 'diamond', 'mmseqs', 'Rscript']
    missing_tools = []
    
    for tool in required_tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        raise OxyMetaGError(f"Missing required tools: {', '.join(missing_tools)}")

def run_kraken2_setup():
    """Download and set up standard Kraken2 database without fungi"""
    logger.info("Setting up Kraken2 database (bacteria, archaea, viral)...")
    
    db_path = Path.cwd() / "kraken2_db"
    db_path.mkdir(exist_ok=True)
    
    try:
        cmd = ['kraken2-build', '--download-taxonomy', '--db', str(db_path)]
        logger.info("Downloading taxonomy...")
        subprocess.run(cmd, check=True)
        logger.info("Taxonomy downloaded successfully")
        
        libraries = ['bacteria', 'archaea', 'viral']
        for lib in libraries:
            cmd = ['kraken2-build', '--download-library', lib, '--db', str(db_path)]
            logger.info(f"Downloading {lib} library...")
            subprocess.run(cmd, check=True)
            logger.info(f"{lib} library downloaded successfully")
        
        cmd = ['kraken2-build', '--build', '--db', str(db_path), '--threads', '48']
        logger.info("Building Kraken2 database...")
        subprocess.run(cmd, check=True)
        
        cmd = ['kraken2-build', '--clean', '--db', str(db_path)]
        logger.info("Cleaning up temporary files...")
        subprocess.run(cmd, check=True)
        
        logger.info(f"Kraken2 database setup complete: {db_path}")
        logger.info("Database includes: bacteria, archaea, viral (fungi excluded)")
        
    except subprocess.CalledProcessError as e:
        raise OxyMetaGError(f"Failed to setup Kraken2 database: {e}")
