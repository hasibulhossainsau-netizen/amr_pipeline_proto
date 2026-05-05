#!/usr/bin/env python3
"""
run_argannot_blast.py - BLAST against SRST2 ARG-ANNOT database
Using ARGannot_r3.fasta 
"""

import subprocess
import sys
import pandas as pd
from pathlib import Path
import yaml
import logging
from typing import Dict, List, Optional
import argparse
from datetime import datetime
import re

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ARGANNOTPipeline:
    def __init__(self, config_path: str = "config.yaml"):
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # BLAST settings
        argannot_config = self.config.get('argannot', {})
        blast_config = argannot_config.get('blast', {})
        
        self.blast_evalue = blast_config.get('evalue', 1e-30)
        self.max_hits = blast_config.get('max_hits', 1)
        self.threads = blast_config.get('threads', 4)
        self.min_identity = argannot_config.get('min_identity', 80)
        self.min_coverage = argannot_config.get('min_coverage', 70)  # kept but not used strictly
        
        # Paths
        self.assembly_dir = Path(self.config['paths']['assemblies'])
        self.results_dir = Path(self.config['paths']['results']) / "argannot"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # ARG-ANNOT database (SRST2 version)
        self.argannot_db = Path(self.config['databases']['argannot'])
        
        self.log_dir = Path(self.config['paths']['logs'])
        self.log_dir.mkdir(exist_ok=True)
        
        self._validate_database()
    
    def _validate_database(self):
        db_files = list(self.argannot_db.parent.glob(f"{self.argannot_db.name}.*"))
        if not db_files:
            logger.error(f"ARG-ANNOT BLAST database not found at {self.argannot_db}")
            logger.error("Please run: makeblastdb -in ARGannot_r3.fasta -dbtype nucl -out ARGannot_r3_db")
            sys.exit(1)
        logger.info(f"ARG-ANNOT DB (SRST2): {self.argannot_db}")
    
    def extract_gene_info(self, header: str) -> tuple:
        """
        Parse ARG-ANNOT header format:
        >(AGly) AadA1:M95287:3320-4111:792 | Cluster_1
        >(Bla) blaCTX-M-15:JF911294:191-1003:813 | Cluster_45
        >(Col) mcr-1.1:KP347127:1-1626:1626 | Cluster_156
        """
        # Extract antibiotic class
        class_match = re.search(r'\(([^)]+)\)', header)
        drug_class = class_match.group(1) if class_match else 'Unknown'
        
        # Map to standard drug classes
        class_mapping = {
            'AGly': 'Aminoglycoside',
            'Bla': 'Beta-lactam',
            'Col': 'Colistin',
            'Fos': 'Fosfomycin',
            'Flq': 'Fluoroquinolone',
            'Gly': 'Glycopeptide',
            'MLS': 'Macrolide',
            'Phe': 'Phenicol',
            'Rif': 'Rifampicin',
            'Sul': 'Sulfonamide',
            'Tet': 'Tetracycline',
            'Tmt': 'Trimethoprim'
        }
        drug_class = class_mapping.get(drug_class, drug_class)
        
        # Extract gene name
        gene_match = re.search(r'\)\s+([^:]+)', header)
        gene_name = gene_match.group(1) if gene_match else 'Unknown'
        
        # Clean gene name (remove suffixes)
        gene_name = re.sub(r'_\d+$', '', gene_name)  # Remove _1, _2 etc
        gene_name = gene_name.split('|')[0].strip()
        
        return gene_name, drug_class
    
    def parse_blast_output(self, blast_file: Path, sample_id: str, argannot_fasta: Path) -> List[Dict]:
        """Parse BLAST output with ARG-ANNOT header parsing"""
        results = []
    
        if not blast_file.exists() or blast_file.stat().st_size == 0:
            return results
    
        # Load ARG-ANNOT headers
        headers = {}
        with open(argannot_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip().lstrip('>')
                    gene_name, drug_class = self.extract_gene_info(header)
                    headers[header.split('|')[0].strip()] = (gene_name, drug_class)
    
        # Parse BLAST output (outfmt 6)
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                'evalue', 'bitscore']
    
        try:
            df = pd.read_csv(blast_file, sep='\t', names=columns, header=None)
            
            for _, row in df.iterrows():
                subject = row['sseqid'].split('|')[0].strip()
                
                if subject in headers:
                    gene_name, drug_class = headers[subject]
                else:
                    gene_name = subject
                    drug_class = 'Unknown'
            
                identity = row['pident']
            
                # Apply identity filter only (coverage optional)
                if identity >= self.min_identity:
                    result = {
                        'gene': gene_name,
                        'identity': identity,
                        'coverage': row['length'],  # alignment length
                        'contig': row['qseqid'],
                        'position_start': int(row['qstart']),
                        'position_end': int(row['qend']),
                        'drug_class': drug_class,
                        'evalue': row['evalue'],
                        'bitscore': row['bitscore'],
                        'database': 'ARG-ANNOT',
                        'detection_tool': 'BLASTn'
                    }
                    results.append(result)
        
            logger.info(f"Parsed {len(results)} hits for {sample_id}")
        
        except Exception as e:
            logger.error(f"Failed to parse BLAST output: {e}")
    
        return results
    
    def run_blast(self, sample_id: str, contigs_file: Path, argannot_fasta: Path) -> Optional[List[Dict]]:
        """Run BLASTn against ARG-ANNOT database"""
        
        if not contigs_file.exists():
            logger.error(f"Contigs not found: {contigs_file}")
            return None
        
        output_dir = self.results_dir / sample_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        blast_output = output_dir / "blast_results.txt"
        
        if blast_output.exists() and blast_output.stat().st_size > 0:
            logger.info(f"BLAST already run for {sample_id}")
            return self.parse_blast_output(blast_output, sample_id, argannot_fasta)
        
        cmd = [
            "blastn",
            "-query", str(contigs_file),
            "-db", str(self.argannot_db),
            "-out", str(blast_output),
            "-outfmt", "6",
            "-evalue", str(self.blast_evalue),
            "-max_target_seqs", str(self.max_hits),
            "-num_threads", str(self.threads),
            "-perc_identity", str(self.min_identity)
        ]
        
        logger.info(f"Running BLAST for {sample_id}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            
            # Log output
            log_file = self.log_dir / f"{sample_id}_argannot.log"
            with open(log_file, 'w') as f:
                f.write(' '.join(cmd) + '\n\n')
                f.write(result.stdout)
                f.write(result.stderr)
            
            if result.returncode != 0:
                logger.error(f"BLAST failed for {sample_id}")
                return None
            
            return self.parse_blast_output(blast_output, sample_id, argannot_fasta)
            
        except subprocess.TimeoutExpired:
            logger.error(f"BLAST timed out for {sample_id}")
            return None
        except Exception as e:
            logger.error(f"Error for {sample_id}: {e}")
            return None
    
    def process_sample(self, sample_id: str, force: bool = False) -> Optional[pd.DataFrame]:
        contigs_file = self.assembly_dir / sample_id / "contigs.fasta"
        # If contigs don't exist, return an empty DataFrame with correct columns
        if not contigs_file.exists():
            logger.warning(f"Contigs file not found for {sample_id}, returning empty result.")
            return pd.DataFrame(columns=[
                'sample_id', 'gene', 'identity', 'coverage', 'contig', 
                'position_start', 'position_end', 'drug_class', 'evalue', 
                'bitscore', 'database', 'detection_tool', 'type'
            ])
        
        argannot_fasta = Path(self.config['databases']['argannot_fasta'])
        
        results = self.run_blast(sample_id, contigs_file, argannot_fasta)
        
        if results is None:
            # BLAST failed or no results, return empty DataFrame
            return pd.DataFrame(columns=[
                'sample_id', 'gene', 'identity', 'coverage', 'contig', 
                'position_start', 'position_end', 'drug_class', 'evalue', 
                'bitscore', 'database', 'detection_tool', 'type'
            ])
        
        df = pd.DataFrame(results)
        if not df.empty:
            df.insert(0, 'sample_id', sample_id)
            df['type'] = 'acquired_gene'
        else:
            # If results list empty, provide an empty DataFrame with columns
            df = pd.DataFrame(columns=[
                'sample_id', 'gene', 'identity', 'coverage', 'contig', 
                'position_start', 'position_end', 'drug_class', 'evalue', 
                'bitscore', 'database', 'detection_tool', 'type'
            ])
        return df

def save_results(all_results: Dict[str, pd.DataFrame], config: Dict):
    if not all_results:
        return
    
    dfs = [df for df in all_results.values() if not df.empty]
    if not dfs:
        return
    
    combined_df = pd.concat(dfs, ignore_index=True)
    
    output_file = Path(config['paths']['results']) / "argannot_blast_results.csv"
    combined_df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(combined_df)} ARG-ANNOT hits to {output_file}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', help='Sample ID')
    parser.add_argument('--manifest', help='Manifest CSV')
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--output', help='Single sample output TSV file (for Snakemake)')
    args = parser.parse_args()
    
    pipeline = ARGANNOTPipeline(args.config)
    
    # If --output is given, we must be in single-sample mode
    if args.output:
        if not args.sample:
            logger.error("--output requires --sample")
            sys.exit(1)
        df = pipeline.process_sample(args.sample, force=args.force)
        # Write TSV (even if empty, it will have headers)
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Single sample result written to {output_path}")
        return

    # Original batch mode (--manifest, --sample without --output, or directory scan)
    if args.sample:
        samples = [args.sample]
    elif args.manifest:
        df = pd.read_csv(args.manifest)
        samples = df['sample_id'].tolist()
    else:
        samples = [d.name for d in pipeline.assembly_dir.iterdir() if d.is_dir()]
    
    all_results = {}
    for idx, sample_id in enumerate(samples, 1):
        logger.info(f"Processing {idx}/{len(samples)}: {sample_id}")
        df = pipeline.process_sample(sample_id, force=args.force)
        if df is not None and not df.empty:
            all_results[sample_id] = df
            logger.info(f"  ✓ {len(df)} hits")
    
    if all_results:
        save_results(all_results, pipeline.config)
    
    logger.info(f"Complete: {len(all_results)} samples processed")

if __name__ == "__main__":
    main()