#!/usr/bin/env python3
"""
run_resfinder.py - Run ResFinder v4.x with PointFinder on assembled contigs
Compatible with ResFinder 4.7.2
Debugged for Snakemake rule compatibility:
- Single-sample mode with --output TSV
- Removed batch processing, manifest updates
"""

import subprocess
import sys
import pandas as pd
from pathlib import Path
import yaml
import json
import logging
from typing import Dict, List, Optional, Set
import argparse
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ResFinderPipeline:
    def __init__(self, config_path: str = "config.yaml"):
        """Initialize ResFinder pipeline with configuration"""
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)
        
        # ResFinder settings - Convert percentage to decimal (0-1 range)
        self.min_identity = self.config['resfinder']['min_identity'] / 100.0
        self.min_coverage = self.config['resfinder']['min_coverage'] / 100.0
        self.min_length = self.config['resfinder']['min_length']
        self.threads = self.config['resfinder']['threads']
        self.point_mutations = self.config['resfinder']['point_mutations']
        self.species = self.config['resfinder']['species']
        
        # Paths
        self.assembly_dir = Path(self.config['paths']['assemblies'])
        self.results_dir = Path(self.config['paths']['results']) / "resfinder"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Database paths from config
        self.resfinder_db = Path(self.config['databases']['resfinder'])
        self.pointfinder_db = Path(self.config['databases']['pointfinder'])
        
        self.log_dir = Path(self.config['paths']['logs'])
        self.log_dir.mkdir(exist_ok=True)
        
        # Validate databases
        self._validate_databases()
    
    def _validate_databases(self):
        """Validate database paths before running"""
        if not self.resfinder_db.exists():
            logger.error(f"ResFinder database not found: {self.resfinder_db}")
            logger.error("Please check path in config.yaml")
            sys.exit(1)
        
        if self.point_mutations and not self.pointfinder_db.exists():
            logger.warning(f"PointFinder database not found: {self.pointfinder_db}")
            logger.warning("Point mutation detection will be disabled")
            self.point_mutations = False
        
        logger.info(f"ResFinder DB: {self.resfinder_db}")
        logger.info(f"PointFinder DB: {self.pointfinder_db if self.point_mutations else 'Disabled'}")
    
    def parse_resfinder_json(self, json_file: Path, sample_id: str) -> pd.DataFrame:
        """
        Parse ResFinder/PointFinder JSON (compatible with your structure)
        Properly extracts gene names and distinguishes between acquired genes and point mutations
        """
        
        if not json_file.exists():
            logger.debug(f"JSON file not found: {json_file}")
            return pd.DataFrame()
        
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Handle different JSON structures
            if 'seq_regions' in data:
                seq_regions = data['seq_regions']
            elif 'resfinder' in data and 'results' in data['resfinder']:
                seq_regions = data['resfinder']['results'].get('seq_regions', {})
            elif 'results' in data and 'seq_regions' in data['results']:
                seq_regions = data['results']['seq_regions']
            else:
                logger.warning(f"Unexpected JSON structure in {json_file}")
                logger.debug(f"JSON keys: {data.keys()}")
                return pd.DataFrame()
            
            results = []
            
            for region_id, region in seq_regions.items():
                logger.debug(f"Processing region: {region_id}")
                logger.debug(f"Region data: {region}")
                
                #  Try multiple ways to get gene name
                gene = None
                
                # Method 1: Check 'name' field
                if 'name' in region and region['name']:
                    gene = region['name']
                # Method 2: Check 'gene' field
                elif 'gene' in region and region['gene']:
                    gene = region['gene']
                # Method 3: Parse from region_id
                elif ';;' in region_id:
                    gene = region_id.split(';;')[0]
                # Method 4: Use region_id as is
                else:
                    gene = region_id
                
                # Skip if no gene name
                if not gene or gene == '':
                    logger.warning(f"No gene name found for region: {region_id}")
                    continue
                
                # Clean up gene name (remove extra whitespace)
                gene = gene.strip()
                
                #  Detect source (ResFinder vs PointFinder)
                database = region.get('ref_database', [])
                if isinstance(database, list) and database:
                    database = database[0]
                elif isinstance(database, str):
                    database = database
                else:
                    database = 'Unknown'
                
                if 'ResFinder' in str(database):
                    result_type = 'acquired_gene'
                elif 'PointFinder' in str(database):
                    result_type = 'point_mutation'
                else:
                    result_type = 'unknown'
                
                #  Extract values safely
                identity = float(region.get('identity', 0.0))
                coverage = float(region.get('coverage', 0.0))
                
                # Skip garbage low coverage hits (important for 23S spam)
                if coverage < 60.0 and result_type == 'acquired_gene':
                    logger.debug(f"Skipping {gene} due to low coverage: {coverage}")
                    continue
                
                # Extract phenotypes/drug classes
                phenotypes = region.get('phenotypes', [])
                if isinstance(phenotypes, list):
                    drug_class = '; '.join(phenotypes) if phenotypes else ''
                else:
                    drug_class = str(phenotypes) if phenotypes else ''
                
                # Extract position information
                query_start = region.get('query_start_pos', 0)
                query_end = region.get('query_end_pos', 0)
                contig = region.get('query_id', '')
                
                # For point mutations, extract mutation details
                mutation = ''
                if result_type == 'point_mutation':
                    # Try different fields for mutation info
                    mutation = region.get('mutation', '')
                    if not mutation and 'aa_mutation' in region:
                        mutation = region.get('aa_mutation', '')
                    if mutation:
                        gene = f"{gene}_{mutation}"
                
                result = {
                    'sample_id': sample_id,
                    'gene': gene,
                    'identity': identity,
                    'coverage': coverage,
                    'contig': contig,
                    'position_start': query_start,
                    'position_end': query_end,
                    'drug_class': drug_class,
                    'database': str(database),
                    'type': result_type,
                    'mutation': str(mutation) if mutation else ''
                }
                
                results.append(result)
                logger.debug(f"Added result: {gene} ({result_type})")
            
            if results:
                logger.info(f"Parsed {len(results)} hits from JSON for {sample_id} ({sum(1 for r in results if r['type']=='acquired_gene')} genes, {sum(1 for r in results if r['type']=='point_mutation')} mutations)")
            else:
                logger.warning(f"No valid hits found in JSON for {sample_id}")
            
            if not results:
                return pd.DataFrame()
            
            return pd.DataFrame(results)
            
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON for {sample_id}: {e}")
            return pd.DataFrame()
        except Exception as e:
            logger.error(f"Unexpected error parsing JSON for {sample_id}: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return pd.DataFrame()
    
    def parse_resfinder_output(self, output_dir: Path, sample_id: str) -> List[Dict]:
        """Parse ResFinder v4.x output - primary method uses JSON"""
        results = []
        
        # Try JSON first (preferred method)
        json_file = output_dir / "ResFinder_results.json"
        if json_file.exists():
            logger.info(f"Found JSON file for {sample_id}, attempting to parse...")
            df = self.parse_resfinder_json(json_file, sample_id)
            if not df.empty:
                logger.info(f"Successfully parsed {len(df)} results from JSON")
                return df.to_dict(orient='records')
            else:
                logger.warning(f"JSON parsing returned empty DataFrame for {sample_id}")
        
        # Fallback to TSV if JSON parsing failed or no JSON
        tsv_file = output_dir / "ResFinder_results_tab.txt"
        if tsv_file.exists():
            logger.info(f"Falling back to TSV parsing for {sample_id}")
            try:
                # Read TSV with proper handling
                df = pd.read_csv(tsv_file, sep='\t')
                logger.debug(f"TSV columns: {df.columns.tolist()}")
                logger.debug(f"TSV shape: {df.shape}")
                
                for idx, row in df.iterrows():
                    # Try different column names for gene
                    gene = ''
                    for col in ['Gene', 'gene', 'Name', 'name', 'Resistance gene']:
                        if col in df.columns and pd.notna(row[col]):
                            gene = str(row[col]).strip()
                            break
                    
                    if not gene:
                        logger.warning(f"No gene name found in row {idx}")
                        continue
                    
                    result = {
                        'sample_id': sample_id,
                        'gene': gene,
                        'identity': float(row.get('Identity', row.get('identity', 0))),
                        'coverage': float(row.get('Coverage', row.get('coverage', 0))),
                        'contig': str(row.get('Contig', row.get('contig', ''))),
                        'position_start': int(row.get('Position start', row.get('position_start', 0))),
                        'position_end': int(row.get('Position end', row.get('position_end', 0))),
                        'drug_class': str(row.get('Phenotype', row.get('phenotype', 'Unknown'))),
                        'database': 'ResFinder',
                        'type': 'acquired_gene',
                        'mutation': '',
                    }
                    results.append(result)
                
                logger.info(f"Parsed {len(results)} AMR genes from TSV for {sample_id}")
                
            except Exception as e:
                logger.error(f"Failed to parse TSV: {e}")
                import traceback
                logger.debug(traceback.format_exc())
        
        return results
    
    def parse_pointfinder_output(self, output_dir: Path, sample_id: str) -> List[Dict]:
        """Parse PointFinder output for point mutations - uses JSON if available"""
        results = []
        
        # Try JSON first
        json_file = output_dir / "ResFinder_results.json"
        if json_file.exists():
            df = self.parse_resfinder_json(json_file, sample_id)
            if not df.empty:
                # Filter for point mutations only
                point_df = df[df['type'] == 'point_mutation']
                if not point_df.empty:
                    logger.info(f"Found {len(point_df)} point mutations from JSON")
                    return point_df.to_dict(orient='records')
        
        # Fallback to legacy PointFinder_results.txt
        point_file = output_dir / "PointFinder_results.txt"
        if not point_file.exists():
            return results
        
        logger.info(f"Falling back to legacy PointFinder output for {sample_id}")
        try:
            with open(point_file, 'r') as f:
                content = f.read()
            
            lines = content.strip().split('\n')
            if len(lines) > 1:
                headers = lines[0].split('\t')
                for line in lines[1:]:
                    values = line.split('\t')
                    if len(values) >= len(headers):
                        row = dict(zip(headers, values))
                        
                        # Extract gene name
                        gene = row.get('Gene', '')
                        mutation = row.get('Mutation', '')
                        
                        if not gene:
                            continue
                        
                        result = {
                            'sample_id': sample_id,
                            'gene': gene,
                            'mutation': mutation,
                            'identity': float(row.get('Identity', 0)),
                            'coverage': float(row.get('Coverage', 0)),
                            'drug_class': row.get('Phenotype', 'Unknown'),
                            'database': 'PointFinder',
                            'type': 'point_mutation',
                            'contig': '',
                            'position_start': 0,
                            'position_end': 0,
                        }
                        # Append mutation to gene name for clarity
                        if mutation:
                            result['gene'] = f"{gene}_{mutation}"
                        results.append(result)
            
            if results:
                logger.info(f"Parsed {len(results)} point mutations from legacy PointFinder output for {sample_id}")
                
        except Exception as e:
            logger.error(f"Failed to parse PointFinder output: {e}")
        
        return results
    
    def run_resfinder(self, sample_id: str, contigs_file: Path) -> Optional[Dict]:
        """Run ResFinder v4.x with PointFinder on assembled contigs"""
        
        if not contigs_file.exists():
            logger.error(f"Contigs file not found: {contigs_file}")
            return None
        
        output_dir = self.results_dir / sample_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if output already exists
        json_file = output_dir / "ResFinder_results.json"
        if json_file.exists() and json_file.stat().st_size > 100:
            logger.info(f"Valid ResFinder output already exists for {sample_id}")
            amr_results = self.parse_resfinder_output(output_dir, sample_id)
            point_results = self.parse_pointfinder_output(output_dir, sample_id)
            
            return {
                'amr_genes': amr_results,
                'point_mutations': point_results,
                'success': True
            }
        
        # Build command with correct flags
        cmd = [
            sys.executable, "-m", "resfinder",
            "-ifa", str(contigs_file),
            "-o", str(output_dir),
            "-db_res", str(self.resfinder_db),
            "-acq",
            "-l", str(self.min_coverage),
            "-t", str(self.min_identity),
            "-s", self.species
        ]
        
        # Add PointFinder if enabled
        if self.point_mutations:
            cmd.extend([
                "-c",
                "-db_point", str(self.pointfinder_db)
            ])
        
        logger.info(f"Running ResFinder for {sample_id}")
        logger.debug(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600
            )
            
            # Log output
            log_file = self.log_dir / f"{sample_id}_resfinder.log"
            with open(log_file, 'w') as f:
                f.write("=== COMMAND ===\n")
                f.write(' '.join(cmd))
                f.write("\n\n=== STDOUT ===\n")
                f.write(result.stdout)
                f.write("\n\n=== STDERR ===\n")
                f.write(result.stderr)
            
            if result.returncode != 0:
                logger.error(f"ResFinder failed for {sample_id} (exit code {result.returncode})")
                return None
            
            # Parse results using improved JSON parser
            amr_results = self.parse_resfinder_output(output_dir, sample_id)
            point_results = self.parse_pointfinder_output(output_dir, sample_id)
            
            if amr_results:
                logger.info(f"Found {len(amr_results)} AMR genes in {sample_id}")
                # Log first few genes
                genes_list = [r.get('gene', '') for r in amr_results[:5]]
                logger.info(f"  Genes found: {', '.join(genes_list)}")
            else:
                logger.warning(f"No AMR genes found in {sample_id}")
            
            if point_results:
                logger.info(f"Found {len(point_results)} point mutations in {sample_id}")
                mutations_list = [r.get('gene', '') for r in point_results[:5]]
                logger.info(f"  Mutations: {', '.join(mutations_list)}")
            
            return {
                'amr_genes': amr_results,
                'point_mutations': point_results,
                'success': True
            }
            
        except subprocess.TimeoutExpired:
            logger.error(f"ResFinder timed out for {sample_id}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error for {sample_id}: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def process_sample(self, sample_id: str, force: bool = False) -> Optional[pd.DataFrame]:
        """Process a single sample through ResFinder"""
        
        # Find contigs file
        contigs_file = self.assembly_dir / sample_id / "contigs.fasta"
        if not contigs_file.exists():
            logger.warning(f"Assembly not found for {sample_id}")
            return None
        
        # Run ResFinder
        results = self.run_resfinder(sample_id, contigs_file)
        
        if not results or not results['success']:
            return None
        
        # Combine results
        all_results = []
        
        for gene in results.get('amr_genes', []):
            if 'sample_id' not in gene:
                gene['sample_id'] = sample_id
            all_results.append(gene)
        
        for mutation in results.get('point_mutations', []):
            if 'sample_id' not in mutation:
                mutation['sample_id'] = sample_id
            all_results.append(mutation)
        
        if not all_results:
            logger.info(f"No AMR determinants found for {sample_id}")
            return pd.DataFrame()
        
        df = pd.DataFrame(all_results)
        logger.info(f"Created DataFrame with {len(df)} rows for {sample_id}")
        logger.debug(f"DataFrame head:\n{df[['gene', 'type']].head()}")
        
        return df

def main():
    parser = argparse.ArgumentParser(description='Run ResFinder on assembled contigs')
    parser.add_argument('--sample', help='Specific sample ID to process')
    parser.add_argument('--output', help='Output TSV file path (required for Snakemake single-sample mode)')
    parser.add_argument('--manifest', help='Input manifest CSV with sample IDs (for batch processing)')
    parser.add_argument('--config', default='config.yaml', help='Config file path')
    parser.add_argument('--force', action='store_true', help='Force reprocessing even if results exist')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    pipeline = ResFinderPipeline(args.config)
    
    # Single-sample mode (Snakemake integration)
    if args.sample and args.output:
        logger.info(f"Running in single-sample mode for {args.sample}")
        sample_df = pipeline.process_sample(args.sample, force=args.force)
        
        # Create output directory if needed
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Predefined column order for TSV (consistent with downstream tools)
        columns_order = ['sample_id', 'gene', 'type', 'identity', 'coverage',
                         'drug_class', 'database', 'contig', 'position_start',
                         'position_end', 'mutation']
        
        if sample_df is not None:
            # Ensure all columns present, fill missing
            for col in columns_order:
                if col not in sample_df.columns:
                    sample_df[col] = ''
            # Reorder and save as TSV
            sample_df = sample_df[columns_order]
            sample_df.to_csv(output_path, sep='\t', index=False)
            logger.info(f"Written {len(sample_df)} hits to {output_path}")
        else:
            # Write empty dataframe with headers
            empty_df = pd.DataFrame(columns=columns_order)
            empty_df.to_csv(output_path, sep='\t', index=False)
            logger.warning(f"No results, written empty TSV to {output_path}")
        
        sys.exit(0)
    
    # Collect samples to process
    if args.sample:
        samples = [args.sample]
    elif args.manifest:
        df = pd.read_csv(args.manifest)
        if 'assembly_success' in df.columns:
            df = df[df['assembly_success'] == True]
        samples = df['sample_id'].tolist()
        logger.info(f"Loaded {len(samples)} assembled samples from {args.manifest}")
    else:
        assembly_dirs = [d for d in pipeline.assembly_dir.iterdir() if d.is_dir()]
        samples = [d.name for d in assembly_dirs]
        logger.info(f"Found {len(samples)} assembled samples")
    
    if not samples:
        logger.error("No samples to process!")
        sys.exit(1)
    
    logger.info(f"Processing {len(samples)} samples in batch mode...")
    
    # Process each sample
    all_results = {}
    failed_samples = []
    
    for idx, sample_id in enumerate(samples, 1):
        logger.info(f"Processing {idx}/{len(samples)}: {sample_id}")
        sample_df = pipeline.process_sample(sample_id, force=args.force)
        if sample_df is not None:
            all_results[sample_id] = sample_df
            if not sample_df.empty:
                logger.info(f"✓ Completed for {sample_id} ({len(sample_df)} hits)")
                genes_found = sample_df['gene'].tolist()[:5]
                logger.info(f"  Genes: {', '.join(genes_found)}{'...' if len(sample_df) > 5 else ''}")
            else:
                logger.info(f"✓ Completed for {sample_id} (no hits)")
        else:
            failed_samples.append(sample_id)
            logger.error(f"✗ Failed for {sample_id}")
    
    # Save combined results for batch mode 
    if all_results:
        # Combine into a single CSV for manual inspection, but do NOT update manifest automatically
        combined_df = pd.concat([df for df in all_results.values() if not df.empty], ignore_index=True)
        detailed_file = Path(pipeline.config['paths']['results']) / "resfinder_detailed_results.csv"
        combined_df.to_csv(detailed_file, index=False)
        logger.info(f"Batch results saved to {detailed_file}")
    
    # Summary
    logger.info("=" * 50)
    logger.info(f"Complete: {len(all_results)} succeeded, {len(failed_samples)} failed")
    if failed_samples:
        logger.info(f"Failed: {', '.join(failed_samples[:5])}")
    total_hits = sum(len(df) for df in all_results.values() if not df.empty)
    logger.info(f"Total AMR determinants found: {total_hits}")
    logger.info("=" * 50)

if __name__ == "__main__":
    main()