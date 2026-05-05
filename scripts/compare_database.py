#!/usr/bin/env python3
"""
compare_databases.py - Compare AMR gene detection across three databases
Produces Venn-diagram-style overlap table and visualizations
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import seaborn as sns
import logging
from typing import Dict, Set, List, Tuple
import argparse
from datetime import datetime
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class DatabaseComparator:
    def __init__(self, results_dir: str = "results", min_identity: float = 80):
        # Initialize comparator with results directory and minimum identity threshold
        self.results_dir = Path(results_dir)
        self.min_identity = min_identity
        
        # File paths
        self.rgi_file = None
        self.resfinder_file = None
        self.argannot_file = None
        
        # Load data
        self.rgi_df = None
        self.resfinder_df = None
        self.argannot_df = None
        
    def normalize_gene_name(self, name: str, database: str = "rgi") -> str:
        """Normalize gene names for comparison across databases
        common suffix/prefix removes so that gene family-level comparison becomes easier
        """
        name = str(name).upper().strip()
        
        # Remove common suffixes and prefixes
        name = re.sub(r'_\d+$', '', name)  # Remove _1, _2 suffixes
        name = re.sub(r'\([^)]*\)', '', name)  # Remove parentheses content
        name = re.sub(r'[-\d]+$', '', name)  # Remove trailing numbers
        name = re.sub(r'[A-Z]+\d+$', '', name)  # Remove trailing alphanumeric
        
        # Map common gene families to canonical names
        mappings = {
            'TEM': 'TEM',
            'CTXM': 'CTX-M',
            'SHV': 'SHV',
            'OXA': 'OXA',
            'KPC': 'KPC',
            'NDM': 'NDM',
            'VIM': 'VIM',
            'IMP': 'IMP',
            'APH': 'APH',
            'AAC': 'AAC',
            'AAD': 'AAD',
            'QNR': 'QNR',
            'TET': 'TET',
            'SUL': 'SUL',
            'DFR': 'DFR',
            'MCR': 'MCR',
            'FLO': 'FLO',
            'CMLA': 'CMLA',
            'MPH': 'MPH',
            'ERM': 'ERM'
        }
        # return mapped value if found
        for key, value in mappings.items():
            if key in name:
                return value
        
        # take first token or truncate long names
        name = name.split()[0] if ' ' in name else name
        
        # Keep only first 10 characters if still too long
        if len(name) > 15:
            name = name[:15]
        
        return name if name else 'UNKNOWN'
    
    def load_per_sample_rgi(self, file_list):
        """Load RGI results from per-sample files"""
        dfs = []
        for f in file_list:
            sample = Path(f).stem  # .txt থাকলে stem দিয়ে
            try:
                df = pd.read_csv(f, sep='\t')
                df['sample_id'] = sample
                
                # Apply same filtering as before
                if 'Cut_Off' in df.columns:
                    df = df[df['Cut_Off'].isin(['Perfect', 'Strict'])]
                
                # Normalize gene names
                if 'Best_Hit_ARO' in df.columns:
                    df['gene_norm'] = df['Best_Hit_ARO'].apply(
                        lambda x: self.normalize_gene_name(x, 'rgi')
                    )
                
                dfs.append(df)
            except Exception as e:
                logger.warning(f"Could not read RGI file {f}: {e}")
        
        if dfs:
            self.rgi_df = pd.concat(dfs, ignore_index=True)
            logger.info(f"  RGI: {len(self.rgi_df)} hits from {len(file_list)} files")
            logger.info(f"  RGI: {self.rgi_df['sample_id'].nunique()} samples loaded")
        else:
            logger.warning("No RGI data loaded")

    def load_per_sample_resfinder(self, file_list):
        """Load ResFinder results from per-sample files"""
        dfs = []
        for f in file_list:
            sample = Path(f).stem.replace('_results', '')  # sample_results.tsv -> sample
            try:
                df = pd.read_csv(f, sep='\t')
                df['sample_id'] = sample
                
                # Filter for acquired genes
                if 'type' in df.columns:
                    df = df[df['type'] == 'acquired_gene']
                
                # Normalize gene names
                if 'gene' in df.columns:
                    df['gene_norm'] = df['gene'].apply(
                        lambda x: self.normalize_gene_name(x, 'resfinder')
                    )
                
                dfs.append(df)
            except Exception as e:
                logger.warning(f"Could not read ResFinder file {f}: {e}")
        
        if dfs:
            self.resfinder_df = pd.concat(dfs, ignore_index=True)
            logger.info(f"  ResFinder: {len(self.resfinder_df)} hits from {len(file_list)} files")
            logger.info(f"  ResFinder: {self.resfinder_df['sample_id'].nunique()} samples loaded")
        else:
            logger.warning("No ResFinder data loaded")

    def load_per_sample_argannot(self, file_list):
        """Load ARG-ANNOT results from per-sample files"""
        dfs = []
        for f in file_list:
            sample = Path(f).stem.replace('_results', '')  # sample_results.tsv -> sample
            try:
                df = pd.read_csv(f, sep='\t')
                df['sample_id'] = sample
                
                # Extract clean gene names from ARG-ANNOT format
                if 'gene' in df.columns:
                    df['gene_clean'] = df['gene'].apply(
                        lambda x: x.split('__')[-2] if '__' in str(x) else str(x).split('_')[0]
                    )
                    df['gene_norm'] = df['gene_clean'].apply(
                        lambda x: self.normalize_gene_name(x, 'argannot')
                    )
                
                # Filter by identity threshold
                if 'identity' in df.columns:
                    df = df[df['identity'] >= self.min_identity]
                
                dfs.append(df)
            except Exception as e:
                logger.warning(f"Could not read ARG-ANNOT file {f}: {e}")
        
        if dfs:
            self.argannot_df = pd.concat(dfs, ignore_index=True)
            logger.info(f"  ARG-ANNOT: {len(self.argannot_df)} hits from {len(file_list)} files")
            logger.info(f"  ARG-ANNOT: {self.argannot_df['sample_id'].nunique()} samples loaded")
        else:
            logger.warning("No ARG-ANNOT data loaded")
    
    def load_data(self, rgi_files=None, resfinder_files=None, argannot_files=None):
        """Load all three database results from per-sample files"""
        logger.info("Loading database results from per-sample files...")
        
        # Load RGI (CARD)
        if rgi_files:
            self.load_per_sample_rgi(rgi_files)
        else:
            logger.warning("No RGI files provided")
        
        # Load ResFinder
        if resfinder_files:
            self.load_per_sample_resfinder(resfinder_files)
        else:
            logger.warning("No ResFinder files provided")
        
        # Load ARG-ANNOT
        if argannot_files:
            self.load_per_sample_argannot(argannot_files)
        else:
            logger.warning("No ARG-ANNOT files provided")
    
    def create_master_table(self, samples=None):
        """Create a master table combining all database results"""
        frames = []
        
        if self.rgi_df is not None:
            # RGI: Best_Identities থেকে identity ও coverage parse
            rgi = self.rgi_df[['sample_id', 'Best_Hit_ARO', 'Best_Identities', 'Cut_Off']].copy()
            
            # Parse Best_Identities (format: "95.5|88.2")
            rgi['identity'] = rgi['Best_Identities'].apply(
                lambda x: float(str(x).replace('%','').split('|')[0]) if pd.notna(x) else None
            )
            rgi['coverage'] = rgi['Best_Identities'].apply(
                lambda x: float(str(x).replace('%','').split('|')[1]) if pd.notna(x) and '|' in str(x) else None
            )
            
            rgi.rename(columns={'Best_Hit_ARO': 'gene', 'Cut_Off': 'cutoff'}, inplace=True)
            rgi['database'] = 'RGI_CARD'
            frames.append(rgi[['sample_id', 'gene', 'identity', 'coverage', 'cutoff', 'database']])
        
        if self.resfinder_df is not None:
            res = self.resfinder_df[['sample_id', 'gene', 'identity', 'coverage']].copy()
            res['database'] = 'ResFinder'
            res['cutoff'] = None
            frames.append(res)
        
        if self.argannot_df is not None:
            arg = self.argannot_df[['sample_id', 'gene_clean', 'identity', 'coverage']].copy()
            arg.columns = ['sample_id', 'gene', 'identity', 'coverage']
            arg['database'] = 'ARG-ANNOT'
            arg['cutoff'] = None
            frames.append(arg)
        
        master = pd.concat(frames, ignore_index=True)
        return master
    
    def get_gene_sets(self, sample_id: str = None) -> Dict[str, Set[str]]:
        """Get sets of normalized genes detected by each database"""
        gene_sets = {}
        
        if sample_id:
            if self.rgi_df is not None:
                genes = self.rgi_df[self.rgi_df['sample_id'] == sample_id]['gene_norm'].dropna().unique()
                gene_sets['RGI (CARD)'] = set(genes)
            
            if self.resfinder_df is not None:
                genes = self.resfinder_df[self.resfinder_df['sample_id'] == sample_id]['gene_norm'].dropna().unique()
                gene_sets['ResFinder'] = set(genes)
            
            if self.argannot_df is not None:
                genes = self.argannot_df[self.argannot_df['sample_id'] == sample_id]['gene_norm'].dropna().unique()
                gene_sets['ARG-ANNOT'] = set(genes)
        else:
            if self.rgi_df is not None:
                gene_sets['RGI (CARD)'] = set(self.rgi_df['gene_norm'].dropna().unique())
            
            if self.resfinder_df is not None:
                gene_sets['ResFinder'] = set(self.resfinder_df['gene_norm'].dropna().unique())
            
            if self.argannot_df is not None:
                gene_sets['ARG-ANNOT'] = set(self.argannot_df['gene_norm'].dropna().unique())
        
        return gene_sets
    
    def create_overlap_table(self, samples: List[str] = None) -> pd.DataFrame:
        """Create a table showing gene detection overlap across all samples"""
        
        # Get all samples from all databases
        all_samples = set()
        if self.rgi_df is not None:
            all_samples.update(self.rgi_df['sample_id'].unique())
        if self.resfinder_df is not None:
            all_samples.update(self.resfinder_df['sample_id'].unique())
        if self.argannot_df is not None:
            all_samples.update(self.argannot_df['sample_id'].unique())
        
        # If specific samples provided, use those
        if samples:
            all_samples = set(samples) & all_samples
        
        if not all_samples:
            logger.error("No samples found")
            return pd.DataFrame()
        
        logger.info(f"Processing {len(all_samples)} samples...")
        
        results = []
        for sample in sorted(all_samples):
            gene_sets = self.get_gene_sets(sample)
            
            # Get sets (empty if missing)
            rgi_set = gene_sets.get('RGI (CARD)', set())
            resfinder_set = gene_sets.get('ResFinder', set())
            argannot_set = gene_sets.get('ARG-ANNOT', set())
            
            # Calculate overlaps
            rgi_only = rgi_set - resfinder_set - argannot_set
            resfinder_only = resfinder_set - rgi_set - argannot_set
            argannot_only = argannot_set - rgi_set - resfinder_set
            
            rgi_resfinder = (rgi_set & resfinder_set) - argannot_set
            rgi_argannot = (rgi_set & argannot_set) - resfinder_set
            resfinder_argannot = (resfinder_set & argannot_set) - rgi_set
            
            all_three = rgi_set & resfinder_set & argannot_set
            
            results.append({
                'sample_id': sample,
                'RGI_only': len(rgi_only),
                'ResFinder_only': len(resfinder_only),
                'ARG-ANNOT_only': len(argannot_only),
                'RGI+ResFinder': len(rgi_resfinder),
                'RGI+ARG-ANNOT': len(rgi_argannot),
                'ResFinder+ARG-ANNOT': len(resfinder_argannot),
                'All_three': len(all_three),
                'Total_unique_genes': len(rgi_set | resfinder_set | argannot_set)
            })
        
        return pd.DataFrame(results)
    
    def plot_venn(self, save_path: Path = None):
        """Create Venn diagram for overall gene overlap"""
        gene_sets = self.get_gene_sets()
        
        if len(gene_sets) != 3:
            logger.error("Need exactly 3 databases for Venn diagram")
            return
        
        sets = list(gene_sets.values())
        labels = list(gene_sets.keys())
        
        # Calculate sizes
        only_rgi = len(sets[0] - sets[1] - sets[2])
        only_resfinder = len(sets[1] - sets[0] - sets[2])
        only_argannot = len(sets[2] - sets[0] - sets[1])
        
        rgi_resfinder = len((sets[0] & sets[1]) - sets[2])
        rgi_argannot = len((sets[0] & sets[2]) - sets[1])
        resfinder_argannot = len((sets[1] & sets[2]) - sets[0])
        
        all_three = len(sets[0] & sets[1] & sets[2])
        
        plt.figure(figsize=(10, 8))
        
        # Create Venn diagram
        venn = venn3(subsets=(only_rgi, only_resfinder, rgi_resfinder, only_argannot, 
                              rgi_argannot, resfinder_argannot, all_three),
                     set_labels=labels)
        
        # Color patches
        color_map = {
            '100': '#66c2a5',
            '010': '#fc8d62',
            '001': '#8da0cb',
            '110': '#e78ac3',
            '101': '#a6d854',
            '011': '#ffd92f',
            '111': '#e5c494'
        }
        
        for patch_id, color in color_map.items():
            patch = venn.get_patch_by_id(patch_id)
            if patch is not None:
                patch.set_color(color)
        
        plt.title(f'AMR Gene Detection Overlap (Normalized Names)\n(Total: {len(sets[0] | sets[1] | sets[2])} unique gene families)', 
                  fontsize=14, fontweight='bold')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Venn diagram saved to {save_path}")
        plt.show()
    
    def plot_heatmap(self, overlap_df: pd.DataFrame, save_path: Path = None):
        """Create heatmap of overlap patterns"""
        if overlap_df.empty:
            logger.warning("No data for heatmap")
            return
        
        heatmap_data = overlap_df.set_index('sample_id')
        heatmap_data = heatmap_data[['RGI_only', 'ResFinder_only', 'ARG-ANNOT_only', 
                                      'RGI+ResFinder', 'RGI+ARG-ANNOT', 'ResFinder+ARG-ANNOT', 
                                      'All_three']]
        
        plt.figure(figsize=(14, max(8, len(heatmap_data) * 0.3)))
        sns.heatmap(heatmap_data, annot=True, fmt='d', cmap='YlOrRd', 
                    cbar_kws={'label': 'Number of Genes'})
        plt.title('AMR Gene Detection Patterns by Sample', fontsize=14, fontweight='bold')
        plt.xlabel('Detection Category', fontsize=12)
        plt.ylabel('Sample ID', fontsize=12)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Heatmap saved to {save_path}")
        plt.show()
    
    def print_agreement_summary(self, overlap_df: pd.DataFrame):
        """Print summary of agreement levels"""
        print("\n" + "="*60)
        print("DATABASE AGREEMENT SUMMARY")
        print("="*60)
        
        total_samples = len(overlap_df)
        print(f"\nTotal samples analyzed: {total_samples}")
        
        # Count samples with any agreement
        samples_with_common = overlap_df[overlap_df['All_three'] > 0]
        samples_with_pair = overlap_df[(overlap_df['RGI+ResFinder'] > 0) | 
                                        (overlap_df['RGI+ARG-ANNOT'] > 0) | 
                                        (overlap_df['ResFinder+ARG-ANNOT'] > 0)]
        
        print(f"Samples with genes found by all 3 databases: {len(samples_with_common)}")
        print(f"Samples with genes found by at least 2 databases: {len(samples_with_pair)}")
        
        # Average detection per sample
        print(f"\nAverage detection per sample:")
        print(f"  RGI only:        {overlap_df['RGI_only'].mean():.1f}")
        print(f"  ResFinder only:  {overlap_df['ResFinder_only'].mean():.1f}")
        print(f"  ARG-ANNOT only:  {overlap_df['ARG-ANNOT_only'].mean():.1f}")
        print(f"  All three:       {overlap_df['All_three'].mean():.1f}")
        
        # Total
        total_unique = overlap_df['Total_unique_genes'].sum()
        total_all_three = overlap_df['All_three'].sum()
        
        if total_unique > 0:
            print(f"\nOverall concordance: {total_all_three}/{total_unique} gene detections ({100*total_all_three/total_unique:.1f}%)")
    
    def export_comparison_table(self, output_file: Path):
        """Export detailed comparison table"""
        gene_sets = self.get_gene_sets()
        
        all_genes = sorted(gene_sets.get('RGI (CARD)', set()) | 
                          gene_sets.get('ResFinder', set()) | 
                          gene_sets.get('ARG-ANNOT', set()))
        
        comparison = []
        for gene in all_genes:
            comparison.append({
                'gene_family': gene,
                'RGI_CARD': gene in gene_sets.get('RGI (CARD)', set()),
                'ResFinder': gene in gene_sets.get('ResFinder', set()),
                'ARG_ANNOT': gene in gene_sets.get('ARG-ANNOT', set()),
                'databases_count': sum([
                    gene in gene_sets.get('RGI (CARD)', set()),
                    gene in gene_sets.get('ResFinder', set()),
                    gene in gene_sets.get('ARG-ANNOT', set())
                ])
            })
        
        df = pd.DataFrame(comparison)
        df.to_csv(output_file, index=False)
        logger.info(f"Gene-level comparison saved to {output_file}")
        
        # Print summary
        print("\n" + "="*60)
        print("GENE-LEVEL AGREEMENT (Normalized Names)")
        print("="*60)
        
        all_three_genes = df[df['databases_count'] == 3]['gene_family'].tolist()
        if all_three_genes:
            print(f"\nGene families found by all 3 databases ({len(all_three_genes)}):")
            print(f"  {', '.join(all_three_genes[:20])}")
        else:
            print("\nNo gene families found by all 3 databases")
        
        two_way = df[df['databases_count'] == 2]['gene_family'].tolist()
        if two_way:
            print(f"\nGene families found by 2 databases ({len(two_way)}):")
            print(f"  {', '.join(two_way[:15])}")

def main():
    parser = argparse.ArgumentParser(description='Compare AMR gene detection across databases')
    parser.add_argument('--samples', nargs='+', help='Specific sample IDs to analyze')
    parser.add_argument('--rgi_files', nargs='+', help='List of per-sample RGI output files')
    parser.add_argument('--resfinder_files', nargs='+', help='List of per-sample ResFinder output files')
    parser.add_argument('--argannot_files', nargs='+', help='List of per-sample ARG-ANNOT output files')
    parser.add_argument('--output', help='Path to output master AMR table (TSV)')
    parser.add_argument('--output_dir', default='results/comparison', help='Directory for additional outputs')
    parser.add_argument('--min_identity', type=float, default=80, help='Minimum identity for filtering')
    parser.add_argument('--no_venn', action='store_true', help='Skip Venn diagram')
    parser.add_argument('--no_heatmap', action='store_true', help='Skip heatmap')
    args = parser.parse_args()
    
    # Initialize comparator
    comparator = DatabaseComparator(min_identity=args.min_identity)
    
    # Load data from per-sample files
    comparator.load_data(
        rgi_files=args.rgi_files,
        resfinder_files=args.resfinder_files,
        argannot_files=args.argannot_files
    )
    
    # Snakefile integration mode (only output master table)
    if args.output:
        logger.info("Creating master table for Snakemake integration...")
        master_df = comparator.create_master_table(samples=args.samples)
        master_df.to_csv(args.output, sep='\t', index=False)
        logger.info(f"Master table saved to {args.output}")
        return
    
    # Create output directory for standalone mode
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create overlap table
    logger.info("Creating overlap table for all samples...")
    overlap_df = comparator.create_overlap_table(samples=args.samples)
    
    if overlap_df.empty:
        logger.error("No samples found for comparison")
        return
    
    # Save overlap table
    overlap_file = output_dir / "database_overlap_table.csv"
    overlap_df.to_csv(overlap_file, index=False)
    logger.info(f"Overlap table saved to {overlap_file} ({len(overlap_df)} samples)")
    
    # Print summary
    comparator.print_agreement_summary(overlap_df)
    
    # Plot Venn diagram
    if not args.no_venn:
        venn_file = output_dir / "venn_diagram.png"
        comparator.plot_venn(save_path=venn_file)
    
    # Plot heatmap
    if not args.no_heatmap and len(overlap_df) > 0:
        heatmap_file = output_dir / "agreement_heatmap.png"
        comparator.plot_heatmap(overlap_df, save_path=heatmap_file)
    
    # Export gene-level comparison
    gene_comparison_file = output_dir / "gene_level_comparison.csv"
    comparator.export_comparison_table(gene_comparison_file)
    
    print(f"\n All outputs saved to: {output_dir}")
    print(f"\n Quick stats:")
    print(f"   Total samples: {len(overlap_df)}")
    print(f"   Samples with All_three > 0: {(overlap_df['All_three'] > 0).sum()}")
    print(f"   Total All_three genes: {overlap_df['All_three'].sum()}")

if __name__ == "__main__":
    main()