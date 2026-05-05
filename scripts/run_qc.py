#!/usr/bin/env python3
"""
run_qc.py - AMR Pipeline QC Stage
Processes FASTQ pairs, generates QC metrics, flags failing samples
"""

import subprocess
import json
import pandas as pd
import glob
import yaml
import sys
from pathlib import Path
import argparse

PROJECT_ROOT = Path(__file__).parent.parent

def load_config():
    """Load QC thresholds from config.yaml"""
    config_path = PROJECT_ROOT / "config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def run_fastp(r1, r2, sample_id, qc_config):
    """Run fastp, return (success, json_path, metrics)"""
    
    trimmed_dir = PROJECT_ROOT / "data" / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    
    json_out = trimmed_dir / f"{sample_id}.json"
    html_out = trimmed_dir / f"{sample_id}.html"
    out_r1 = trimmed_dir / f"{sample_id}_R1_clean.fastq.gz"
    out_r2 = trimmed_dir / f"{sample_id}_R2_clean.fastq.gz"
    
    cmd = [
        "fastp",
        "-i", str(r1), "-I", str(r2),
        "-o", str(out_r1), "-O", str(out_r2),
        "-j", str(json_out), "-h", str(html_out),
        "--detect_adapter_for_pe",
        "--length_required", str(qc_config.get('min_length', 50))
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True, json_out
    except subprocess.CalledProcessError as e:
        print(f"   ❌ fastp failed: {e.stderr[:100]}")
        return False, None

def parse_metrics(json_path):
    """Extract Q30%, adapter content%, duplication% from fastp JSON"""
    with open(json_path) as f:
        data = json.load(f)
    
    before = data['summary']['before_filtering']
    after = data['summary']['after_filtering']
    
    q30 = after.get('q30_rate', 0) * 100
    duplication = data.get('duplication', {}).get('rate', 0) * 100
    
    adapter_bases = data.get('adapter_cutting', {}).get('adapter_trimmed_bases', 0)
    total_bases = before.get('total_bases', 1)
    adapter_pct = (adapter_bases / total_bases) * 100
    
    return {
        'q30_rate': round(q30, 2),
        'adapter_content_pct': round(adapter_pct, 3),
        'duplication_rate': round(duplication, 2)
    }
def resolve_path(path_str, default):
    """Resolve a user-supplied path relative to PROJECT_ROOT if not absolute"""
    if path_str:
        p = Path(path_str)
        if not p.is_absolute():
            p = PROJECT_ROOT / p
        return p
    return PROJECT_ROOT / default
def main():
    # ---------- 1. argument parsing ----------
    parser = argparse.ArgumentParser(description="AMR QC Script")
    parser.add_argument('--summarize', action='store_true',
                        help='Only summarize existing JSON files')
    parser.add_argument('--input_dir', type=str,
                        help='Directory containing fastp JSON files')
    parser.add_argument('--output', type=str,
                        help='Output TSV file path')
    parser.add_argument('--multiqc', action='store_true',
                    help='Generate MultiQC report (default mode only)')
    args = parser.parse_args()

    # ---------- 2. config loading ----------
    config = load_config()
    qc_config = config.get('qc', {})
    q30_cutoff = qc_config.get('min_q30', 80)
    adapter_cutoff = qc_config.get('adapter_content_threshold', 10)
    
    # ============================
    #    SUMMARIZE MODE
    # ============================
    if args.summarize:
        print("="*60)
        print(" QC Summary Mode - Aggregating fastp JSONs")
        print("="*60)
        
        # path resolution both inpur and output
        input_dir = resolve_path(args.input_dir, "data/trimmed")
        output_path = resolve_path(args.output, "data/master_qc_report.tsv")
        
        # searching json file
        json_files = list(input_dir.glob("*.json"))
        if not json_files:
            print(f"❌ No JSON files found in {input_dir}")
            sys.exit(1)
        
        print(f" Found {len(json_files)} JSON files\n")
        
        results = []
        for jf in json_files:
            sample_id = jf.stem   
            print(f" Processing {jf.name}...")
            
            metrics = parse_metrics(jf)
            metrics['sample_id'] = sample_id
            metrics['qc_pass'] = (
                metrics['q30_rate'] >= q30_cutoff and 
                metrics['adapter_content_pct'] <= adapter_cutoff
            )
            results.append(metrics)
            
            status = " PASS" if metrics['qc_pass'] else "❌ FAIL"
            print(f"   {status} | Q30: {metrics['q30_rate']}% | Adapter: {metrics['adapter_content_pct']}%")
        
        # saving TSV
        df = pd.DataFrame(results)
        output_path.parent.mkdir(parents=True, exist_ok=True) 
        df.to_csv(output_path, sep='\t', index=False)
        
        # printing report
        print("\n" + "="*60)
        print(" MASTER QC REPORT")
        print("="*60)
        print(f"Total processed: {len(df)}")
        print(f" PASS: {df['qc_pass'].sum()}")
        print(f" FAIL: {(~df['qc_pass']).sum()}")
        
        failed_df = df[~df['qc_pass']]
        if not failed_df.empty:
            print("\nFailed samples:")
            for _, row in failed_df.iterrows():
                print(f"  ✗ {row['sample_id']}: Q30={row['q30_rate']}%, Adapter={row['adapter_content_pct']}%")
        
        print(f"\n Master report: {output_path}")
        print("\n QC Summary Complete!")
        sys.exit(0)  
    
    # ============================
    #    DEFAULT MODE (fastp)
    # ============================
    else:
        print("="*60)
        print(" AMR Pipeline - QC Stage (fastp)")
        print(f"   Q30 cutoff: {q30_cutoff}%")
        print(f"   Adapter cutoff: {adapter_cutoff}%")
        print("="*60)
        
        fastq_dir = PROJECT_ROOT / "data" / "fastq"
        r1_files = list(fastq_dir.glob("*_1.fastq.gz"))
        
        if not r1_files:
            print(f"❌ No FASTQ files found in {fastq_dir}")
            sys.exit(1)
        
        print(f"\n Found {len(r1_files)} samples\n")
        
        results = []
        for r1 in r1_files:
            sample_id = r1.name.replace("_1.fastq.gz", "")
            r2 = fastq_dir / f"{sample_id}_2.fastq.gz"
            
            if not r2.exists():
                print(f"  {sample_id}: missing R2, skipping")
                continue
            
            print(f" {sample_id}...")
            success, json_path = run_fastp(r1, r2, sample_id, qc_config)
            
            if success and json_path and json_path.exists():
                metrics = parse_metrics(json_path)
                metrics['sample_id'] = sample_id
                metrics['qc_pass'] = (
                    metrics['q30_rate'] >= q30_cutoff and 
                    metrics['adapter_content_pct'] <= adapter_cutoff
                )
                results.append(metrics)
                
                status = " PASS" if metrics['qc_pass'] else "❌ FAIL"
                print(f"   {status} | Q30: {metrics['q30_rate']}% | Adapter: {metrics['adapter_content_pct']}%")
            else:
                print(f"    Failed to process {sample_id}")
        
        if results:
            df = pd.DataFrame(results)
            output_path = resolve_path(args.output, "data/master_qc_report.tsv")
            df.to_csv(output_path, sep='\t', index=False)
            
            print("\n" + "="*60)
            print(" MASTER QC REPORT")
            print("="*60)
            print(f"Total processed: {len(df)}")
            print(f" PASS: {df['qc_pass'].sum()}")
            print(f" FAIL: {(~df['qc_pass']).sum()}")
            
            failed_df = df[~df['qc_pass']]
            if not failed_df.empty:
                print("\nFailed samples:")
                for _, row in failed_df.iterrows():
                    print(f"  ✗ {row['sample_id']}: Q30={row['q30_rate']}%, Adapter={row['adapter_content_pct']}%")
            
            print(f"\n Master report: {output_path}")
            
            # MultiQC report 
            if args.multiqc:
                print("\n Generating MultiQC report...")
                trimmed_dir = PROJECT_ROOT / "data" / "trimmed"
                subprocess.run(["multiqc", str(trimmed_dir), "-m", "fastp", "-f", "-o", str(trimmed_dir)], 
                               capture_output=True)
                print(f" MultiQC report: {trimmed_dir}/multiqc_report.html")
        else:
            print("\n No successful QC results to report!")
        
        print("\n QC Stage Complete!")

if __name__ == "__main__":
    main()