#!/usr/bin/env python3
"""
fetch_all.py - AMR Pipeline Universal Downloader
Supports both NCBI (prefetch) and ENA (wget) sources
Auto-detects from manifest 'source' column
"""

import pandas as pd
import subprocess
import os
import time
import json
import requests
import warnings
from pathlib import Path
from typing import List, Optional
import argparse 
warnings.filterwarnings("ignore", category=FutureWarning)

MAX_RETRIES = 3
BASE_DELAY = 5

PROJECT_ROOT = Path(__file__).parent.parent
FASTQ_DIR = PROJECT_ROOT / "data" / "fastq"
MANIFEST_PATH = PROJECT_ROOT / "data" / "master_manifest.csv"

ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
def run_with_retry(cmd: List[str], max_retries: int = MAX_RETRIES, base_delay: int = BASE_DELAY):
    """Run shell command with retry logic"""
    for attempt in range(max_retries):
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            return result
        except subprocess.CalledProcessError as e:
            print(f"   Attempt {attempt + 1}/{max_retries} failed")
            if e.stderr:
                print(f"   Error: {e.stderr[:200]}")
            if attempt == max_retries - 1:
                raise
            delay = base_delay * (2 ** attempt)
            print(f"   Retrying in {delay} seconds...")
            time.sleep(delay)
def get_ena_ftp_link(srr_id: str) -> Optional[str]:
    params = {
        "accession": srr_id,
        "result": "read_run",
        "format": "json",
        "fields": "fastq_ftp"
    }
    try:
        response = requests.get(ENA_API_URL, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        if data and data[0].get("fastq_ftp"):
            ftp_paths = data[0]["fastq_ftp"].split(";")
            if ftp_paths:
                return f"ftp://{ftp_paths[0]}"
        print(f"     No FTP link found for {srr_id} in ENA")
        return None
    except Exception as e:
        print(f"     ENA API error: {e}")
        return None
def download_from_ena(srr_id: str):
    ftp_link = get_ena_ftp_link(srr_id)
    if not ftp_link:
        raise Exception(f"Could not get ENA FTP link for {srr_id}")
    filename = ftp_link.split("/")[-1]
    output_path = FASTQ_DIR / filename
    
    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"    Already downloaded: {filename}")
        return
    
    run_with_retry(["wget", "-c", ftp_link, "-P", str(FASTQ_DIR)])
    
    # For paired-end: get second file if exists
    if "_1.fastq.gz" in filename:
        second_link = ftp_link.replace("_1.fastq.gz", "_2.fastq.gz")
        second_file = FASTQ_DIR / second_link.split("/")[-1]
        if not (second_file.exists() and second_file.stat().st_size > 0):
            try:
                run_with_retry(["wget", "-c", second_link, "-P", str(FASTQ_DIR)])
            except:
                print(f"     No second file (probably single-end)")
    
    # Rename to standard {sample}_1.fastq.gz and {sample}_2.fastq.gz
    if "_1.fastq.gz" in filename:
        expected_1 = FASTQ_DIR / f"{srr_id}_1.fastq.gz"
        actual_1 = FASTQ_DIR / filename
        if actual_1.exists():
            actual_1.rename(expected_1)
            print(f"       Renamed: {filename} → {srr_id}_1.fastq.gz")
        
        # Handle second file for paired-end
        second_actual = FASTQ_DIR / filename.replace("_1.fastq.gz", "_2.fastq.gz")
        expected_2 = FASTQ_DIR / f"{srr_id}_2.fastq.gz"
        if second_actual.exists():
            second_actual.rename(expected_2)
            print(f"       Renamed: {second_actual.name} → {srr_id}_2.fastq.gz")
    elif ".fastq.gz" in filename and "_1" not in filename and "_2" not in filename:
        # Single-end file - rename to _1 format for consistency
        actual = FASTQ_DIR / filename
        expected = FASTQ_DIR / f"{srr_id}_1.fastq.gz"
        if actual.exists():
            actual.rename(expected)
            print(f"       Renamed single-end: {filename} → {srr_id}_1.fastq.gz")
def download_from_ncbi(srr_id: str):
    run_with_retry(["prefetch", srr_id])
    run_with_retry([
        "fasterq-dump",
        srr_id,
        "--split-files",
        "--outdir", str(FASTQ_DIR)
    ])
    for suffix in ["_1.fastq", "_2.fastq", ".fastq"]:
        fq_file = FASTQ_DIR / f"{srr_id}{suffix}"
        if fq_file.exists() and fq_file.stat().st_size > 0:
            try:
                subprocess.run(["pigz", "-f", str(fq_file)], check=True)
                print(f"       Compressed (pigz): {fq_file.name}.gz")
            except (FileNotFoundError, subprocess.CalledProcessError):
                try:
                    subprocess.run(["gzip", "-f", str(fq_file)], check=True)
                    print(f"       Compressed (gzip): {fq_file.name}.gz")
                except Exception as e:
                    print(f"        Compression failed: {e}")
def download_sample(srr_id: str, source: str):
    if source.upper() == "ENA":
        print(f"    Downloading from ENA...")
        download_from_ena(srr_id)
    elif source.upper() == "NCBI":
        print(f"    Downloading from NCBI...")
        download_from_ncbi(srr_id)
    else:
        raise Exception(f"Unknown source: {source}. Use 'NCBI' or 'ENA'")
def fetch_samples():
    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    if not MANIFEST_PATH.exists():
        print(f" Manifest not found: {MANIFEST_PATH}")
        return
    df = pd.read_csv(MANIFEST_PATH)
    if 'sample_id' not in df.columns:
        print(" Manifest missing 'sample_id' column")
        return
    if 'source' not in df.columns:
        print("  No 'source' column found. Assuming NCBI for all.")
        df['source'] = 'NCBI'
    if 'download_status' not in df.columns:
        df['download_status'] = 'Pending'
    if 'error_log' not in df.columns:
        df['error_log'] = pd.Series(dtype='str')
    df['error_log'] = df['error_log'].fillna('').astype(str)
    
    pending = df[df['download_status'] == 'Pending']
    if pending.empty:
        print(" All samples already downloaded!")
        return
    print(f"\n Downloading {len(pending)} samples...")
    for idx, (row_idx, row) in enumerate(pending.iterrows(), 1):
        srr_id = row['sample_id']
        source = row.get('source', 'NCBI')
        print(f"[{idx}/{len(pending)}] {srr_id} ({source}) ...")
        try:
            download_sample(srr_id, source)
            df.at[row_idx, 'download_status'] = 'Downloaded'
            df.at[row_idx, 'error_log'] = ''
            print(f"    Done!\n")
        except Exception as e:
            df.at[row_idx, 'download_status'] = 'Failed'
            df.at[row_idx, 'error_log'] = str(e)[:500]
            print(f"    Failed: {e}\n")
        tmp_path = MANIFEST_PATH.with_suffix(".tmp")
        df.to_csv(tmp_path, index=False)
        os.replace(tmp_path, MANIFEST_PATH)
    # summary
    df_final = pd.read_csv(MANIFEST_PATH)
    print("="*50)
    print(f" Downloaded: {len(df_final[df_final['download_status'] == 'Downloaded'])}")
    print(f" Failed: {len(df_final[df_final['download_status'] == 'Failed'])}")
    print("="*50)
def main():
    global MANIFEST_PATH
    parser = argparse.ArgumentParser(description="Download FASTQ for one or all samples")
    parser.add_argument("--sample", help="Specific sample ID to download")
    parser.add_argument("--manifest", help="Path to manifest CSV file")
    args = parser.parse_args()
      
    if args.manifest:
        MANIFEST_PATH = Path(args.manifest)

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)

    if args.sample:
        df = pd.read_csv(MANIFEST_PATH)
        row = df[df['sample_id'] == args.sample]
        if row.empty:
            print(f" Sample {args.sample} not found in manifest")
            return
        source = row.iloc[0].get('source', 'NCBI')
        print(f"Downloading single sample: {args.sample} ({source})")
        try:
            download_sample(args.sample, source)
            print(" Done")
        except Exception as e:
            print(f" Failed: {e}")
    else:
        fetch_samples()

if __name__ == "__main__":
    print("=" * 50)
    print(" AMR Pipeline Universal Downloader")
    print("   Supports: NCBI (prefetch) + ENA (wget)")
    print("=" * 50)
    print(f"Project: {PROJECT_ROOT}")
    print(f"Output: {FASTQ_DIR}")
    print("=" * 50)
    main()