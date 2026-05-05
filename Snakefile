# Snakefile for AMR Pipeline
# Usage: snakemake --use-conda --cores 4


configfile: "config.yaml"

import pandas as pd

manifest_df = pd.read_csv(config["paths"]["manifest"])
SAMPLES = manifest_df['sample_id'].tolist()

rule all:
    input:
        "results/master_amr_table.tsv",
        "data/master_qc_report.tsv"

# ── 1. FETCH ────────────────────────────────────────────────────────
rule fetch:
    output:
        r1 = "data/fastq/{sample}_1.fastq.gz",
        r2 = "data/fastq/{sample}_2.fastq.gz"
    log:
        "logs/fetch/{sample}.log"
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        "python scripts/fetch_all.py --sample {wildcards.sample} --manifest {config[paths][manifest]} > {log} 2>&1"

# ── 2. QC ────────────────────────────────────────────────────────────
rule qc:
    input:
        r1 = "data/fastq/{sample}_1.fastq.gz",
        r2 = "data/fastq/{sample}_2.fastq.gz"
    output:
        trimmed_r1 = "data/trimmed/{sample}_R1_clean.fastq.gz",
        trimmed_r2 = "data/trimmed/{sample}_R2_clean.fastq.gz",
        json       = "data/trimmed/{sample}.json",
        html       = "data/trimmed/{sample}.html"
    log:
        "logs/qc/{sample}.log"
    params:
        min_length  = config['qc']['min_length'],
        min_quality = config['qc']['min_mean_quality']
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.trimmed_r1} -O {output.trimmed_r2} \
              -j {output.json} -h {output.html} \
              --detect_adapter_for_pe \
              --length_required {params.min_length} \
              --average_qual {params.min_quality} \
              > {log} 2>&1
        """

# ── 2b. QC SUMMARY ───────────────────────────────────────────────────
rule qc_summary:
    input:
        jsons = expand("data/trimmed/{sample}.json", sample=SAMPLES)
    output:
        "data/master_qc_report.tsv"
    log:
        "logs/qc_summary.log"
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        "python scripts/run_qc.py --summarize --input_dir data/trimmed/ --output {output} > {log} 2>&1"

# ── 3. ASSEMBLY ──────────────────────────────────────────────────────
rule assemble:
    input:
        r1 = "data/trimmed/{sample}_R1_clean.fastq.gz",
        r2 = "data/trimmed/{sample}_R2_clean.fastq.gz"
    output:
        "data/assemblies/{sample}/contigs.fasta"
    threads:
        config['assembly']['threads']
    resources:
        mem_mb = config['resources']['high']['mem_mb']
    log:
        "logs/assemble/{sample}.log"
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        """
        spades.py --isolate \
                  -1 {input.r1} -2 {input.r2} \
                  -o data/assemblies/{wildcards.sample} \
                  --threads {threads} \
                  --memory {config[assembly][memory_gb]} \
                  > {log} 2>&1
        """

# ── 4. RGI ───────────────────────────────────────────────────────────
rule rgi_detect:
    input:
        contigs = "data/assemblies/{sample}/contigs.fasta"
    output:
        txt  = "results/rgi/{sample}.txt",
        json = "results/rgi/{sample}.json"
    log:
        "logs/rgi/{sample}.log"
    conda:
        "envs/rgi_env.yaml"
    shell:
        """
        rgi main -i {input.contigs} \
                 -o results/rgi/{wildcards.sample} \
                 -t contig -a DIAMOND --clean -n 4 \
                 > {log} 2>&1
        """

# ── 5. RESFINDER ─────────────────────────────────────────────────────
rule resfinder_detect:
    input:
        contigs = "data/assemblies/{sample}/contigs.fasta"
    output:
        "results/resfinder/{sample}_results.tsv"
    log:
        "logs/resfinder/{sample}.log"
    conda:
        "envs/resfinder_env.yaml"
    shell:
        """
        python scripts/run_resfinder.py \
               --sample {wildcards.sample} \
               --output {output} \
               > {log} 2>&1
        """

# ── 6. ARG-ANNOT ─────────────────────────────────────────────────────
rule argannot_blast:
    input:
        contigs = "data/assemblies/{sample}/contigs.fasta"
    output:
        "results/argannot/{sample}_results.tsv"
    log:
        "logs/argannot/{sample}.log"
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        """
        python scripts/run_argannot_blast.py \
               --sample {wildcards.sample} \
               --output {output} \
               > {log} 2>&1
        """

# ── 7. INTEGRATE ─────────────────────────────────────────────────────
rule integrate_results:
    input:
        rgi       = expand("results/rgi/{sample}.txt", sample=SAMPLES),
        resfinder = expand("results/resfinder/{sample}_results.tsv", sample=SAMPLES),
        argannot  = expand("results/argannot/{sample}_results.tsv", sample=SAMPLES),
        manifest  = config["paths"]["manifest"]
    output:
        "results/master_amr_table.tsv"
    log:
        "logs/integrate.log"
    conda:
        "envs/amr_pipeline.yaml"
    shell:
        "python scripts/compare_database.py --rgi_files {input.rgi} --resfinder_files {input.resfinder} --argannot_files {input.argannot} --output {output} > {log} 2>&1"