#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samplesheet = null
params.manifest = null
params.genome = null
params.reference_config = null
params.aligner = 'bwa'
params.umi_length = 8
params.mapq = 20
params.map_threads = 8
params.sort_threads = 1
params.sort_memory = '768M'
params.dedup_method = 'directional'
params.outdir = 'results'

process EXTRACT_UMI {
    tag "${sample}"
    label 'small'
    cpus 2
    publishDir "${params.outdir}/prepared", mode: 'copy'

    input:
    tuple val(sample), val(barcode), path(reads)

    output:
    tuple val(sample), path("${sample}.umi.fastq.gz"), emit: prepared

    script:
    """
    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern='(?P<umi_1>.{${params.umi_length}})${barcode}' \\
        --stdin '${reads}' \\
        --stdout '${sample}.umi.fastq.gz' \\
        --log '${sample}.extract.log'
    """
}

process MAP_READS {
    tag "${sample} (${params.genome}, ${params.aligner})"
    label 'mapping'
    cpus params.map_threads + params.sort_threads
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    val genome
    path reference_config

    output:
    tuple val(sample), path("${sample}.q${params.mapq}.bam"), emit: filtered_bam
    path "${sample}.mapping.json", emit: mapping_qc
    path "${sample}.all.bam", emit: all_bam
    path "${sample}.all.bam.bai", emit: all_bai
    path "${sample}.q${params.mapq}.bam.bai", emit: filtered_bai
    path "${sample}.*.log", optional: true, emit: logs

    script:
    """
    python ${projectDir}/scripts/nextflow_stage.py map \\
        --sample '${sample}' \\
        --read1 '${reads}' \\
        --genome '${genome}' \\
        --reference-config '${reference_config}' \\
        --aligner '${params.aligner}' \\
        --mapq '${params.mapq}' \\
        --threads '${params.map_threads}' \\
        --sort-threads '${params.sort_threads}' \\
        --sort-memory '${params.sort_memory}' \\
        --outdir .
    """
}

process DEDUPLICATE {
    tag "${sample}"
    label 'small'
    cpus 1
    publishDir "${params.outdir}/deduplication", mode: 'copy'

    input:
    tuple val(sample), path(bam)

    output:
    path "${sample}.families.tsv"
    path "${sample}.counts.bed"
    path "${sample}.molecules.bed"
    path "${sample}.sites.tsv"
    path "${sample}.dedup.json"

    script:
    """
    python ${projectDir}/scripts/nextflow_stage.py dedup \\
        --sample '${sample}' \\
        --bam '${bam}' \\
        --method '${params.dedup_method}' \\
        --umi-length '${params.umi_length}' \\
        --mapq '${params.mapq}' \\
        --outdir .
    """
}

workflow {
    if (!params.manifest || !params.genome || !params.reference_config) {
        error 'Provide --manifest, --genome and --reference_config (see docs/nextflow.md)'
    }

    manifest = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.sample as String, row.barcode as String, file(row.fastq as String))
        }

    prepared = EXTRACT_UMI(manifest)
    mapped = MAP_READS(prepared, params.genome, file(params.reference_config, checkIfExists: true))
    DEDUPLICATE(mapped.filtered_bam)
}
