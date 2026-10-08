#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samplesheet = null
params.manifest = null
params.genome = null
params.reference_config = null
params.aligner = 'bwa'
params.umi_length = 8
params.barcode_orientation = 'forward'
params.minimum_insert_length = 20
params.technical_prefix = ''
params.mapq = 20
params.map_threads = 8
params.sort_threads = 1
params.sort_memory = '768M'
params.dedup_method = 'directional'
params.outdir = 'results'
params.coverage_bin_size = 50

process CHECK_ENVIRONMENT {
    label 'small'
    cpus 1
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path reference_config
    val genome
    val aligner

    output:
    path 'environment.json', emit: ready

    script:
    """
    python ${projectDir}/scripts/post_mapping.py check \\
        --reference-config '${reference_config}' --genome '${genome}' --aligner '${aligner}'
    """
}

process EXTRACT_UMI {
    tag "${meta.sample}"
    label 'small'
    cpus 2
    publishDir "${params.outdir}/prepared", mode: 'copy'
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/qc" }, mode: 'copy', pattern: '*.preparation*.json'

    input:
    tuple val(meta), val(barcode), path(reads)
    path environment_check

    output:
    tuple val(meta), path("${meta.sample}.umi.fastq.gz"), emit: prepared
    path "${meta.sample}.preparation.json", emit: preparation_qc
    path "${meta.sample}.preparation_mqc.json", emit: preparation_multiqc

    script:
    """
    python ${projectDir}/scripts/prepare_bliss_reads.py \\
        --reads '${reads}' --sample '${meta.sample}' --barcode '${barcode}' \\
        --barcode-orientation '${params.barcode_orientation}' --umi-length '${params.umi_length}' \\
        --minimum-insert-length '${params.minimum_insert_length}' \\
        --technical-prefix '${params.technical_prefix}' --outdir .
    """
}

process MAP_READS {
    tag "${meta.sample} (${params.genome}, ${params.aligner})"
    label 'mapping'
    // samtools sort -@ N adds N worker threads plus its main thread.
    cpus { (params.map_threads as Integer) + (params.sort_threads as Integer) + 1 }
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/mapping" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)
    val genome
    path reference_config

    output:
    tuple val(meta), path("${meta.sample}.q${params.mapq}.bam"), path("${meta.sample}.q${params.mapq}.bam.bai"), emit: filtered_bam
    tuple val(meta), path("${meta.sample}.mapping.json"), emit: mapping_qc
    tuple val(meta), path("${meta.sample}.all.bam"), path("${meta.sample}.all.bam.bai"), emit: all_bam
    path "${meta.sample}.*.log", optional: true, emit: logs

    script:
    """
    python ${projectDir}/scripts/nextflow_stage.py map \\
        --sample '${meta.sample}' \\
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
    tag "${meta.sample}"
    label 'small'
    cpus 1
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/deduplication" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.sample}.families.tsv"
    tuple val(meta), path("${meta.sample}.counts.bed"), path("${meta.sample}.dedup.json"), emit: end_counts
    path "${meta.sample}.molecules.bed"
    path "${meta.sample}.sites.tsv"

    script:
    """
    python ${projectDir}/scripts/nextflow_stage.py dedup \\
        --sample '${meta.sample}' \\
        --bam '${bam}' \\
        --method '${params.dedup_method}' \\
        --umi-length '${params.umi_length}' \\
        --mapq '${params.mapq}' \\
        --outdir .
    """
}

process SAMPLE_QC {
    tag "${meta.sample}"
    label 'reporting'
    cpus 2
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/qc" }, mode: 'copy', pattern: '*.{txt,json,tsv}'
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/bigwig" }, mode: 'copy', pattern: '*.bw'
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/bigwig" }, mode: 'copy', pattern: '*.ends.counts.bed'

    input:
    tuple val(meta), path(filtered), path(filtered_bai), path(all_bam), path(all_bai), path(mapping_qc), path(counts), path(dedup_qc)

    output:
    tuple val(meta), path("${meta.sample}.summary.json"), path("${meta.sample}.chrom.sizes.json"), path(counts), path(filtered), path(filtered_bai), emit: condition_input
    path "sample__${meta.sample}.*.txt", emit: samtools_qc
    path "${meta.sample}.qc.tsv"
    path "${meta.sample}.*.bw"
    path "${meta.sample}.*.provenance.json"
    path "${meta.sample}.ends.counts.bed"

    script:
    """
    python ${projectDir}/scripts/post_mapping.py sample \\
        --sample '${meta.sample}' --group '${meta.group}' --replicate '${meta.replicate}' \\
        --mapping '${mapping_qc}' --dedup '${dedup_qc}' \\
        --all-bam '${all_bam}' --filtered-bam '${filtered}' --counts '${counts}' \\
        --threads ${task.cpus} --bin-size ${params.coverage_bin_size}
    """
}

process CONDITION_QC {
    tag "${group}"
    label 'mapping'
    cpus 2
    publishDir { "${params.outdir}/MergedReplicate/${group}/qc" }, mode: 'copy', pattern: '*.{json,tsv}'
    publishDir { "${params.outdir}/MergedReplicate/${group}/bigwig" }, mode: 'copy', pattern: '*.{bw,bed}'

    publishDir { "${params.outdir}/MergedReplicate/${group}/mapping" }, mode: 'copy', pattern: '*.{bam,bai}'
    publishDir { "${params.outdir}/MergedReplicate/${group}/qc" }, mode: 'copy', pattern: '*.txt'

    input:
    tuple val(group), path(summaries, arity: '1..*'), path(counts, arity: '1..*'), path(headers, arity: '1..*'), path(bams, arity: '1..*'), path(bais, arity: '1..*')

    output:
    path "${group}.condition.summary.json", emit: summary
    path "${group}.*.tsv"
    path "${group}.*.bw"
    path "${group}.ends.counts.bed"
    path "${group}.*.provenance.json"
    tuple path("${group}.condition.filtered.bam"), path("${group}.condition.filtered.bam.bai"), emit: merged_bam
    path "condition__${group}.filtered.*.txt", emit: samtools_qc

    script:
    def summaryArgs = summaries.collect { "'${it}'" }.join(' ')
    def countArgs = counts.collect { "'${it}'" }.join(' ')
    def headerArgs = headers.collect { "'${it}'" }.join(' ')
    def bamArgs = bams.collect { "'${it}'" }.join(' ')
    """
    python ${projectDir}/scripts/post_mapping.py group --group '${group}' \\
        --summaries ${summaryArgs} \\
        --counts ${countArgs} \\
        --headers ${headerArgs} \\
        --bams ${bamArgs} --threads ${task.cpus} --bin-size ${params.coverage_bin_size}
    """
}

process MULTIQC {
    label 'reporting'
    cpus 1
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path samples, stageAs: 'sample_summaries/*', arity: '1..*'
    path conditions, stageAs: 'condition_summaries/*', arity: '1..*'
    path samtools_reports, arity: '1..*'
    path preparation_reports, arity: '1..*'

    output:
    path 'multiqc_report.html'
    path 'multiqc_data'
    path 'helixbusters_*'
    path 'reporting_versions.txt'

    script:
    def sampleArgs = samples.collect { "'${it}'" }.join(' ')
    def conditionArgs = conditions.collect { "'${it}'" }.join(' ')
    """
    python ${projectDir}/scripts/post_mapping.py report \\
        --samples ${sampleArgs} --conditions ${conditionArgs} --genome '${params.genome}'
    multiqc . --filename multiqc_report.html --outdir . --data-dir --cl-config 'data_dir_name: multiqc_data' --config helixbusters_multiqc_config.json
    multiqc --version > reporting_versions.txt
    samtools --version >> reporting_versions.txt
    """
}

workflow {
    if (!params.manifest || !params.genome || !params.reference_config) {
        error 'Provide --manifest, --genome and --reference_config (see docs/nextflow.md)'
    }
    ['map_threads', 'umi_length', 'coverage_bin_size', 'minimum_insert_length'].each { key ->
        if (!(params[key].toString() ==~ /[1-9][0-9]*/) || params[key].toLong() > Integer.MAX_VALUE) {
            error "--${key} must be a positive integer"
        }
    }
    if (!(params.sort_threads.toString() ==~ /0|[1-9][0-9]*/) || params.sort_threads.toLong() > Integer.MAX_VALUE) {
        error '--sort_threads must be a nonnegative integer'
    }
    if (!(params.mapq.toString() ==~ /0|[1-9][0-9]*/) || params.mapq.toLong() > 254) {
        error '--mapq must be an integer between 0 and 254'
    }
    if (!(params.aligner in ['bwa', 'bowtie2'])) { error '--aligner must be bwa or bowtie2' }
    if (!(params.dedup_method in ['exact', 'directional'])) { error '--dedup_method must be exact or directional' }
    if (!(params.barcode_orientation in ['forward', 'reverse_complement'])) {
        error '--barcode_orientation must be forward or reverse_complement'
    }
    if (params.technical_prefix && !(params.technical_prefix ==~ /[ACGT]{12,}/)) {
        error '--technical_prefix must contain at least 12 A/C/G/T bases'
    }

    // Validate the complete manifest before emitting any sample for processing.
    manifest = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .collect()
        .flatMap { rows ->
            if (!rows) { error 'Manifest contains no samples' }
            def seenSamples = [] as Set
            def seenReplicates = [] as Set
            rows.collect { row ->
                ['sample', 'group', 'replicate'].each { key ->
                    if (!(row[key] ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                        error "Invalid ${key} in manifest: ${row[key]}; use letters, digits, _, . or -"
                    }
                }
                if (!seenSamples.add(row.sample)) { error "Duplicate sample: ${row.sample}" }
                if (!seenReplicates.add([row.group, row.replicate])) {
                    error "Multiple libraries for group/replicate ${row.group}/${row.replicate}; expected one biological library"
                }
                if (!(row.barcode ==~ /[ACGT]+/)) { error "Invalid barcode for ${row.sample}" }
                if (!row.fastq) { error "Missing FASTQ for ${row.sample}" }
                tuple([sample: row.sample, group: row.group, replicate: row.replicate],
                      row.barcode, file(row.fastq, checkIfExists: true))
            }
        }

    reference = file(params.reference_config, checkIfExists: true)
    environment = CHECK_ENVIRONMENT(reference, params.genome, params.aligner)
    prepared = EXTRACT_UMI(manifest, environment.ready)
    mapped = MAP_READS(prepared.prepared, params.genome, reference)
    dedup = DEDUPLICATE(mapped.filtered_bam)
    sample_inputs = mapped.filtered_bam
        .join(mapped.all_bam)
        .join(mapped.mapping_qc)
        .join(dedup.end_counts)
    qc = SAMPLE_QC(sample_inputs)
    grouped = qc.condition_input
        .map { meta, summary, header, counts, bam, bai -> tuple(meta.group, meta.sample, summary, counts, header, bam, bai) }
        .groupTuple()
        .map { group, samples, summaries, counts, headers, bams, bais ->
            def order = (0..<samples.size()).toList().sort { i, j -> samples[i] <=> samples[j] }
            tuple(group, order.collect { summaries[it] }, order.collect { counts[it] }, order.collect { headers[it] },
                  order.collect { bams[it] }, order.collect { bais[it] })
        }
    conditions = CONDITION_QC(grouped)
    MULTIQC(qc.condition_input.map { meta, summary, header, counts, bam, bai -> summary }.collect(),
            conditions.summary.collect(), qc.samtools_qc.mix(conditions.samtools_qc).flatten().collect(),
            prepared.preparation_multiqc.collect())
}
