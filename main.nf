#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samplesheet = null
params.design = 'unspecified'
params.manifest = null
params.genome = null
params.reference_config = null
params.aligner = 'bwa'
params.umi_length = 8
params.barcode_orientation = 'forward'
params.minimum_insert_length = 20
params.technical_prefix = ''
params.technical_prefix_max_errors = 0
params.mapq = 20
params.map_threads = 8
params.sort_threads = 1
params.sort_memory = '768M'
params.dedup_method = 'directional'
params.outdir = 'results'
params.coverage_bin_size = 50
params.run_windows = false
params.run_peak_calling = false
params.window_sizes = '1000,5000,10000'
params.min_reps_consensus = 2
params.peak_width = 100
params.peak_qvalue = 0.01
params.peak_nolambda = false
params.effective_genome_size = 'auto'

process CHECK_ENVIRONMENT {
    label 'small'
    cpus 1
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path reference_config
    val genome
    val aligner
    path manifest_file

    output:
    path 'environment.json', emit: ready
    path 'design.*', emit: design_metadata

    script:
    """
    python ${projectDir}/scripts/validate_experimental_design.py \\
        --manifest '${manifest_file}' --design '${params.design}'
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
        --technical-prefix '${params.technical_prefix}' \\
        --technical-prefix-max-errors '${params.technical_prefix_max_errors}' --outdir .
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
    tuple val(meta), path("${meta.sample}.molecules.bed"), emit: molecules
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

process MACS3_CALLPEAK {
    tag "${meta.sample}"
    label 'reporting'
    cpus 1
    conda "${projectDir}/environment.yml"
    publishDir { "${params.outdir}/SingleReplicate/${meta.sample}/peaks" }, mode: 'copy'

    input:
    tuple val(meta), path(counts), path(header), path(molecules)

    output:
    tuple val(meta), path("${meta.sample}_peaks.narrowPeak"), path("${meta.sample}.provenance.json"), emit: peaks
    tuple val(meta), path("${meta.sample}_peaks.xls"), emit: xls
    tuple val(meta), path("${meta.sample}_summits.bed"), emit: bed
    path "${meta.sample}.versions.json", emit: versions
    path 'macs3.log', emit: log

    script:
    def background = params.peak_nolambda.toString() == 'true' ? '--nolambda' : ''
    def effectiveSize = params.effective_genome_size.toString() == 'auto' ?
        ([hg38: 2913022398L, mm10: 2652783500L][params.genome] ?: 1) : params.effective_genome_size
    """
    python ${projectDir}/scripts/call_bliss_peaks.py \\
        --sample '${meta.sample}' --counts '${counts}' --header '${header}' --molecules '${molecules}' \\
        --peak-width '${params.peak_width}' --peak-qvalue '${params.peak_qvalue}' \\
        --effective-genome-size '${effectiveSize}' ${background}
    """
}

process ANALYZE_REGIONS {
    label 'reporting'
    cpus 1
    publishDir "${params.outdir}", mode: 'copy', pattern: '{SingleReplicate,MergedReplicate}/**'
    publishDir "${params.outdir}/Analysis", mode: 'copy', pattern: '*.{tsv,bed,json}'

    input:
    tuple val(samples), path(counts, arity: '1..*'), path(headers, arity: '1..*'), path(molecules, arity: '1..*'), path(peak_files), path(peak_provenance)
    path design_files, arity: '1..*'

    output:
    path '*.{tsv,bed,json}', emit: tables
    path 'analysis_mqc.json', emit: multiqc
    path 'SingleReplicate/*/peaks/*', optional: true
    path 'MergedReplicate/*/peaks/*', optional: true

    script:
    def sampleArgs = samples.collect { "'${it}'" }.join(' ')
    def countArgs = counts.collect { "'${it}'" }.join(' ')
    def headerArgs = headers.collect { "'${it}'" }.join(' ')
    def moleculeArgs = molecules.collect { "'${it}'" }.join(' ')
    def externalPeaks = peak_files ? '--peak-files ' + peak_files.collect { "'${it}'" }.join(' ') +
        ' --peak-provenance ' + peak_provenance.collect { "'${it}'" }.join(' ') : ''
    def windows = params.run_windows.toString() == 'true' ? params.window_sizes : ''
    def peaks = params.run_peak_calling.toString() == 'true' ? '--peaks' : ''
    def background = params.peak_nolambda.toString() == 'true' ? '--nolambda' : ''
    def effectiveSize = params.effective_genome_size.toString() == 'auto' ?
        ([hg38: 2913022398L, mm10: 2652783500L][params.genome] ?: 1) : params.effective_genome_size
    """
    python ${projectDir}/scripts/analyze_regions.py \\
        --samples ${sampleArgs} --counts ${countArgs} --headers ${headerArgs} --molecules ${moleculeArgs} \\
        --design-file design.summary.json --windows '${windows}' ${peaks} ${background} ${externalPeaks} \\
        --min-reps-consensus '${params.min_reps_consensus}' --peak-width '${params.peak_width}' \\
        --peak-qvalue '${params.peak_qvalue}' --effective-genome-size '${effectiveSize}'
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
    path design_reports, arity: '1..*'
    path analysis_reports

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
    if (params.effective_genome_size.toString() == 'auto') {
        if (params.run_peak_calling.toString() == 'true' && !(params.genome in ['hg38', 'mm10'])) {
            error 'Provide --effective_genome_size for peak calling with this genome assembly'
        }
    } else if (!(params.effective_genome_size.toString() ==~ /[1-9][0-9]*/)) {
        error '--effective_genome_size must be auto or a positive integer'
    }
    ['run_windows', 'run_peak_calling', 'peak_nolambda'].each { key ->
        if (!(params[key].toString() in ['true', 'false'])) { error "--${key} must be true or false" }
    }
    if (!(params.window_sizes.toString() ==~ /[1-9][0-9]*(,[1-9][0-9]*)*/)) {
        error '--window_sizes must be comma-separated positive integers'
    }
    ['min_reps_consensus', 'peak_width'].each { key ->
        if (!(params[key].toString() ==~ /[1-9][0-9]*/) || params[key].toLong() > Integer.MAX_VALUE) {
            error "--${key} must be a positive integer"
        }
    }
    if ((params.peak_width as Integer) % 2 != 0) { error '--peak_width must be even' }
    if (!(params.peak_qvalue.toString() ==~ /[0-9.eE+-]+/) || (params.peak_qvalue as Double) <= 0 || (params.peak_qvalue as Double) >= 1) {
        error '--peak_qvalue must be between 0 and 1'
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
    if (!(params.design in ['paired', 'unpaired', 'unspecified'])) {
        error '--design must be paired, unpaired or unspecified'
    }
    if (!(params.dedup_method in ['exact', 'directional'])) { error '--dedup_method must be exact or directional' }
    if (!(params.barcode_orientation in ['forward', 'reverse_complement'])) {
        error '--barcode_orientation must be forward or reverse_complement'
    }
    if (params.technical_prefix && !(params.technical_prefix ==~ /[ACGT]{12,}/)) {
        error '--technical_prefix must contain at least 12 A/C/G/T bases'
    }
    if (!(params.technical_prefix_max_errors.toString() in ['0', '1', '2']) ||
        (params.technical_prefix_max_errors.toString() != '0' && !params.technical_prefix)) {
        error '--technical_prefix_max_errors must be 0, 1 or 2; a technical_prefix is required for tolerant matching'
    }

    // Validate the complete manifest before emitting any sample for processing.
    manifest = Channel
        .fromPath(params.manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .collect()
        .flatMap { rows ->
            if (!rows) { error 'Manifest contains no samples' }
            if (params.run_peak_calling.toString() == 'true' &&
                rows.groupBy { it.group }.any { group, members -> members.size() < (params.min_reps_consensus as Integer) }) {
                error '--min_reps_consensus exceeds the biological replicate count of a condition'
            }
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
                tuple([sample: row.sample, group: row.group, replicate: row.replicate, donor: row.donor ?: ''],
                      row.barcode, file(row.fastq, checkIfExists: true))
            }
        }

    reference = file(params.reference_config, checkIfExists: true)
    environment = CHECK_ENVIRONMENT(reference, params.genome, params.aligner,
                                    file(params.manifest, checkIfExists: true))
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
    analysis_reports = Channel.empty()
    if (params.run_windows.toString() == 'true' || params.run_peak_calling.toString() == 'true') {
        region_inputs = qc.condition_input
            .map { meta, summary, header, counts, bam, bai -> tuple(meta, counts, header) }
            .join(dedup.molecules)
        if (params.run_peak_calling.toString() == 'true') {
            calls = MACS3_CALLPEAK(region_inputs)
            region_inputs = region_inputs.join(calls.peaks)
        }
        region_inputs = region_inputs
            .collect(flat: false)
            .map { rows ->
                def ordered = rows.sort { a, b -> a[0].sample <=> b[0].sample }
                tuple(ordered.collect { it[0].sample }, ordered.collect { it[1] },
                      ordered.collect { it[2] }, ordered.collect { it[3] },
                      ordered.collect { it.size() > 4 ? it[4] : null }.findAll { it != null },
                      ordered.collect { it.size() > 5 ? it[5] : null }.findAll { it != null })
            }
        analysis = ANALYZE_REGIONS(region_inputs, environment.design_metadata.flatten().collect())
        analysis_reports = analysis.multiqc
    }
    MULTIQC(qc.condition_input.map { meta, summary, header, counts, bam, bai -> summary }.collect(),
            conditions.summary.collect(), qc.samtools_qc.mix(conditions.samtools_qc).flatten().collect(),
            prepared.preparation_multiqc.collect(), environment.design_metadata.flatten().collect(), analysis_reports.collect().ifEmpty([]))
}
