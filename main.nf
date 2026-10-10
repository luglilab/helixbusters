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
params.exact20_rescue = false
params.mapq = 20
params.map_threads = 8
params.sort_threads = 1
params.sort_memory = '768M'
params.dedup_method = 'directional'
params.outdir = 'results'
params.coverage_bin_size = 50
params.run_windows = false
params.run_pca = true
params.gtf = null
params.gtf_genome = null
params.promoter_upstream = 2000
params.promoter_downstream = 500
params.gene_min_reps = 2
params.gene_min_molecules = 2
params.window_sizes = '1000,5000,10000'
params.run_window_upset = true
params.window_min_molecules = 5
params.window_min_reps = 2
params.run_peak_calling = false // Legacy false is accepted; true is rejected.
params.run_differential = false
params.differential_mode = 'contrast'
params.timepoints = null
params.contrast = null
params.differential_min_count = 5
params.differential_min_samples = 2
params.differential_fdr = 0.05
params.robustness_iterations = 50
params.robustness_seed = 1729

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
    def differentialCheck = params.run_differential.toString() == 'true' ?
        "Rscript --vanilla -e 'if (!requireNamespace(\"DESeq2\", quietly=TRUE)) stop(\"DESeq2 is required for --run_differential\")'" : ''
    """
    ${differentialCheck}
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
        --technical-prefix-max-errors '${params.technical_prefix_max_errors}' ${params.exact20_rescue.toString() == 'true' ? '--exact20-rescue' : ''} --outdir .
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

process ANALYZE_REGIONS {
    label 'reporting'
    cpus 1
    publishDir "${params.outdir}/Analysis", mode: 'copy', pattern: '*.{tsv,bed,json}'
    publishDir "${params.outdir}/Analysis", mode: 'copy', pattern: 'PCA'
    publishDir "${params.outdir}/Analysis", mode: 'copy', pattern: 'Annotation'
    publishDir "${params.outdir}/Analysis", mode: 'copy', pattern: 'WindowOverlap'

    input:
    tuple val(samples), path(counts, arity: '1..*'), path(headers, arity: '1..*'), path(molecules, arity: '1..*')
    path design_files, arity: '1..*'
    path annotation_gtf, stageAs: 'annotation_reference/*'
    path annotation_environment

    output:
    path '*.{tsv,bed,json}', emit: tables
    path 'analysis_mqc.json', emit: multiqc
    path 'PCA', optional: true, emit: pca
    path 'PCA/pca_mqc.json', optional: true, emit: pca_multiqc
    path 'Annotation', optional: true, emit: annotation
    path 'WindowOverlap', optional: true, emit: window_overlap
    path 'WindowOverlap/window_upset_mqc.json', optional: true, emit: window_overlap_multiqc
    path 'Annotation/*_mqc.json', optional: true, emit: annotation_multiqc

    script:
    def sampleArgs = samples.collect { "'${it}'" }.join(' ')
    def countArgs = counts.collect { "'${it}'" }.join(' ')
    def headerArgs = headers.collect { "'${it}'" }.join(' ')
    def moleculeArgs = molecules.collect { "'${it}'" }.join(' ')
    def windows = params.run_windows.toString() == 'true' ? params.window_sizes : ''
    def pcaCommand = windows && params.run_pca.toString() == 'true' ?
        "env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python ${projectDir}/scripts/pca_windows.py --analysis-dir . --outdir PCA --windows ${windows.toString().split(',').join(' ')}" : ''
    def geneMinimumReps = params.gene_min_reps != null ? params.gene_min_reps : 2
    def annotationCommand = annotation_gtf ?
        "env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python ${projectDir}/scripts/annotate_regions.py --samples ${sampleArgs} --counts ${countArgs} --headers ${headerArgs} --design-file design.summary.json --gtf '${annotation_gtf}' --gtf-genome '${params.gtf_genome}' --genome '${params.genome}' --environment-file '${annotation_environment}' --promoter-upstream '${params.promoter_upstream}' --promoter-downstream '${params.promoter_downstream}' --gene-min-reps '${geneMinimumReps}' --gene-min-molecules '${params.gene_min_molecules}' --analysis-dir . --outdir Annotation" : ''
    def overlapCommand = windows && params.run_window_upset.toString() == 'true' ?
        "env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python ${projectDir}/scripts/window_upset.py --analysis-dir . --outdir WindowOverlap --windows ${windows.toString().split(',').join(' ')} --minimum-molecules '${params.window_min_molecules}' --minimum-replicates '${params.window_min_reps}' --iterations '${params.robustness_iterations}' --seed '${params.robustness_seed}'" : ''
    """
    python ${projectDir}/scripts/analyze_regions.py \\
        --samples ${sampleArgs} --counts ${countArgs} --headers ${headerArgs} --molecules ${moleculeArgs} \\
        --design-file design.summary.json --windows '${windows}'
    ${overlapCommand}
    ${pcaCommand}
    ${annotationCommand}
    """
}

process DIFFERENTIAL_DSB {
    label 'reporting'
    cpus 1
    publishDir "${params.outdir}/Analysis", mode: 'copy'

    input:
    path region_tables, arity: '1..*'
    path annotation_files

    output:
    path 'Differential', emit: results
    path 'Differential/differential_mqc.json', emit: multiqc

    script:
    def timecourse = params.differential_mode == 'timecourse'
    def contrastArgs = (timecourse ? params.timepoints : params.contrast).toString().split(',').collect { "'${it}'" }.join(' ')
    def differentialScript = timecourse ? 'timecourse_dsb.py' : 'differential_dsb.py'
    def modelArgs = timecourse ? "--timepoints ${contrastArgs}" : "--design '${params.design}' --contrast ${contrastArgs}"
    def robustnessArgs = timecourse ? '' : "--iterations '${params.robustness_iterations}' --seed '${params.robustness_seed}'"
    """
    env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
        python ${projectDir}/scripts/${differentialScript} \\
        --analysis-dir . --outdir Differential ${modelArgs} \\
        --minimum-count '${params.differential_min_count}' --minimum-samples '${params.differential_min_samples}' \\
        --fdr '${params.differential_fdr}' ${robustnessArgs}
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
    if (params.gtf) {
        def aliases = [hg19: 'hg19', grch37: 'hg19', hg38: 'hg38', grch38: 'hg38', mm10: 'mm10', grcm38: 'mm10', mm39: 'mm39', grcm39: 'mm39']
        if (!params.gtf_genome || !aliases[params.gtf_genome.toString().toLowerCase()] ||
            aliases[params.gtf_genome.toString().toLowerCase()] != aliases[params.genome.toString().toLowerCase()]) {
            error '--gtf requires --gtf_genome matching --genome'
        }
    } else if (params.gtf_genome) {
        error '--gtf_genome requires --gtf'
    }
    if (!(params.promoter_upstream.toString() ==~ /0|[1-9][0-9]*/) ||
        !(params.promoter_downstream.toString() ==~ /[1-9][0-9]*/)) {
        error '--promoter_upstream must be nonnegative and --promoter_downstream positive'
    }
    annotation_reference = params.gtf ? file(params.gtf, checkIfExists: true) : []
    ['gene_min_molecules', 'gene_min_reps'].each { key ->
        if (params[key] != null && !(params[key].toString() ==~ /[1-9][0-9]*/)) {
            error "--${key} must be a positive integer"
        }
    }
    if (params.run_peak_calling.toString() != 'false') {
        error 'MACS3 peak calling has been removed from main.nf; use window support and --run_window_upset instead'
    }
    ['run_windows', 'run_pca', 'run_window_upset', 'run_differential'].each { key ->
        if (!(params[key].toString() in ['true', 'false'])) { error "--${key} must be true or false" }
    }
    if (!(params.window_sizes.toString() ==~ /[1-9][0-9]*(,[1-9][0-9]*)*/)) {
        error '--window_sizes must be comma-separated positive integers'
    }
    ['window_min_molecules', 'window_min_reps'].each { key ->
        if (!(params[key].toString() ==~ /[1-9][0-9]*/) || params[key].toLong() > Integer.MAX_VALUE) {
            error "--${key} must be a positive integer"
        }
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
    if (!(params.differential_mode in ['contrast', 'timecourse'])) { error '--differential_mode must be contrast or timecourse' }
    if (params.run_differential.toString() == 'true') {
        if (params.design == 'unspecified') { error '--run_differential requires --design paired or unpaired' }
        if (params.differential_mode == 'timecourse') {
            if (params.design != 'paired' || params.contrast) { error 'Time course requires --design paired and no --contrast' }
            if (!(params.timepoints ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*(,[A-Za-z0-9][A-Za-z0-9_.-]*){2,}/)) {
                error '--timepoints requires >=3 comma-separated labels, baseline first'
            }
            def times = params.timepoints.toString().split(',')
            if (times.toList().toSet().size() != times.size()) { error '--timepoints contains duplicates' }
        } else {
            if (!(params.contrast ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*,[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error '--contrast must be NUMERATOR,DENOMINATOR with two explicit condition labels'
            }
            def contrastLabels = params.contrast.toString().split(',')
            if (contrastLabels[0] == contrastLabels[1]) { error 'Contrast conditions must differ' }
        }
        if (params.run_windows.toString() != 'true' && !params.gtf) { error '--run_differential requires --run_windows or --gtf' }
    }
    ['differential_min_count', 'differential_min_samples', 'robustness_iterations'].each { key ->
        if (!(params[key].toString() ==~ /[1-9][0-9]*/) || params[key].toLong() > Integer.MAX_VALUE) {
            error "--${key} must be a positive integer"
        }
    }
    if ((params.differential_min_samples as Integer) < 2 || (params.robustness_iterations as Integer) < 2) {
        error '--differential_min_samples and --robustness_iterations must be >=2'
    }
    if (!(params.robustness_seed.toString() ==~ /0|[1-9][0-9]*/)) { error '--robustness_seed must be nonnegative' }
    if (!(params.differential_fdr.toString() ==~ /[0-9.eE+-]+/) ||
        (params.differential_fdr as Double) <= 0 || (params.differential_fdr as Double) >= 1) {
        error '--differential_fdr must be between 0 and 1'
    }
    if (!(params.dedup_method in ['exact', 'directional'])) { error '--dedup_method must be exact or directional' }
    if (!(params.barcode_orientation in ['forward', 'reverse_complement'])) {
        error '--barcode_orientation must be forward or reverse_complement'
    }
    if (params.technical_prefix && !(params.technical_prefix ==~ /[ACGT]{12,}/)) {
        error '--technical_prefix must contain at least 12 A/C/G/T bases'
    }
    if (!(params.exact20_rescue.toString() in ['true', 'false'])) {
        error '--exact20_rescue must be true or false'
    }
    if (params.exact20_rescue.toString() == 'true' && params.technical_prefix != 'CCCTATAGTGAGTCGTAT') {
        error '--exact20_rescue requires --technical_prefix CCCTATAGTGAGTCGTAT'
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
            if (params.run_differential.toString() == 'true') {
                def contrasts = (params.differential_mode == 'timecourse' ? params.timepoints : params.contrast).toString().split(',')
                if (params.differential_mode == 'timecourse' && rows.collect { it.group }.toSet() != contrasts.toList().toSet()) {
                    error '--timepoints must match all manifest conditions exactly'
                }
                if (contrasts.any { group -> rows.count { it.group == group } < 3 }) {
                    error '--run_differential requires >=3 biological samples per contrast condition'
                }
            }
            if (params.run_windows.toString() == 'true' && params.run_window_upset.toString() == 'true' &&
                rows.groupBy { it.group }.any { group, members -> members.size() < (params.window_min_reps as Integer) }) {
                error '--window_min_reps exceeds the biological replicate count of a condition'
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
    if (params.run_windows.toString() == 'true' || params.gtf) {
        region_inputs = qc.condition_input
            .map { meta, summary, header, counts, bam, bai -> tuple(meta, counts, header) }
            .join(dedup.molecules)
        region_inputs = region_inputs
            .collect(flat: false)
            .map { rows ->
                def ordered = rows.sort { a, b -> a[0].sample <=> b[0].sample }
                tuple(ordered.collect { it[0].sample }, ordered.collect { it[1] },
                      ordered.collect { it[2] }, ordered.collect { it[3] })
            }
        analysis = ANALYZE_REGIONS(region_inputs, environment.design_metadata.flatten().collect(), annotation_reference, environment.ready)
        analysis_reports = analysis.multiqc.mix(analysis.pca_multiqc, analysis.annotation_multiqc, analysis.window_overlap_multiqc).flatten()
        if (params.run_differential.toString() == 'true') {
            differential = DIFFERENTIAL_DSB(analysis.tables, analysis.annotation.collect().ifEmpty([]))
            analysis_reports = analysis_reports.mix(differential.multiqc)
        }
    }
    MULTIQC(qc.condition_input.map { meta, summary, header, counts, bam, bai -> summary }.collect(),
            conditions.summary.collect(), qc.samtools_qc.mix(conditions.samtools_qc).flatten().collect(),
            prepared.preparation_multiqc.collect(), environment.design_metadata.flatten().collect(), analysis_reports.collect().ifEmpty([]))
}
