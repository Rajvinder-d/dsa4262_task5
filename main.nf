nextflow.enable.dsl = 2

def inputPaths(value) {
    if (!value) return []
    def entries = value instanceof List ? value : value.toString().split(',').toList()
    return entries.collect { entry ->
        def p = file(entry.toString().trim(), checkIfExists: true)
        if (p instanceof List || !java.nio.file.Files.isRegularFile(p)) {
            error 'Provide explicit file paths, separated by commas; not folders or glob patterns.'
        }
        p.toAbsolutePath().normalize()
    }
}

def shellQuote(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

// Reference preparation is cached across the two annotation scenarios.
process PREPARE_REFERENCE {
    tag 'genome'
    cpus 3
    publishDir "${params.outdir}/reference/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/reference", mode: 'symlink', pattern: 'genome.*'

    input:
    path fasta_input, stageAs: 'source/*'

    output:
    path 'genome.fa', emit: fasta
    path 'genome.fa.fai', emit: fai
    path 'genome.directRNA.k14.I1G.mmi', emit: drna_index
    path 'genome.cDNA.k15.I1G.mmi', emit: cdna_index
    path '*.log', emit: logs

    script:
    """
    case ${shellQuote(fasta_input)} in
        *.gz) gzip -dc ${shellQuote(fasta_input)} > genome.fa ;;
        *) cp ${shellQuote(fasta_input)} genome.fa ;;
    esac
    samtools faidx genome.fa
    minimap2 -x splice -k14 -I 1G -t ${task.cpus} \
        -d genome.directRNA.k14.I1G.mmi genome.fa 2> directRNA.index.log
    minimap2 -x splice -k15 -I 1G -t ${task.cpus} \
        -d genome.cDNA.k15.I1G.mmi genome.fa 2> cDNA.index.log
    """
}

process MINIMAP2_ALIGN {
    tag "${sample}"
    cpus 4
    publishDir "${params.outdir}/alignment/sam", mode: 'symlink', pattern: '*.sam'
    publishDir "${params.outdir}/alignment/logs", mode: 'copy', pattern: '*.minimap2.*'

    input:
    tuple val(sample), val(technology), path(reads)
    path drna_index
    path cdna_index

    output:
    tuple val(sample), path("${sample}.sam"), emit: sam
    path "${sample}.minimap2.*", emit: logs

    script:
    def rna_options = technology == 'directRNA' ? '-uf -k14' : ''
    def index = technology == 'directRNA' ? drna_index : cdna_index
    """
    mkdir -p tmp
    minimap2 --version > '${sample}.minimap2.version.txt'
    minimap2 -ax splice ${rna_options} -t ${task.cpus} \
        -K 100M --cap-kalloc 500m --split-prefix=tmp/${sample} \
        ${shellQuote(index)} ${shellQuote(reads)} \
        > '${sample}.sam' 2> '${sample}.minimap2.log'
    """
}

process SAM_TO_BAM {
    tag "${sample}"
    cpus 2
    publishDir "${params.outdir}/alignment/bam", mode: 'copy', pattern: '*.bam'
    publishDir "${params.outdir}/alignment/logs", mode: 'copy', pattern: '*.samtools.*'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}.bam"), emit: bam
    path "${sample}.samtools.*", emit: logs

    script:
    """
    samtools --version > '${sample}.samtools.version.txt'
    samtools view -@ 1 -b -o '${sample}.bam' '${sam}' \
        2> '${sample}.samtools.log'
    """
}

process BAMBU_DISCOVERY {
    tag "${annotation_mode}"
    cpus 2
    publishDir { "${params.outdir}/scenarios/${annotation_mode}" }, mode: 'copy'

    input:
    path bams
    path fasta, stageAs: 'genome.fa'
    path fai, stageAs: 'genome.fa.fai'
    path gtf, stageAs: 'reference.gtf'
    val annotation_mode
    val rscript
    path r_script

    output:
    path 'bambu', emit: results

    script:
    def bam_args = bams.collect { bam -> "'${bam}'" }.join(' ')
    def gtf_arg = annotation_mode == 'annotated' ? "'${gtf}'" : 'NONE'
    """
    mkdir -p bambu
    env \
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        ${shellQuote(rscript)} --vanilla '${r_script}' \
        '${fasta}' ${gtf_arg} '${annotation_mode}' ${task.cpus} ${bam_args} \
        2>&1 | tee bambu/run.log
    """
}

process QC_REPORT {
    tag "${annotation_mode}"
    cpus 2
    publishDir { "${params.outdir}/scenarios/${annotation_mode}" }, mode: 'copy'

    input:
    path se_file, stageAs: 'final_se.rds'
    path fastqs, stageAs: 'reads/*'
    path bams, stageAs: 'bams/*'
    val annotation_mode
    val rscript
    path qc_r
    path qc_py

    output:
    path 'qc', emit: report

    script:
    """
    mkdir -p qc
    env \
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        ${shellQuote(rscript)} --vanilla '${qc_r}' '${se_file}' qc \
        2>&1 | tee qc/r.log
    python3 '${qc_py}' --reads-dir reads --bam-dir bams --outdir qc \
        --annotation-mode '${annotation_mode}' \
        --min-primary-reads ${params.qc_min_primary_reads} \
        --min-mapping-pct ${params.qc_min_mapping_pct} \
        --short-length ${params.qc_short_length} \
        --max-short-pct ${params.qc_max_short_pct} \
        2>&1 | tee qc/run.log
    """
}

workflow {
    workflow.onComplete = {
        log.info "Workflow: ${workflow.success ? 'SUCCESS' : 'FAILED'}; annotation mode: ${params.without_annotations ? 'unannotated' : 'annotated'}"
        log.info "Outputs: ${params.outdir}"
    }
    if (!params.outdir) error 'Provide --outdir /path/to/output_directory'
    def output_root = file(params.outdir).toAbsolutePath().normalize()
    def saved_path = output_root.resolve('run_inputs.json')
    def saved = [:]
    if (java.nio.file.Files.exists(saved_path)) {
        if (!workflow.resume) {
            error 'This output directory already contains a run. Choose a new --outdir for a fresh run, or use -resume.'
        }
        saved = new groovy.json.JsonSlurper().parseText(saved_path.text)
        if (saved.session_id != workflow.sessionId.toString()) {
            error "The selected cache belongs to another run. Use -resume ${saved.session_id} for this output directory."
        }
    }

    def fasta_arg = params.fasta ?: saved.fasta
    def gtf_arg = params.gtf ?: saved.gtf
    def drna_arg = params.direct_rna ?: saved.direct_rna
    def cdna_arg = params.cdna ?: saved.cdna
    def rscript_arg = params.rscript ?: saved.rscript ?: 'Rscript'
    def annotation_mode = params.without_annotations ? 'unannotated' : 'annotated'
    if (!fasta_arg || (!drna_arg && !cdna_arg)) {
        error 'First run requires --fasta and explicit FASTQs via --direct_rna and/or --cdna. A resumed run can load these from --outdir/run_inputs.json.'
    }
    if (annotation_mode == 'annotated' && !gtf_arg) error 'Annotated discovery requires --gtf /path/to/annotation.gtf'
    def genome_inputs = inputPaths(fasta_arg)
    if (genome_inputs.size() != 1) error 'Supply exactly one genome FASTA with --fasta'
    def fasta_file = genome_inputs[0]
    // No GTF is opened or staged in the unannotated scenario.
    def annotations = annotation_mode == 'annotated' ? inputPaths(gtf_arg) : []
    if (annotation_mode == 'annotated' && annotations.size() != 1) error 'Supply exactly one annotation GTF with --gtf'

    def drna_files = inputPaths(drna_arg)
    def cdna_files = inputPaths(cdna_arg)
    def sample_definitions = [[drna_files, 'directRNA'], [cdna_files, 'cDNA']].collectMany { entry ->
        entry[0].collect { reads ->
            def sample = reads.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            if (sample == reads.name) error "Expected FASTQ (.fastq, .fq, optionally .gz): ${reads}"
            if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error 'FASTQ basenames must contain only letters, digits, dot, underscore or hyphen.'
            tuple(sample, entry[1], reads)
        }
    }
    sample_definitions = sample_definitions.sort { a, b -> a[0] <=> b[0] }
    if (sample_definitions.collect { sample -> sample[0] }.unique().size() != sample_definitions.size()) {
        error 'FASTQ basenames must be unique across all samples and technology groups.'
    }
    def read_files = sample_definitions.collect { sample -> sample[2] }

    // Persist the first run's resolved inputs. The second command needs only
    // --outdir, --without_annotations and -resume, from the same launch directory.
    java.nio.file.Files.createDirectories(output_root)
    if (!java.nio.file.Files.exists(saved_path)) {
        saved_path.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
            fasta: fasta_file.toString(),
            gtf: annotation_mode == 'annotated' ? annotations[0].toString() : null,
            direct_rna: drna_files.collect { reads -> reads.toString() },
            cdna: cdna_files.collect { reads -> reads.toString() },
            rscript: rscript_arg.toString(),
            session_id: workflow.sessionId.toString()
        ])) + '\n'
    }
    log.info "Input: ${read_files.size()} FASTQs; annotation mode: ${annotation_mode}"
    def samples = channel.fromList(sample_definitions)
    PREPARE_REFERENCE(fasta_file)
    MINIMAP2_ALIGN(samples, PREPARE_REFERENCE.out.drna_index,
                   PREPARE_REFERENCE.out.cdna_index)
    SAM_TO_BAM(MINIMAP2_ALIGN.out.sam)
    def qc_bams = SAM_TO_BAM.out.bam.map { sample, bam -> bam }
        .collect()
        .map { files -> files.sort { a, b -> a.name <=> b.name } }
    BAMBU_DISCOVERY(
        qc_bams,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai,
        annotation_mode == 'annotated' ? annotations[0] : [],
        annotation_mode,
        rscript_arg,
        file("${projectDir}/bambu.R", checkIfExists: true)
    )
    def qc_se = BAMBU_DISCOVERY.out.results.map { directory -> directory.resolve('se.rds') }
    QC_REPORT(
        qc_se,
        read_files,
        qc_bams,
        annotation_mode,
        rscript_arg,
        file("${projectDir}/qc.R", checkIfExists: true),
        file("${projectDir}/qc.py", checkIfExists: true)
    )
}
