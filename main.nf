/*
 * pipeline input parameters
 */

params.csv = ""
params.outdir = "output"
params.adapter_seq = "GTGGAAAGGACGAAACACCG"
params.guide_ref = ""
params.screen_name = "mini_screen"
params.normalization = ""
params.control_guides = "$projectDir/assets/NO_FILE"
params.design_matrix = ""

params.min_counts_threshold = 50
params.sample_frac_with_min_counts = 0.80
params.sample_to_remove = ""

/*
 * PROCESSES
 */

process GUIDE_REF_TO_FASTA {
    tag "guide_ref"
    label 'process_single'
    publishDir "${params.outdir}/reference", enabled: false, mode:'copy'

    input:
    path gRNA_table

    output:
    path "guides.fasta", emit: gRNA_fasta

    script:
    """
    awk -F, '(NR > 1) {print ">"\$1"\\n"\$2}' ${gRNA_table} > guides.fasta
    """
}

process REF_INDEX {
    tag "guide_ref"
    label 'process_single'
    publishDir "${params.outdir}/reference", enabled: false, mode:'copy'
    container 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0'

    input:
    path gRNA_fasta

    output:
    path "bowtie2", emit: index

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads $task.cpus $gRNA_fasta bowtie2/${gRNA_fasta.baseName}
    """
}

process TRIM_ADAPTER {
    tag "${sample_id}"
    label 'process_single'
    publishDir "${params.outdir}/fastq_trimming", enabled: false, mode:'copy'
    container "https://depot.galaxyproject.org/singularity/cutadapt:4.3--py39hbf8eff0_0"

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${fastq_file.getSimpleName()}_trim.fastq.gz"), emit: trimmed_fastq

    script:
    """
       cutadapt -e 0.25 -O 15 -l 20 \
            -g ${params.adapter_seq} \
            -o ${fastq_file.getSimpleName()}_trim.fastq.gz \
            ${fastq_file}
    """
}

process ALIGNMENT {
    tag "${sample_id}"
    label 'process_single'
    publishDir "${params.outdir}/alignment", enabled: false, mode:'copy'
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0'

    input:
    path(bowtie2_index)
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    tuple val(sample_id), path("*.log"), emit: log

    script:
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    bowtie2 -5 0 -3 0 --norc -x \$INDEX -p $task.cpus -U ${fastq_file} 2> ${sample_id}.bowtie2.log \
        | samtools sort -@ $task.cpus -o ${sample_id}.bam -
    """
}

process MAGECK_COUNT {
    tag "mageck_count"
    label 'process_single'
    publishDir "${params.outdir}/mageck_count", enabled: true, mode:'copy'
    container "docker://davidliwei/mageck:latest"

    input:
    val sample_list
    path fastq_list
    path guide_reference_file
    path control_guides

    output:
    path("*.count_normalized.txt"), emit: norm_counts
    path("*.count.txt"), emit: raw_counts
    path("${params.screen_name}*")

    script:
    def ctrl_guides_opt = (params.normalization == "control") ? "--control-sgrna ${ control_guides }" : ""
    """
    mageck count -l ${guide_reference_file} \
                --fastq ${fastq_list} \
                --sample-label  ${ (sample_list instanceof List ? sample_list.join(',') : sample_list) } \
                --norm-method ${ params.normalization } ${ ctrl_guides_opt } \
                -n ${params.screen_name}
    """
}

process COUNT_FILTER {
    tag "filter_count"
    label 'process_single'
    publishDir "${params.outdir}/mageck_count_filtered", enabled: true, mode:'copy'
    container "docker://sddcunit/samplesheet_parsing:v1.0.0"

    input:
    path raw_count_table
    path norm_count_table

    output:
    path("*.count_normalized.filtered.txt"), emit: filter_norm_counts
    path("*.count.filtered.txt"), emit: filter_raw_counts

    """
    #!/usr/bin/env python3

    # import statements
    import pandas as pd

    # read raw and normalize count tables
    raw_counts_df = pd.read_csv("${ raw_count_table }", sep="\t", index_col=0)
    norm_counts_df = pd.read_csv("${ norm_count_table }", sep="\t", index_col=0, dtype="str")

    # optionally filter out columns
    columns_to_remove = "${ params.sample_to_remove }"
    if columns_to_remove:
        raw_counts_df = raw_counts_df.drop(labels = columns_to_remove.split(","), axis = 1)
        norm_counts_df = norm_counts_df.drop(labels = columns_to_remove.split(","), axis = 1)

    # select genes with less than ${params.min_counts_threshold} counts in
    # more than ${ params.sample_frac_with_min_counts } of the samples
    sample_columns = raw_counts_df.columns[1:]
    drop = (raw_counts_df[sample_columns] < ${params.min_counts_threshold}).sum(axis=1) / len(sample_columns) \
                > ${ params.sample_frac_with_min_counts }

    # print filtered tables to tsv files
    raw_counts_df[~drop].to_csv("${ raw_count_table.getSimpleName() }.count.filtered.txt", sep="\t", index=True)
    norm_counts_df[~drop].to_csv("${ norm_count_table.getSimpleName() }.count_normalized.filtered.txt", sep="\t", index=True)

    """
}

process MAGECK_MLE {
    tag "${design_id}"
    label 'process_single_highmem'
    publishDir "${params.outdir}/mageck_mle", enabled: true, mode:'copy'
    container "docker://davidliwei/mageck:latest"

    input:
    path counts
    tuple val(design_id), path(design_matrix)
    path control_guides

    output:
    path("${design_id}*")

    script:
    def ctrl_guides_opt = (params.normalization == "control") ? "--control-sgrna ${ control_guides }" : ""
    """
    mageck mle -k ${ counts } -d ${ design_matrix } \
               -n ${ design_id } --norm-method none ${ ctrl_guides_opt }
    """
}

/*
 * WORKFLOW
 */

workflow {

// prepare guide index for alingment
GUIDE_REF_TO_FASTA( params.guide_ref )
REF_INDEX( GUIDE_REF_TO_FASTA.out.gRNA_fasta )

// read list of fastq files from input csv file
Channel
    .fromPath(params.csv)
    .splitCsv(header: true)
    .map { it -> [it.sample_id, it.fastq] }
    .groupTuple()
    .set { fastq_files_ch }

// trim the fastq files with cutadapt
TRIM_ADAPTER( fastq_files_ch )

// align FASTQ files to bwa index
ALIGNMENT(
    REF_INDEX.out.index,
    TRIM_ADAPTER.out.trimmed_fastq
)

ALIGNMENT.out.bam
    .multiMap {
        sample_list: it[0]
        fastq_list:  it[1]
    }
    .set { aligned_bam_ch }

MAGECK_COUNT(
    aligned_bam_ch.sample_list.collect(),
    aligned_bam_ch.fastq_list.collect(),
    params.guide_ref,
    params.control_guides
)

COUNT_FILTER(
    MAGECK_COUNT.out.raw_counts,
    MAGECK_COUNT.out.norm_counts,
)

Channel
    .fromPath( params.design_matrix.tokenize(",") )
    .map { it -> [it.getSimpleName(), it] }
    .set{ design_matrix_ch }

MAGECK_MLE(
    COUNT_FILTER.out.filter_norm_counts,
    design_matrix_ch,
    params.control_guides
)

}
