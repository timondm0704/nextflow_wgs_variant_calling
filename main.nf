// nextflow_wgs_variant_calling
// Paired-end WGS: trim reads -> align -> check coverage (12x minimum) -> call variants

params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.reference   = "${projectDir}/data/reference/SHH_locus_hg38.fa"
params.min_coverage   = 12
params.outdir      = "results"


// Step 1: trim low-quality 3' ends (below Q20) and drop reads shorter than 50 bp
process TRIM {
    tag "${id}"

    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_R1.trimmed.fastq.gz"), path("${id}_R2.trimmed.fastq.gz")

    script:
    """
    cutadapt -q 20 -m 50 \\
        -o ${id}_R1.trimmed.fastq.gz -p ${id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} > ${id}.cutadapt.log
    """
}


// Step 2: index the reference so reads can be aligned to it
process INDEX {
    input:
    path fasta

    output:
    path "ref"

    script:
    """
    mkdir ref
    cp ${fasta} ref/genome.fa
    bwa index ref/genome.fa
    samtools faidx ref/genome.fa
    """
}


// Step 3: align each individual's read pairs, then sort and index the BAM
process ALIGN {
    tag "${id}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(id), path(r1), path(r2)
    path ref

    output:
    tuple val(id), path("${id}.bam"), path("${id}.bam.bai")

    script:
    """
    bwa mem -R "@RG\\tID:${id}\\tSM:${id}" ref/genome.fa ${r1} ${r2} | samtools sort -o ${id}.bam -
    samtools index ${id}.bam
    """
}


// Step 4: mean coverage per individual, PASS if >= 12x
process COVERAGE {
    tag "${id}"
    publishDir "${params.outdir}/coverage", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path(bam), path(bai), env(STATUS), emit: result
    path "${id}.coverage.txt",                       emit: report

    script:
    """
    COV=\$(samtools coverage ${bam} | awk 'NR == 2 { print \$7 }')
    STATUS=\$(awk -v d=\$COV -v m=${params.min_coverage} 'BEGIN { print (d >= m ? "PASS" : "FAIL") }')
    echo "${id} mean_coverage=\$COV \$STATUS" > ${id}.coverage.txt
    """
}


// Step 5: call SNPs jointly across all individuals that passed
process CALL_VARIANTS {
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path bams
    path bais
    path ref

    output:
    path "cohort.vcf.gz"

    script:
    """
    bcftools mpileup -f ref/genome.fa ${bams} | bcftools call -mv -Oz -o cohort.vcf.gz
    """
}


workflow {
    // Read the sample sheet: one row per individual
    samples = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id,
                            file("${projectDir}/${row.fastq_1}"),
                            file("${projectDir}/${row.fastq_2}")) }

    ref = INDEX(file(params.reference))
    TRIM(samples)
    ALIGN(TRIM.out, ref)
    COVERAGE(ALIGN.out)

    // Keep only individuals with enough coverage
    passed = COVERAGE.out.result.filter { id, bam, bai, status -> status == 'PASS' }

    CALL_VARIANTS(passed.map { it[1] }.collect(), passed.map { it[2] }.collect(), ref)
}
