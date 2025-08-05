# Align the fastq  file to the reference fastq file and sort with samtools to create the bam file
# See http://bio-bwa.sourceforge.net/bwa.shtml
rule bwa_mem:
    input:
        config["ref"]["fasta"],
        '%s/data/fastp/{sample}_R1.fastq.gz' % R_FOLDER,
        '%s/data/fastp/{sample}_R2.fastq.gz' % R_FOLDER
    output:
        temp("%s/mapping/{sample}.bam" % R_FOLDER)
    benchmark:
        "%s/mapping/bwa_mem/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bwa_mem/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        extra="-M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:truseq'"
    threads: 8
    resources:
        queue = "shortq",
        mem_mb = 50000,
        time_min = 360
    shell:
        """(bwa mem -t {threads} {params.extra} {input} \
            | samtools view -@ {threads} -b - \
            | samtools sort -@ {threads} - -o {output}) 2> {log}"""


# Remove duplicated reads from the bam file with MarkDuplicates (picard)
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
rule mark_duplicates:
    input:
        "%s/mapping/{sample}.bam" % R_FOLDER
    output:
        bam=temp("%s/mapping/{sample}.nodup.bam" % R_FOLDER),
        bai=temp("%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER),
        metrics="%s/qc/mark_duplicates/{sample}.metrics.txt" % R_FOLDER
    benchmark:
        "%s/mapping/mark_duplicates/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/mark_duplicates/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    params:
        java="'-Xmx120g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"],
        remove_duplicates=config["params"]["picard"]["MarkDuplicates"]["remove_duplicates"]
    resources:
        queue="shortq",
        mem_mb=120000,
        time_min=240
    shell:
        """
        gatk --java-options {params.java} MarkDuplicates \
            --INPUT {input} \
            --REMOVE_DUPLICATES {params.remove_duplicates} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} &> {log}
        samtools index \
            -@ {threads} \
            -b {output.bam} {output.bai} &>> {log}
        """


# Recalibrate base qualities using gatk BaseRecalibrator
rule bqsr_step_1:
    input:
        ref=config["ref"]["fasta"],
        ref_dict=config["ref"]["dict"],
        known_sites=config["params"]["gatk"]["BaseRecalibrator"]["known_sites"],
        bam="%s/mapping/{sample}.nodup.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_sample(w, file="bed_padded")
    output:
        temp("%s/mapping/{sample}_bqsr.table" % R_FOLDER)
    benchmark:
        "%s/mapping/bqsr_step_1/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bqsr_step_1/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=50000,
        time_min=120
    shell:
        """
        gatk --java-options {params.java} BaseRecalibrator \
            --reference {input.ref} \
            --known-sites {input.known_sites} \
            --intervals {input.bed} \
            --input {input.bam} \
            --output {output} &> {log}
        """


# Recalibrate base qualities using gatk ApplyBQSR
rule bqsr_step_2:
    input:
        bam="%s/mapping/{sample}.nodup.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER,
        ref=config["ref"]["fasta"],
        bqsr_recal_file="%s/mapping/{sample}_bqsr.table" % R_FOLDER
    output:
       temp("%s/mapping/{sample}.nodup.recal.beforeReformat.bam" % R_FOLDER)
    benchmark:
        "%s/mapping/bqsr_step_2/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bqsr_step_2/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=%s'" % config["general"]["tmpdir"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=15000,
        time_min=180
    shell:
        """
        gatk --java-options {params.java} ApplyBQSR \
            --create-output-bam-index false \
            --reference {input.ref} \
            --bqsr-recal-file {input.bqsr_recal_file} \
            --input {input.bam} \
            --output {output} &> {log}
        """


# Reformat bam files after base quality score recalibration with samtools in order to save space
# and generate bam index with samtools after recalibration
rule bqsr_postprocess:
    input:
        "%s/mapping/{sample}.nodup.recal.beforeReformat.bam" % R_FOLDER
    output:
        bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
    benchmark:
        "%s/mapping/bqsr_postprocess/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bqsr_postprocess/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 4
    resources:
        queue="shortq",
        mem_mb=10000,
        time_min=60
    shell:
        """
        samtools view \
            -@ {threads} \
            -b \
            -h \
            -o {output.bam} {input} &> {log} && \
        samtools index \
            -@ {threads} \
            -b {output.bam} {output.bai}  &>> {log}
        """


# Copy stored BAM file from archive
if len(samples_copy_bam)>0:
    rule archive_copy_bam:
        wildcard_constraints:
            sample = "|".join([re.escape(x) for x in samples_copy_bam])
        output:
            bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
            bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
        benchmark:
            "%s/mapping/archive_copy_bam/{sample}.tsv" % B_FOLDER
        log:
            "%s/mapping/archive_copy_bam/{sample}.log" % L_FOLDER
        threads: 1
        params:
            bam_archive="%s/results/mapping/{sample}.nodup.recal.bam" % config["general"]["archive_folder"],
            bai_archive="%s/results/mapping/{sample}.nodup.recal.bam.bai" % config["general"]["archive_folder"]
        resources:
            local_load=1
        shell:
            """
            cp {params.bai_archive} {output.bai}
            cp {params.bam_archive} {output.bam}
            """


# Merge multiple normal bam files into one in order to use for unmatched normal
# CNV calling.
rule normal_pool_bam:
    input:
        expand("%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER, sample=nsamples_normal_pool)
    output:
        bam="%s/mapping/%s.nodup.recal.bam" % (R_FOLDER, nsample_normal_pool),
        bai="%s/mapping/%s.nodup.recal.bam.bai" % (R_FOLDER, nsample_normal_pool)
    benchmark:
        "%s/mapping/normal_pool_bam/%s.tsv" % (B_FOLDER, nsample_normal_pool)
    log:
        "%s/mapping/normal_pool_bam/%s.log" % (L_FOLDER, nsample_normal_pool)
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="mediumq",
        mem_mb=40000,
        time_min=240
    shell:
        """
        samtools merge {output.bam} {input} && \
        samtools index \
            -@ {threads} \
            -b {output.bam} {output.bai} &> {log}
        """
