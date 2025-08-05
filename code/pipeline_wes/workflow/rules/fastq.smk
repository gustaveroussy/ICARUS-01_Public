if not config["general"]["automatic_irods_check"]["disable"]:
    rule irods_get_content:
        log:
            "%s/fastq/irods_get_content.log" % L_FOLDER
        benchmark:
            "%s/fastq/irods_get_content.log" % B_FOLDER
        params:
            project_names = config["general"]["automatic_irods_check"]["project_names"]
        output:
            irods = "%s/irods/irods_get_content.tsv" % R_FOLDER
        resources:
            queue = "shortq"
        threads: 1
        shell:
            """
            module unload anaconda3
            pns=( {params.project_names} )
            for pn in "${{pns[@]}}"; do
                imeta qu -C projectName like "%${{pn}}%" >> {output.irods} 2> {log}
            done
            """

    rule irods_check_status:
        log:
            "%s/fastq/irods_check_status.log" % L_FOLDER
        benchmark:
            "%s/fastq/irods_check_status.log" % B_FOLDER
        conda:
            "../envs/python.yaml"
        input:
            table = config["general"]["samples"],
            irods = "%s/irods/irods_get_content.tsv" % R_FOLDER
        output:
            status = "%s/irods/irods_check_status.pass" % R_FOLDER
        resources:
            queue = "shortq"
        threads: 1
        shell:
            """
            python workflow/scripts/01.1_irods_check_status.py \
                --table {input.table} \
                --irods {input.irods} \
                --output {output.status} &> {log}
            """
else:
    rule irods_check_status:
        log:
            "%s/fastq/irods_check_status.log" % L_FOLDER
        benchmark:
            "%s/fastq/irods_check_status.log" % B_FOLDER
        output:
            status = "%s/irods/irods_check_status.pass" % R_FOLDER
        resources:
            queue = "shortq"
        threads: 1
        shell:
            """
            touch {output}
            """


rule irods_get_fastq:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples_redo_bam])
    log:
        "%s/fastq/irods_get_fastq_{sample}.log" % L_FOLDER
    benchmark:
        "%s/fastq/irods_get_fastq_{sample}.log" % B_FOLDER
    input:
        table = config["general"]["samples"],
        check = "%s/irods/irods_check_status.pass" % R_FOLDER
    params:
        dir = "%s/data/fastq" % R_FOLDER,
        irods_1 = lambda w: get_fastqs_irods(w)["r1"],
        irods_2 = lambda w: get_fastqs_irods(w)["r2"],
        local_1 = lambda w: get_fastqs_local(w)["r1"],
        local_2 = lambda w: get_fastqs_local(w)["r2"],
        name_1 = lambda w: get_fastqs_names(w)["r1"],
        name_2 = lambda w: get_fastqs_names(w)["r2"]
    output:
        r1 = temp("%s/data/fastq/{sample}_R1.fastq.gz" % R_FOLDER),
        r2 = temp("%s/data/fastq/{sample}_R2.fastq.gz" % R_FOLDER)
    resources:
        queue = "shortq",
        time_min = 120
    threads: 1
    shell:
        """
        module unload anaconda3
        bash workflow/scripts/01.2_irods_get_fastq.sh \
            -d {params.dir} \
            -a {params.irods_1} \
            -b {params.irods_2} \
            -x {params.local_1} \
            -y {params.local_2} \
            -m {params.name_1} \
            -n {params.name_2} &> {log}
        """
