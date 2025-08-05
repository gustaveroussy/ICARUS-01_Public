rule table_summary:
    input:
        config["data"]["inputs"]["daisy_pfs_response"],
        config["data"]["inputs"]["wes_ids"],
        config["data"]["inputs"]["wes_qc_annotated"],
        config["data"]["inputs"]["slide_ids"],
        config["data"]["inputs"]["geomx_ids"]
    output:
        config["data"]["outputs"]["summary_data"]
    benchmark:
        "%s/table_summary.tsv" % B_FOLDER
    log:
        "%s/table_summary.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    threads: 1
    shell:
        """
        python -u scripts/01.1_table_summary.py \
            --daisy_pfs_response {input[0]} \
            --wes_ids {input[1]} \
            --wes_qc_annotated {input[2]} \
            --slide_ids {input[3]} \
            --geomx_ids {input[4]} \
            --output_summary {output}  &> {log}
        """


rule table_ega_metadata:
    input:
        config["data"]["outputs"]["summary_data"],
        config["data"]["inputs"]["wes_ids"],
        config["data"]["inputs"]["ids_anonym"],
        config["data"]["inputs"]["sam_used"]
    output:
        ega_samples = config["data"]["outputs"]["ega_samples"],
        ega_metadata = config["data"]["outputs"]["ega_metadata"]
    benchmark:
        "%s/table_ega_metadata.tsv" % B_FOLDER
    log:
        "%s/table_ega_metadata.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    threads: 1
    shell:
        """
        Rscript scripts/01.2_table_ega_metadata.R \
            --sum {input[0]} \
            --ids {input[1]} \
            --ids_anonym {input[2]} \
            --sam_used {input[3]} \
            --output_ega_samples {output.ega_samples} \
            --output_ega_metadata {output.ega_metadata} \
            --log {log}
        """
