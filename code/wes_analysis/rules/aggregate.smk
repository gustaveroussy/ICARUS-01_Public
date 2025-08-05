rule aggregate_alterations_across_modalities:
    input:
        config["data"]["inputs"]["ids_curated"],
        config["data"]["inputs"]["cln_curated"],
        config["data"]["inputs"]["cna_annotated"],
        config["data"]["inputs"]["cna_pass"],
        config["data"]["inputs"]["mut_annotated"],
        config["data"]["inputs"]["drug"]
    output:
        config["data"]["outputs"]["alterations_best"],
        config["data"]["outputs"]["alterations_all"]
    benchmark:
        "%s/aggregate_alterations_across_modalities.tsv" % B_FOLDER
    log:
        "%s/aggregate_alterations_across_modalities.log" % L_FOLDER
    conda:
        config["setup"]["r-icarus"]
    threads: 1
    shell:
        """
        Rscript scripts/01.1_aggregate_alterations_across_modalities.R \
            --ids {input[0]} \
            --cln {input[1]} \
            --cna {input[2]} \
            --cna_pass {input[3]} \
            --mut {input[4]} \
            --drug {input[5]} \
            --output_best {output[0]} \
            --output_all {output[1]} \
            --log {log}
        """


rule aggregate_alterations_across_cohorts:
    input:
        config["data"]["outputs"]["alterations_best"]
    output:
        config["data"]["outputs"]["alterations"]
    benchmark:
        "%s/aggregate_alterations_across_cohorts.tsv" % B_FOLDER
    log:
        "%s/aggregate_alterations_across_cohorts.log" % L_FOLDER
    conda:
        config["setup"]["r-icarus"]
    params:
        cohorts=["icarus"]
    threads: 1
    shell:
        """
        Rscript scripts/01.2_aggregate_alterations_across_cohorts.R \
            --cohorts {params.cohorts} \
            --alts {input[0]} \
            --output {output} \
            --log {log}
        """


rule detail_alterations_gene:
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["data"]["inputs"]["mut_pass"],
        config["data"]["inputs"]["cna_pass"]
    output:
        "%s/alterations/alterations_{gene}_all.xlsx" % R_FOLDER
    benchmark:
        "%s/detail_alterations_gene_{gene}.tsv" % B_FOLDER
    log:
        "%s/detail_alterations_gene_{gene}.log" % L_FOLDER
    conda:
        config["setup"]["py-icarus"]
    params:
        cohorts=["daisy"]
    threads: 1
    shell:
        """
        python -u scripts/01.3_detail_alterations_gene.py \
            --cln {input[0]} \
            --alt {input[1]} \
            --mut {input[2]} \
            --cna {input[3]} \
            --gen {wildcards.gene} \
            --output {output} &> {log}
        """
