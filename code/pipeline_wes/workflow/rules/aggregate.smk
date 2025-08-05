####
#### Somatic mutations ####
####

# Extract a tab-delimited format file from the VCF with minimal information about the filters applied on variants.
rule somatic_maf_filters_aggregate:
    input:
        expand("%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        temp("%s/aggregate/somatic_maf/somatic_calls_filters_prefinal.tsv.gz" % R_FOLDER)
    benchmark:
        "%s/aggregate/somatic_maf/somatic_maf_filters_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_maf/somatic_maf_filters_aggregate.log" % L_FOLDER
    threads: 2
    resources:
        queue="shortq",
        mem_mb=1000,
        time_min=60
    shell:
        """
        zcat {input} | sed -n '1p;/^CHROM/ !p' | gzip > {output} 2> {log}
        """


# Aggregate all somatic MAF tables.
if config["params"]["vep"]["annotate_all"]:
    rule somatic_maf_aggregate:
        input:
            expand("%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
                get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
        output:
            all="%s/aggregate/somatic_maf/somatic_calls_prefinal.all.maf.gz" % R_FOLDER,
            subset=temp("%s/aggregate/somatic_maf/somatic_calls_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_aggregate.log" % L_FOLDER
        threads: 1
        params:
            column_subset="FILTER",
            values_subset=["PASS", "common_variant"]
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --column_subset {params.column_subset} \
                --values_subset {params.values_subset} \
                --keep_header \
                --output {output.all} \
                --output_subset {output.subset} &> {log}
            """


    # Aggregate all somatic civic-annotated MAF tables.
    rule somatic_maf_civic_aggregate:
        input:
            lambda w: get_input_concatenate(w, typ="maf", db="civic")
        output:
            all="%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.all.maf.gz" % R_FOLDER,
            subset=temp("%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_civic_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_civic_aggregate.log" % L_FOLDER
        threads: 1
        params:
            column_subset="FILTER",
            values_subset=["PASS", "common_variant"]
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --column_subset {params.column_subset} \
                --values_subset {params.values_subset} \
                --keep_header \
                --output {output.all} \
                --output_subset {output.subset} &> {log}
            """


    # Aggregate all somatic oncokb-annotated MAF tables.
    rule somatic_maf_oncokb_aggregate:
        input:
            lambda w: get_input_concatenate(w, typ="maf", db="oncokb")
        output:
            all="%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.all.maf.gz" % R_FOLDER,
            subset=temp("%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_oncokb_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_oncokb_aggregate.log" % L_FOLDER
        threads: 1
        params:
            column_subset="FILTER",
            values_subset=["PASS", "common_variant"]
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --column_subset {params.column_subset} \
                --values_subset {params.values_subset} \
                --keep_header \
                --output {output.all} \
                --output_subset {output.subset} &> {log}
            """

    # Aggregate oncokb and civic mutation annotations.
    rule somatic_maf_all_union_ann:
        input:
            civ="%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.all.maf.gz" % R_FOLDER,
            okb="%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.all.maf.gz" % R_FOLDER
        output:
            temp("%s/aggregate/somatic_maf/somatic_calls_union_ann_prefinal.all.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_all_union_ann.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_all_union_ann.log" % L_FOLDER
        resources:
            partition="cpu_short",
            mem_mb=8000,
            time="00:15:00"
        threads: 1
        shell:
            """
            python -u workflow/scripts/06.2_concatenate_annotations.py \
                --civ {input.civ} \
                --okb {input.okb} \
                --cat maf \
                --output {output} &> {log}
            """


    # Apply last filtering steps
    rule somatic_maf_all_last_filtering:
        input:
            maf="%s/aggregate/somatic_maf/somatic_calls_prefinal.all.maf.gz" % R_FOLDER,
            maf_ann="%s/aggregate/somatic_maf/somatic_calls_union_ann_prefinal.all.maf.gz" % R_FOLDER,
            maf_civ="%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.all.maf.gz" % R_FOLDER,
            maf_okb="%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.all.maf.gz" % R_FOLDER
        output:
            maf="%s/aggregate/somatic_maf/somatic_calls.all.maf.gz" % R_FOLDER,
            maf_ann="%s/aggregate/somatic_maf/somatic_calls_union_ann.all.maf.gz" % R_FOLDER,
            maf_civ="%s/aggregate/somatic_maf/somatic_calls_civic.all.maf.gz" % R_FOLDER,
            maf_okb="%s/aggregate/somatic_maf/somatic_calls_oncokb.all.maf.gz" % R_FOLDER
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf_all_last_filtering/somatic_maf_all_last_filtering.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf_all_last_filtering/somatic_maf_all_last_filtering.log" % L_FOLDER
        resources:
            partition="cpu_short",
            mem_mb=24000,
            time="00:45:00"
        threads: 1
        shell:
            """
            python -u workflow/scripts/03.3_maf_all_last_filtering_steps.py \
                --inp_maf {input.maf} \
                --inp_maf_ann {input.maf_ann} \
                --inp_maf_civ {input.maf_civ} \
                --inp_maf_okb {input.maf_okb} \
                --out_maf {output.maf} \
                --out_maf_ann {output.maf_ann} \
                --out_maf_civ {output.maf_civ} \
                --out_maf_okb {output.maf_okb} &> {log}
            """


else:
    rule somatic_maf_aggregate:
        input:
            expand("%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
                get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
        output:
            temp("%s/aggregate/somatic_maf/somatic_calls_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_aggregate.log" % L_FOLDER
        threads: 1
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --keep_header \
                --output {output} &> {log}
            """


    # Aggregate all somatic civic-annotated MAF tables.
    rule somatic_maf_civic_aggregate:
        input:
            lambda w: get_input_concatenate(w, typ="maf", db="civic")
        output:
            temp("%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_civic_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_civic_aggregate.log" % L_FOLDER
        threads: 1
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --keep_header \
                --output {output} &> {log}
            """


    # Aggregate all somatic oncokb-annotated MAF tables.
    rule somatic_maf_oncokb_aggregate:
        input:
            lambda w: get_input_concatenate(w, typ="maf", db="oncokb")
        output:
            temp("%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.maf.gz" % R_FOLDER)
        conda:
            "../envs/python.yaml"
        benchmark:
            "%s/aggregate/somatic_maf/somatic_maf_oncokb_aggregate.tsv" % B_FOLDER
        log:
            "%s/aggregate/somatic_maf/somatic_maf_oncokb_aggregate.log" % L_FOLDER
        threads: 1
        resources:
            queue="shortq",
            mem_mb=16000,
            time_min=60
        shell:
            """
            python -u workflow/scripts/06.1_concatenate_tables.py \
                --files {input} \
                --keep_header \
                --output {output} &> {log}
            """


# Aggregate oncokb and civic mutation annotations.
rule somatic_maf_union_ann:
    input:
        civ="%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.maf.gz" % R_FOLDER,
        okb="%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.maf.gz" % R_FOLDER
    output:
        temp("%s/aggregate/somatic_maf/somatic_calls_union_ann_prefinal.maf.gz" % R_FOLDER)
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/somatic_maf/somatic_maf_union_ann.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_maf/somatic_maf_union_ann.log" % L_FOLDER
    resources:
        partition="cpu_short",
        mem_mb=8000,
        time="00:15:00"
    threads: 1
    shell:
        """
        python -u workflow/scripts/06.2_concatenate_annotations.py \
            --civ {input.civ} \
            --okb {input.okb} \
            --cat maf \
            --output {output} &> {log}
        """


# Apply last filtering steps
rule somatic_maf_last_filtering:
    input:
        filters="%s/aggregate/somatic_maf/somatic_calls_filters_prefinal.tsv.gz" % R_FOLDER,
        maf="%s/aggregate/somatic_maf/somatic_calls_prefinal.maf.gz" % R_FOLDER,
        maf_ann="%s/aggregate/somatic_maf/somatic_calls_union_ann_prefinal.maf.gz" % R_FOLDER,
        maf_civ="%s/aggregate/somatic_maf/somatic_calls_civic_prefinal.maf.gz" % R_FOLDER,
        maf_okb="%s/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.maf.gz" % R_FOLDER
    output:
        filters="%s/aggregate/somatic_maf/somatic_calls_filters.tsv.gz" % R_FOLDER,
        maf="%s/aggregate/somatic_maf/somatic_calls.maf.gz" % R_FOLDER,
        maf_ann="%s/aggregate/somatic_maf/somatic_calls_union_ann.maf.gz" % R_FOLDER,
        maf_civ="%s/aggregate/somatic_maf/somatic_calls_civic.maf.gz" % R_FOLDER,
        maf_okb="%s/aggregate/somatic_maf/somatic_calls_oncokb.maf.gz" % R_FOLDER
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/somatic_maf_last_filtering/somatic_maf_last_filtering.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_maf_last_filtering/somatic_maf_last_filtering.log" % L_FOLDER
    resources:
        partition="cpu_short",
        mem_mb=24000,
        time="00:45:00"
    threads: 1
    shell:
        """
        python -u workflow/scripts/03.3_maf_last_filtering_steps.py \
            --inp_filters {input.filters} \
            --inp_maf {input.maf} \
            --inp_maf_ann {input.maf_ann} \
            --inp_maf_civ {input.maf_civ} \
            --inp_maf_okb {input.maf_okb} \
            --out_filters {output.filters} \
            --out_maf {output.maf} \
            --out_maf_ann {output.maf_ann} \
            --out_maf_civ {output.maf_civ} \
            --out_maf_okb {output.maf_okb} &> {log}
        """


# Aggregate tables with pileup data
rule somatic_maf_with_pileup:
    input:
        expand("%s/annotation/somatic_maf_add_pileup/{tsample_1}_vs_{nsample}_w_pileup_{tsample_2}.maf" % R_FOLDER,
            get_allowed_triplets_tumor_tumor_normal(), tsample_1=tsamples, nsample=nsamples_na, tsample_2=tsamples)
    output:
        "%s/aggregate/somatic_maf/somatic_calls_w_pileup.all.maf.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_maf_with_pileup/somatic_maf_with_pileup.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_maf_with_pileup/somatic_maf_with_pileup.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        python -u workflow/scripts/06.1_concatenate_tables.py \
            --files {input} \
            --keep_header \
            --output {output} &> {log}
        """

####
#### Somatic MSI ####
####

# Aggregate all somatic MSI tables.
rule somatic_msi_aggregate:
    input:
        expand("%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples)
    output:
        "%s/aggregate/somatic_msi/somatic_msi.tsv" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_msi/somatic_msi_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_msi/somatic_msi_aggregate.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=15
    shell:
        """
        python -u workflow/scripts/06.3_concatenate_msi.py \
            --input {input} \
            --output {output} &> {log}
        """

####
#### Somatic purity/ploidy ####
####

# Aggregate all somatic purity and ploidy values.
rule somatic_ppy_aggregate:
    input:
        expand("%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_ppy/somatic_ppy.tsv" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_ppy/somatic_ppy_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_ppy/somatic_ppy_aggregate.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=15
    shell:
        """
        python -u workflow/scripts/06.4_concatenate_ppy.py \
            --input {input} \
            --output {output} &> {log}
        """



####
#### Copy number variants ####
####

rule somatic_cna_filters_aggregate:
    input:
        expand("%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_cna/somatic_calls_filters.tsv.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_filters_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_filters_aggregate.log" % L_FOLDER
    params:
        threshold=config["params"]["cnv"]["calls_threshold"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        files=( {input} )
        if [[ ${{#files[@]}} == 1 ]]; then
            zcat ${{files[0]}} | gzip > {output} 2> {log}
        else
            {{ zcat ${{files[@]:0:1}}; zgrep --no-filename -v "^##\|Tumor_Sample_Barcode" ${{files[@]:1}}; }} | gzip > {output} 2> {log}
        fi
        """


rule somatic_cna_calls_aggregate:
    input:
        expand("%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_cna/somatic_calls.tsv.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_calls_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_calls_aggregate.log" % L_FOLDER
    params:
        threshold=config["params"]["cnv"]["calls_threshold"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        files=( {input} )
        if [[ ${{#files[@]}} == 1 ]]; then
            zcat ${{files[0]}} | grep "Tumor_Sample_Barcode\|PASS" | gzip > {output} 2> {log}
        else
            {{ zcat ${{files[@]:0:1}}; zgrep --no-filename -v "^##\|Tumor_Sample_Barcode" ${{files[@]:1}}; }} | grep "Tumor_Sample_Barcode\|PASS" | gzip > {output} 2> {log}
        fi
        """


rule somatic_cna_chr_arm_aggregate:
    input:
        expand("%s/calling/somatic_cnv_chr_arm/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_cna/somatic_calls_per_chr_arm.tsv.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_chr_arm_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_chr_arm_aggregate.log" % L_FOLDER
    threads: 2
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        cat {input} | sed -n '1p;/^Tumor/ !p' | gzip > {output} 2> {log}
        """


rule somatic_cna_sum_aggregate:
    input:
        expand("%s/calling/somatic_cnv_sum/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_cna/somatic_calls_summary_statistics.tsv.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_sum_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_sum_aggregate.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        cat {input} | sed -n '1p;/^Tumor/ !p' | gzip > {output} 2> {log}
        """


rule somatic_cna_table_aggregate:
    input:
        expand("%s/calling/somatic_cnv_table/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
            get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na)
    output:
        "%s/aggregate/somatic_cna/somatic_segments.tsv.gz" % R_FOLDER
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_table_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_table_aggregate.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        cat {input} | sed -n '1p;/^Tumor/ !p' | gzip > {output} 2> {log}
	"""


# Aggregate all somatic civic-annotated MAF tables.
rule somatic_cna_civic_aggregate:
    input:
        lambda w: get_input_concatenate(w, typ="cna", db="civic")
    output:
        "%s/aggregate/somatic_cna/somatic_calls_civic.tsv.gz" % R_FOLDER
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_civic_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_civic_aggregate.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        python -u workflow/scripts/06.1_concatenate_tables.py \
            --files {input} \
            --output {output} &> {log}
        """


rule somatic_cna_oncokb_aggregate:
    input:
        lambda w: get_input_concatenate(w, typ="cna", db="oncokb")
    output:
        "%s/aggregate/somatic_cna/somatic_calls_oncokb.tsv.gz" % R_FOLDER
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_oncokb_aggregate.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_oncokb_aggregate.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        python -u workflow/scripts/06.1_concatenate_tables.py \
            --files {input} \
            --output {output} &> {log}
        """


# Aggregate oncokb and civic mutation annotations.
rule somatic_cna_union_ann:
    input:
        civ="%s/aggregate/somatic_cna/somatic_calls_civic.tsv.gz" % R_FOLDER,
        okb="%s/aggregate/somatic_cna/somatic_calls_oncokb.tsv.gz" % R_FOLDER
    output:
        "%s/aggregate/somatic_cna/somatic_calls_union_ann.tsv.gz" % R_FOLDER
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/somatic_cna/somatic_cna_union_ann.tsv" % B_FOLDER
    log:
        "%s/aggregate/somatic_cna/somatic_cna_union_ann.log" % L_FOLDER
    resources:
        partition="cpu_short",
        mem_mb=8000,
        time="00:15:00"
    threads: 1
    shell:
        """
        python -u workflow/scripts/06.2_concatenate_annotations.py \
            --civ {input.civ} \
            --okb {input.okb} \
            --cat cna \
            --output {output} &> {log}
        """

# Aggregate all WES qc metrics
# This rule should be run after the pipeline is finished using the command
# snakemake -s workflow/Snakefile --profile [your_profile] results/aggregate/wes_analyses_summary.xlsx
rule summary_wes_qc:
    input:
        tnp = config["general"]["tumor_normal_pairs"],
    output:
        "%s/aggregate/wes_analyses_summary.xlsx" % R_FOLDER
    conda:
        "../envs/python.yaml"
    benchmark:
        "%s/aggregate/summary_wes_qc/summary_wes_qc.tsv" % B_FOLDER
    log:
        "%s/aggregate/summary_wes_qc/summary_wes_qc.log" % L_FOLDER
    resources:
        partition="cpu_short",
        mem_mb=8000,
        time="00:60:00"
    params:
        somatic_facets = get_files_recursively("%s/calling/somatic_cnv_facets" % R_FOLDER, ".vcf.gz"),
        somatic_maf_agg = "%s/aggregate/somatic_maf/somatic_calls.maf.gz" % R_FOLDER,
        qc_fastqc_r1 = get_files_recursively("%s/qc/fastqc_1" % R_FOLDER, "_R1_fastqc.zip"),
        qc_fastqc_r2 = get_files_recursively("%s/qc/fastqc_1" % R_FOLDER, "_R2_fastqc.zip"),
        qc_ncm = get_files_recursively("%s/qc/ngs_checkmate" % R_FOLDER, "_all.txt"),
        qc_flagstat = get_files_recursively("%s/qc/samtools_flagstat" % R_FOLDER, "_flagstat.tsv"),
        qc_hsmetrics = get_files_recursively("%s/qc/collect_hs_metrics" % R_FOLDER, "_hs_metrics.tsv"),
        min_cov_t = 40,
        min_cov_n = 40,
        min_10x_pct_t = 60,
        min_10x_pct_n = 60
    threads: 1
    shell:
        """
        python -u workflow/scripts/06.6_concatenate_qc.py \
            --tnp {input.tnp} \
            --somatic_facets {params.somatic_facets} \
            --somatic_maf_agg {params.somatic_maf_agg} \
            --qc_fastqc_r1 {params.qc_fastqc_r1} \
            --qc_fastqc_r2 {params.qc_fastqc_r2} \
            --qc_ncm {params.qc_ncm} \
            --qc_flagstat {params.qc_flagstat} \
            --qc_hsmetrics {params.qc_hsmetrics} \
            --min_cov_t {params.min_cov_t} \
            --min_cov_n {params.min_cov_n} \
            --min_10x_pct_t {params.min_10x_pct_t} \
            --min_10x_pct_n {params.min_10x_pct_n} \
            --output {output} &> {log}
        """
