# Download databases required for annotating with annovar
# See https://annovar.openbioinformatics.org/en/latest/user-guide/filter
rule annovar_downdb:
    output:
        "%s/humandb/hg19_{resource}.txt" % config["params"]["annovar"]["path"],
        "%s/humandb/hg19_{resource}.txt.idx" % config["params"]["annovar"]["path"]
    benchmark:
        "%s/annotation/annovar_downdb/{resource}.tsv" % B_FOLDER
    log:
        "%s/annotation/annovar_downdb/{resource}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        annovar_dir=config["params"]["annovar"]["path"],
        buildver=config["params"]["annovar"]["buildver"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=90
    shell:
        """{params.annovar_dir}/annotate_variation.pl -buildver {params.buildver} \
            -downdb \
            -webfrom annovar {wildcards.resource} {params.annovar_dir}/humandb 2> {log}"""


# Install VEP code
# See https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html
# See http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# See https://githubmemory.com/repo/Ensembl/ensembl-vep/issues/930
rule vep_install:
    output:
        "%s/INSTALL.pl" % config["params"]["vep"]["path"],
        "%s/%s.gz" % (config["params"]["vep"]["cache"], config["params"]["vep"]["fasta"])
    benchmark:
        "%s/annotation/vep_install/vep.tsv" % B_FOLDER
    log:
        "%s/annotation/vep_install/vep.log" % L_FOLDER
    params:
        assembly=config["ref"]["build"],
        species=config["ref"]["species"],
        plugins=",".join(config["params"]["vep"]["plugins"]),
        cache=config["params"]["vep"]["cache"],
        release=config["params"]["vep"]["release"],
        path=config["params"]["vep"]["path"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=240
    shell:
        """
        bash workflow/scripts/02.1_vep_install.sh \
            -a {params.assembly} \
            -s {params.species} \
            -g {params.plugins} \
            -c {params.cache} \
            -r {params.release} \
            -p {params.path} &> {log}
        """


# Download plugins databases required for annotating with vep plugins
# See https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
rule vep_plugins:
    input:
        "%s/INSTALL.pl" % config["params"]["vep"]["path"]
    output:
        "%s/Plugins/dbNSFP4.1a_%s.gz" % (config["params"]["vep"]["cache"],  config["ref"]["build"].lower()),
        "%s/Plugins/LoFtool_scores.txt" % config["params"]["vep"]["cache"],
        "%s/Plugins/Condel/config/condel_SP.conf" % config["params"]["vep"]["cache"]
    benchmark:
        "%s/annotation/vep_plugins/vep.tsv" % B_FOLDER
    log:
        "%s/annotation/vep_plugins/vep.log" % L_FOLDER
    params:
        cache=config["params"]["vep"]["cache"],
        assembly=config["ref"]["build"],
        dbnsfp_version="4.1"
    threads: 6
    resources:
        queue="longq",
        mem_mb=20000,
        time_min=1440
    shell:
        """
        bash workflow/scripts/02.2_vep_plugins.sh \
            -c {params.cache} \
            -a {params.assembly} \
            -v {params.dbnsfp_version} &> {log}
        """


# Download oncokb code and required databases if not done yet.
rule oncokb_install:
    output:
        "%s/CnaAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/FusionAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/MafAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/cancerGeneList_oncokb_annotated.tsv" % config["params"]["oncokb"]["data_dir"]
    benchmark:
        "%s/annotation/oncokb_install/oncokb.tsv" % B_FOLDER
    log:
        "%s/annotation/oncokb_install/oncokb.log" % L_FOLDER
    params:
        code_dir=config["params"]["oncokb"]["code_dir"],
        data_dir=config["params"]["oncokb"]["data_dir"],
        gene_list=config["params"]["oncokb"]["gene_list"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=15
    shell:
        """
        bash workflow/scripts/04.2_oncokb_install.sh \
            -c {params.code_dir} \
            -d {params.data_dir} \
            -g {params.gene_list} &> {log}
        """
