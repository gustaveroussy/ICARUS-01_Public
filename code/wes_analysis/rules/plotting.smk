def build_oncoplot_args(params_dict):
    args = []
    for key, value in params_dict.items():
        # Handle lists for arguments like --gene_list_names
        if isinstance(value, list):
            args.append(f"--{key} {' '.join(map(str, value))}")
        else:
            args.append(f"--{key} {value}")
    return ' '.join(args)


def get_oncoplot_args_from_wildcards(wildcards, config):
    """
    Finds the plot config based on wildcards and returns the argument string.
    """
    plot_name = wildcards.plot_name

    plot_config = next(
        (p for p in config["oncoplots"] if p["name"] == plot_name), None
    )

    if not plot_config:
        # This raises a clear, helpful error if the name isn't found
        raise snakemake.exceptions.WorkflowError(
            f"FATAL: Could not find a plot configuration for name '{plot_name}' in the 'oncoplots' section")

    # Get the nested params dictionary and build the command-line string
    params_dict = plot_config["params"]
    return build_oncoplot_args(params_dict)


rule draw_oncoplot_sample:
    wildcard_constraints:
        plot_name = "|".join(oncoplots)
    output:
        pdf = f"{R_FOLDER}/oncoplot/{{plot_name}}.pdf",
        tab = f"{R_FOLDER}/oncoplot/{{plot_name}}.xlsx",
    input:
        script = "scripts/02.1_draw_oncoplot_samples.py",
        cln = config["data"]["inputs"]["cln_curated"],
        alt = config["data"]["outputs"]["alterations"],
        mut = config["data"]["inputs"]["mut_pass"],
        cna = config["data"]["inputs"]["cna_pass"],
        sam = config["data"]["inputs"]["sam"],
        ids_anonym = config["data"]["inputs"]["ids_anonym"]
    params:
        args_str = lambda wildcards: get_oncoplot_args_from_wildcards(wildcards, config)
    benchmark:
        f"{B_FOLDER}/draw_oncoplot/{{plot_name}}.tsv"
    log:
        f"{L_FOLDER}/draw_oncoplot/{{plot_name}}.log"
    conda:
        config["setup"]["py-icarus"]
    threads: 1
    shell:
        """
        python -u {input.script} \
            --cln {input.cln} \
            --alt {input.alt} \
            --mut {input.mut} \
            --cna {input.cna} \
            --sam {input.sam} \
            --ids_anonym {input.ids_anonym} \
            {params.args_str} \
            --output_pdf {output.pdf} \
            --output_tab {output.tab} &> {log}
        """
