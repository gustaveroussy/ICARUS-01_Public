from snakemake.utils import min_version

min_version("5.4.0")

configfile: "config/config.yaml"

L_FOLDER = "logs"
B_FOLDER = "benchmarks"
R_FOLDER = "../../results/wes_analysis"
D_FOLDER = "../../data"

oncoplots = []
for plot in config["oncoplots"]:
    oncoplots.append(plot['name'])

oncoplots_paired = []
for plot in config["oncoplots_paired"]:
    oncoplots_paired.append(plot['name'])
