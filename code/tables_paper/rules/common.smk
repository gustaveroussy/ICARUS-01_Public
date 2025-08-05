from snakemake.utils import min_version

min_version("5.4.0")

configfile: "config/config.yaml"

L_FOLDER = "logs"
B_FOLDER = "benchmarks"
R_FOLDER = "../../results/tables_paper"
D_FOLDER = "../../data"
