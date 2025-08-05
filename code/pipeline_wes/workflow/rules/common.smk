from itertools import product
import os
import pandas as pd
import re
from snakemake.utils import validate
from snakemake.utils import min_version
import subprocess

min_version("5.4.0")

B_FOLDER = "workflow/benchmarks"
L_FOLDER = "workflow/logs"
R_FOLDER = "results"

###### Config file and sample sheets #####
configfile: "config/config.yaml"

table = pd.read_table(config["general"]["samples"], dtype=str).set_index(["Sample_Id"], drop=False)
if "IRODS_Status" in table:
    mask_ok = table["IRODS_Status"]=="OK"
else:
    mask_ok = pd.Series(True, index=table.index)

samples = table.loc[mask_ok]["Sample_Id"].tolist()
nsamples = table.loc[(mask_ok)&(table["Sample_Type"]=="DNA_N")]["Sample_Id"].tolist()
nsamples_na = nsamples + ["NA"]
tsamples = table.loc[(mask_ok)&(table["Sample_Type"]=="DNA_T")]["Sample_Id"].tolist()
fastqs = ["%s_R%d" % (sample, i) for sample in samples for i in [1,2]]

# Define constraints
wildcard_constraints:
    nsample = "|".join([re.escape(x) for x in nsamples_na]),
    tsample = "|".join([re.escape(x) for x in tsamples+["^$"]]),
    sample = "|".join([re.escape(x) for x in samples+["^$"]]),

if config["general"]["regenerate_bam"]:
    samples_redo_bam = samples
    samples_copy_bam = []
else:
    samples_existing_bam = [x for x in samples if os.path.isfile("results/mapping/%s.nodup.recal.bam" % x)]

    if config["general"]["archive_folder"] is not None:
        archive_bam_folder = os.path.join(config["general"]["archive_folder"], "results/mapping")
        archive_bams = [x for x in os.listdir(archive_bam_folder) if x.endswith(".bam")]
        samples_copy_bam = [re.sub(".nodup.recal.bam", "", x) for x in archive_bams]
        samples_copy_bam = [x for x in samples_copy_bam if x in samples]
    else:
        samples_copy_bam = []
    samples_existing_bam += samples_copy_bam
    samples_redo_bam = [x for x in samples if x not in samples_existing_bam]

fastqs_redo_bam = ["%s_R%d" % (sample, i) for sample in samples_redo_bam for i in [1,2]]

if len(samples_copy_bam)>0:
    print("-INFO: will copy the BAMs from archive for the following samples")
    print("\t" + "\n\t".join(samples_copy_bam))


if "normal_pool" in config["general"] and config["general"]["normal_pool"] is not None:
    table_normal_pool = pd.read_table(config["general"]["normal_pool"], dtype=str)
    nsamples_normal_pool = table_normal_pool["Sample_Id"].tolist()
    nsample_normal_pool = "_".join(table_normal_pool["Sample_Id"].tolist())
else:
    nsamples_normal_pool = []
    nsample_normal_pool = ""


if "bam_coverage" in config["general"] and config["general"]["bam_coverage"] is not None:
    with open(config["general"]["bam_coverage"], "r") as file:
        bam_coverage_regions = [line.strip() for line in file.readlines()]
else:
    bam_coverage_regions = None


def get_intervals_for_pon_target_file(target):
    folder = config["target_files"]["intervals_pon"][target]
    file = "%s_{interval}.bed" % config["target_files"]["target"][target]
    intervals, = glob_wildcards(os.path.join(folder, file))
    return intervals


if "Capture_Kit" in table:
    targets = table["Capture_Kit"].unique().tolist()
    intervals = set().union(*[get_intervals_for_pon_target_file(target) for target in targets])
else:
    targets = []
    intervals = []

if "download_gdc" in config.keys():
    ids_gdc_somatic_mc3 = pd.read_table(config["download_gdc"]["manifests"]["somatic_mc3"])["id"].tolist()


##### Helper functions #####

def filter_combinator(combinator, comblist, white_list=True):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
        # Use frozenset instead of tuple
        # in order to accomodate
        # unpredictable wildcard order
            if white_list:
                if frozenset(wc_comb) in comblist:
                    yield wc_comb
            else:
                if frozenset(wc_comb) not in comblist:
                    yield wc_comb
    return filtered_combinator


def get_allowed_pairs_tumor_normal():
    allowed = []
    df_pairs = pd.read_table(config["general"]["tumor_normal_pairs"]).fillna("NA")
    for (tsample, nsample) in zip(df_pairs["DNA_T"], df_pairs["DNA_N"]):
        allowed.append(frozenset({("tsample", tsample), ("nsample", nsample)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_pairs_tumor_normal():
    allowed = []
    df_pairs = pd.read_table(config["general"]["tumor_normal_pairs"]).fillna("NA")
    for (tsample, nsample) in zip(df_pairs["DNA_T"], df_pairs["DNA_N"]):
        allowed.append(frozenset({("tsample", tsample), ("nsample", nsample)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_triplets_tumor_tumor_normal():
    allowed = []
    df_triplets = pd.read_table(config["general"]["tumor_tumor_pairs"]).fillna("NA")
    for (tsample_1, tsample_2, nsample) in zip(df_triplets["DNA_T1"], df_triplets["DNA_T2"], df_triplets["DNA_N"]):
        allowed.append(frozenset({("tsample_1", tsample_1), ("tsample_2", tsample_2), ("nsample", nsample)}))
        allowed.append(frozenset({("tsample_1", tsample_2), ("tsample_2", tsample_1), ("nsample", nsample)}))
    return filter_combinator(product, allowed, white_list=True)


def get_fastqs_irods(wildcards):
    """Get fastq files irods paths of given sample."""
    fastqs = table.loc[wildcards.sample, ["FASTQ_1", "FASTQ_2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.FASTQ_1, "r2": fastqs.FASTQ_2}
    return {"r1": fastqs.FASTQ_1}


def get_fastqs_local(wildcards):
    """Get fastq files local paths of given sample."""
    return {"r%d" % i: "%s/data/fastq/%s_R%d.fastq.gz" % (R_FOLDER, wildcards.sample, i) for i in [1,2]}


def get_fastqs_names(wildcards):
    """Get fastq files base name of given sample."""
    fastqs = table.loc[wildcards.sample, ["FASTQ_1_Name", "FASTQ_2_Name"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.FASTQ_1_Name, "r2": fastqs.FASTQ_2_Name}
    return {"r1": fastqs.FASTQ_1_Name}



def get_column_table_sample(wildcards, col):
    """Get the value of the column col for the sample"""
    try:
        value = table.loc[wildcards.sample, col]
    except AttributeError:
        try:
            value = table.loc[wildcards.tsample, col]
        except AttributeError:
            try:
                value = table.loc[wildcards.nsample, col]
            except AttributeError:
                if wildcards.sample_pair=="all_samples":
                    value = ""
                else:
                    tsample = wildcards.sample_pair.split("_vs_")[0]
                    value = table.loc[tsample, col]
    return value


def get_target_name_sample(wildcards):
    """Get name of target file used for the sample"""
    return get_column_table_sample(wildcards, "Capture_Kit")


def get_target_file_sample(wildcards, file="bed_padded"):
    """Get path to target files used for the sample"""
    capture_kit = get_target_name_sample(wildcards)
    return config["target_files"][file][capture_kit]


def get_target_name_tumor_normal(wildcards):
    """Get name of target file used for the pair tumor, normal"""
    try:
        n_capture_kit = table.loc[wildcards.nsample, "Capture_Kit"]
        t_capture_kit = table.loc[wildcards.tsample, "Capture_Kit"]
        if n_capture_kit!=t_capture_kit:
            captur_kit = "all_targets_intersect"
        else:
            capture_kit = t_capture_kit
    except AttributeError:
        t_capture_kit = table.loc[wildcards.tsample, "Capture_Kit"]
        capture_kit = t_capture_kit
    return capture_kit


def get_target_file_tumor_normal(wildcards, file="bed_padded"):
    """Get path to target file to be used when calling somatic mutations.
    It may happen that the targets used for the normal and the tumor are different.
    In this case, the intersection of all target files of the project is used
    as the target file.
    """
    capture_kit = get_target_name_tumor_normal(wildcards)
    return config["target_files"][file][capture_kit]


def get_samples_target_file(wildcards, type="DNA_T"):
    """Return the list of samples prepared with a given target file."""
    samples_sel = []
    for sample, sample_type, sample_target in zip(table["Sample_Id"], table["Sample_Type"], table["Capture_Kit"]):
        if sample_type==type and sample_target==wildcards.target:
            samples_sel.append(sample)
    return samples_sel


def get_input_vcfs_pon_genomicsdb(wildcards):
    nsamples_target = get_samples_target_file(wildcards, type="DNA_N")
    return expand("%s/calling/pon_mutect2/{nsample}.vcf.gz" % R_FOLDER, nsample=nsamples_target)


def get_input_bed_pon_genomicsdb(wildcards):
    folder = config["target_files"]["intervals_pon"][wildcards.target]
    file = "%s_%s.bed" % (config["target_files"]["target"][wildcards.target], wildcards.interval)
    return os.path.join(folder, file)


def get_annovar_operation():
    protocols=config["params"]["annovar"]["protocol"]
    operations=[]
    for protocol in protocols:
        if protocol.startswith("refGene"):
            operations.append("g")
        else:
            operations.append("f")
    return ",".join(operations)


def get_target_names():
	return sorted(list(config["target_files"]["target"].keys()))


def get_target_files(names, file="bed"):
    return [config["target_files"][file][name] for name in names]


def get_tumor_type_mskcc_oncotree(wildcards):
    """Get the tumor type MSKCC oncotree of the sample"""
    return get_column_table_sample(wildcards, "MSKCC_Oncotree")


def get_tumor_type_civic(wildcards):
    """Get the tumor type Civic_Disease of the sample"""
    return get_column_table_sample(wildcards, "Civic_Disease")


def get_facets_diplogr(wildcards):
    """Get manually chosen value of dipLogR for FACETS, if any."""
    try:
        df_diplogr = pd.read_table(config["params"]["cnv"]["cnv_facets"]["diplogr"]).set_index(["DNA_P"])
        try:
            dna_p = "%s_vs_%s" % (wildcards.tsample, wildcards.nsample)
        except:
            dna_p = "%s_vs_NA" % wildcards.tsample
        diplogr = df_diplogr.loc[dna_p, "DipLogR"]
        return diplogr
    except:
        return -999


def get_input_vep_vcf(w):
    if config["params"]["vep"]["annotate_all"]:
        return "%s/calling/somatic_maf_filters/%s_vs_%s.vcf.gz" % (R_FOLDER, w.tsample, w.nsample)
    else:
        return "%s/calling/somatic_maf_select_pass/%s_vs_%s.vcf.gz" % (R_FOLDER, w.tsample, w.nsample)


def get_input_concatenate(w, typ, db):
    input_folder = "%s/annotation/somatic_%s_%s" % (R_FOLDER, typ, db)

    if config["params"][db]["run_per_sample"][typ]:
        sample_pairs = expand("{tsample}_vs_{nsample}", get_allowed_pairs_tumor_normal(),
            tsample=tsamples, nsample=nsamples_na)
    else:
        sample_pairs = ["all_samples"]

    if typ=="maf":
        return ["%s/%s.maf" % (input_folder, sample_pair) for sample_pair in sample_pairs]
    elif typ=="cna":
        return ["%s/%s.tsv" % (input_folder, sample_pair) for sample_pair in sample_pairs]


def get_files_recursively(folder, suffix):
    files = []
    for root, subdirs, filenames in os.walk(folder):
            for filename in filenames:
                if filename.endswith(suffix):
                    filepath = os.path.join(root, filename)
                    files.append(filepath)
    return files


def get_coordinates_gene(w):
    cmd="zcat %s " % config["ref"]["gencode"] + \
        "| awk -F '\t' 'NR==1; {if($(3)==\"gene\") print $0}' " + \
        "| grep \"gene_name=%s;\" " % w.gene + \
        "| awk -F '\t' '{printf (\"%5s,%s,%s\", $1, $4, $5)}'"
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = ps.communicate()[0].decode('utf-8')
    chr, start, end = out.split(",")
    chr = re.sub("^chr", "", chr)

    return (chr,start,end)
