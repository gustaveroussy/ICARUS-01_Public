# ENTRY FILE: .bed file

#### how to generate .intervals file ===================================================================================

# .intervals files are Picard-style and should actually have .interval_list extension instead of .intervals
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists

# STEP 1: create conda environment with GATK installed
# for instance this environment
# name: gatk
# channels:
#     - bioconda
#     - conda-forge
#     - defaults
# dependencies:
#     - gatk4


# STEP 2: from your environment run
gatk --java-options '-Xmx2g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp' BedToIntervalList \
  I=target_files/SureSelect_CR2.bed O=target_files_recreate/SureSelect_CR2.intervals \
  SD=/mnt/beegfs/database/bioinfo/metaprism/wes/resources/b37_gatk/human_g1k_v37.dict

# for some reason that escape my understanding, the above command sometimes splits
# an interval into two and skips a base.
# for instance, interval
# 1   8384366   8384906  +  .
# is split into
# 1   8384366   8384786  +  .
# 1   8384787   8384906  +  .
# here the base 8384786 is consequently skipped
#

#### how to generate padded file =======================================================================================


# STEP 1: load commands in your PATH
# We will need samtools and bedtools command. Both are pre-installed on GR cluster, load them via
module load samtools
module load bedtools

# STEP 2: create a genome file from your reference FASTA via
cut -f1,2 /mnt/beegfs/database/bioinfo/metaprism/wes/resources/b37_gatk/human_g1k_v37.fasta.fai > human_g1k_v37.sizes
# If the index does not exist already, you may run
# samtools faidx /mnt/beegfs/database/bioinfo/metaprism/wes/resources/b37_gatk/human_g1k_v37.fasta


# STEP 3: pad your bed via
bedtools slop \
  -i target_files/SureSelect_CR2.bed \
  -g human_g1k_v37.sizes \
  -b 10 > target_files_recreate/SureSelect_CR2_padded_10n.bed


#### how to generate split files for PON ===============================================================================

# The PON is created via a multi-step process starting with GenomicsDBImport command on VCFs produced by running Mutect2
# in tumor-only mode on the normal samples. The GenomicsDBImport expects a list of genomic intervals to operate on (-L
# or --intervals argument). However, it is notably slow when the number of intervals exceeds a couple hundreds.
#
# See https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines
# for more info about the usage of GenomicsDBImport command.
#
# To split a file into files with identical number of lines, you may the split GNU command

mkdir target_files_recreate/SureSelect_CR2_for_pon
cd target_files_recreate/SureSelect_CR2_for_pon
split  --lines 500 \ # split into files of 500 lines
  --additional-suffix .bed \ # self-explanatory
  -d \ # digit suffixes
  -a 3 \ # suffixes over 3 characters i.e from 000 to 999
  ../../target_files/SureSelect_CR2.bed SureSelect_CR2_ # path to bed file to be split and prefix
