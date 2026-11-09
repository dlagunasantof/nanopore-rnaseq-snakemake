#!/usr/bin/env bash
# Find and download SG-NEx nanopore FASTQ files from the public AWS bucket (no AWS account needed).
# Requires the AWS CLI:  conda install -c conda-forge awscli
#
# Step 1 - list the FASTQ files available for a cell line and protocol:
#   bash scripts/download_sgnex.sh list HepG2 cDNA
#
# Step 2 - download one file and save it with a short name:
#   bash scripts/download_sgnex.sh get <s3://path/from/step1> data/HepG2_rep1.fastq.gz
#
# Protocols in SG-NEx: cDNA (PCR-cDNA; also matches cDNAStranded), directcDNA, directRNA.
# File names look like SGNex_K562_cDNAStranded_replicate3_run3.fastq.gz
# Data description and citation: https://registry.opendata.aws/sgnex/
set -euo pipefail

BUCKET="s3://sg-nex-data"
cmd="${1:-}"

case "$cmd" in
  list)
    cell="${2:?cell line, e.g. HepG2}"
    protocol="${3:?protocol, e.g. cDNA}"
    aws s3 ls --no-sign-request --recursive "${BUCKET}/" \
      | grep -iE "\.f(ast)?q(\.gz)?$" \
      | grep -i "_${cell}_" \
      | grep -i "_${protocol}" \
      | awk -v b="${BUCKET}" '{printf "%8.2f GB  %s/%s\n", $3/1e9, b, $4}'
    ;;
  get)
    src="${2:?S3 path from the list command}"
    dest="${3:?local file name, e.g. data/HepG2_rep1.fastq.gz}"
    mkdir -p "$(dirname "$dest")"
    aws s3 cp --no-sign-request "$src" "$dest"
    ;;
  *)
    sed -n '2,13p' "$0"
    exit 1
    ;;
esac
