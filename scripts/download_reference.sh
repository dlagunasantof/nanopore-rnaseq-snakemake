#!/usr/bin/env bash
# Download the GENCODE human transcriptome and annotation into resources/.
# Usage: bash scripts/download_reference.sh [release]   (default: 47)
# Check the latest release at https://www.gencodegenes.org/human/ and keep
# config/config.yaml in sync with the release you download.
set -euo pipefail

RELEASE="${1:-47}"
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${RELEASE}"
mkdir -p resources

for f in "gencode.v${RELEASE}.transcripts.fa.gz" "gencode.v${RELEASE}.annotation.gtf.gz"; do
    if [ ! -s "resources/${f}" ]; then
        echo "Downloading ${f}"
        wget -q --show-progress -O "resources/${f}" "${BASE}/${f}"
    fi
done
echo "Done. Files in resources/:"
ls -lh resources/
