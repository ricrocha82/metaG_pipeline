#!/bin/bash
# make_samplesheet2.sh
# Builds the samplesheet for pipeline2-binning-mags:
#   sample,assembly,short_reads_1,short_reads_2
#
# Usage:
#   bash make_samplesheet2.sh <assembly_dir> <reads_dir> [outfile]
#
# Example:
#   bash make_samplesheet2.sh \
#       /fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/assembly/06_megahit \
#       /fs/project/PAS1117/ricardo/pig_burn/host_depletion/raw_reads/metaG \
#       assets/samplesheet2.csv

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <assembly_dir> <reads_dir> [outfile]" >&2
    echo "  assembly_dir: folder containing one subfolder per sample, each with <sample>.contigs.fa" >&2
    echo "  reads_dir:    folder containing <sample>_S*_R1_001.fastq.gz / _R2_001.fastq.gz" >&2
    echo "  outfile:      where to write the CSV (default: assets/samplesheet2.csv)" >&2
    exit 1
fi

ASSEMBLY_DIR="$1"
READS_DIR="$2"
OUTFILE="${3:-assets/samplesheet2.csv}"

if [[ ! -d "$ASSEMBLY_DIR" ]]; then
    echo "ERROR: assembly_dir does not exist: $ASSEMBLY_DIR" >&2
    exit 1
fi
if [[ ! -d "$READS_DIR" ]]; then
    echo "ERROR: reads_dir does not exist: $READS_DIR" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUTFILE")"
echo "sample,assembly,short_reads_1,short_reads_2" > "$OUTFILE"

# Loop over each sample folder under the assembly directory.
# Assumes one contigs.fa per sample, named <sample>.contigs.fa,
# and paired reads named <sample>_S*_R1_001.fastq.gz / _R2_001.fastq.gz
for sample_dir in "$ASSEMBLY_DIR"/*/; do
    sample=$(basename "$sample_dir")
    assembly="${sample_dir}${sample}.contigs.fa"

    # Find matching R1/R2 files by sample prefix (handles the _S### suffix)
    r1=$(find "$READS_DIR" -maxdepth 1 -name "${sample}_S*_R1_001.fastq.gz" | head -n1)
    r2=$(find "$READS_DIR" -maxdepth 1 -name "${sample}_S*_R2_001.fastq.gz" | head -n1)

    # Sanity checks -- fail loudly rather than writing a bad row
    if [[ ! -f "$assembly" ]]; then
        echo "WARNING: assembly not found for sample '$sample': $assembly" >&2
        continue
    fi
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "WARNING: reads not found for sample '$sample' in $READS_DIR" >&2
        continue
    fi

    echo "${sample},${assembly},${r1},${r2}" >> "$OUTFILE"
done

echo "Wrote $OUTFILE:" >&2
cat "$OUTFILE" >&2


# bash /fs/project/PAS1117/ricardo/pipeline2-binning-mags/bin/make_samplesheet2.sh \
#     /fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/assembly/06_megahit \
#     /fs/project/PAS1117/ricardo/pig_burn/host_depletion/raw_reads/metaG \
#     /fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/samplesheet2.csv