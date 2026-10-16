# Find the exact 2.6.1 tag
curl -s "https://quay.io/api/v1/repository/biocontainers/gtdbtk/tag/?filter_tag_name=like:2.6.1" \
  | grep -o '"name": *"[^"]*"'

cd /fs/ess/PAS0822/modules/singularity/
apptainer pull gtdbtk_2.6.1.sif docker://quay.io/biocontainers/gtdbtk:<tag-from-above>
apptainer exec gtdbtk_2.6.1.sif sh -c 'command -v ps; gtdbtk --version'