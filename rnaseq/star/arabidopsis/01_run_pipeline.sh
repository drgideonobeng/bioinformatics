#!/usr/bin/env bash

#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
set -euo pipefail

echo "########################################################"
echo " INITIALIZING AUTOMATED RNA-SEQ UPSTREAM PIPELINE"
echo "########################################################"

# Run scripts sequentially. The 'bash' command executes them.
bash scripts/01_build_refs.sh
bash scripts/02_download_sra.sh
bash scripts/03_fastqc_raw.sh
bash scripts/04_trim_fastp.sh
bash scripts/05_quantify_salmon.sh

echo "########################################################"
echo " PIPELINE COMPLETE!"
echo " All samples have been downloaded, cleaned, and quantified."
echo " Proceed to R for Differential Expression (DESeq2)."
echo "########################################################"
