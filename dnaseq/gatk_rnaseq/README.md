# GATK RNA-seq Variant Discovery (Practice Project)

Beginner-friendly, step-by-step practice repo for calling germline SNPs/indels from **RNA-seq** using GATK.

## Scope (important)
- Follows the official GATK RNA-seq short variant discovery workflow.
- Intended to be run **per sample** (not joint calling for RNA-seq in this workflow).

## Project layout
- `env/` environment files
- `scripts/` runnable scripts
- `data/` inputs (kept local; not committed)
- `results/` outputs (kept local; not committed)
- `docs/` written guide

## Steps (we will do these one-at-a-time)
1. Repo + environment setup (current step)
2. Get a small practice RNA-seq dataset
3. Reference + known-sites resources
4. Pre-processing: duplicates + SplitNCigarReads (RNA-seq specific)
5. Call variants with HaplotypeCaller
6. Basic filtering + QC
7. (Optional) functional annotation
