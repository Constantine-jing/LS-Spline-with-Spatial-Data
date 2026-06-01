#!/bin/bash
# Submit all 4 n=1500 jobs in parallel. Run from project root.
set -e
sbatch scripts/submit/run_n1500_headline_seed1.sub
sbatch scripts/submit/run_n1500_seed2.sub
sbatch scripts/submit/run_n1500_seed3.sub
sbatch scripts/submit/run_n1500_seed4.sub
echo "All 4 jobs submitted. Check queue with: squeue -u \$USER"
