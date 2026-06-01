# Hellbender HPC notes

University of Missouri's HPC cluster. Used for long-running MCMC.

## Quick facts

- **Cluster:** Hellbender (University of Missouri)
- **User:** mjd4d
- **Email for SBATCH notifications:** mjd4d@umsystem.edu
- **Default partition:** `general`
- **R modules:** `r/4.4.0`, `r/4.5`
- **Wall time cap (typical):** up to 48 hours
- **Project files:** `~/LS_spline/`

## Submitting a job

Standard SBATCH header:

```bash
#!/bin/bash
#SBATCH --partition=general
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjd4d@umsystem.edu
#SBATCH --job-name=lsspline_run
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load r/4.4.0
cd ~/LS_spline
Rscript run_interaction_2x2_v2.R
```

Submit:

```bash
sbatch my_job.sub
squeue -u mjd4d   # check status
scancel <jobid>   # cancel
```

## DOS line-ending gotcha

If you edit `.sub` scripts on Windows, they may have CRLF line endings, which break SBATCH parsing on Linux. Fix:

```bash
sed -i 's/\r$//' my_job.sub
```

Do this **every time** before submitting a script edited on Windows. It's a frequent source of "why won't this submit?" frustration.

## R module choice

- `r/4.4.0` is the safe default (matches what's been used throughout)
- `r/4.5` matches local Windows R 4.5 — useful if you want exact reproducibility with local development

If a package compiles fine locally but fails on Hellbender, try the other R module.

## File transfer

Edit locally on Windows, push to Hellbender via your usual transfer tool (scp / WinSCP / rsync).

## Connects to

- [[R file index]] (what to actually run)
- [[Bug log]] (gotchas including CRLF)
