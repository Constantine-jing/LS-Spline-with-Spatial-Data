# Hellbender Run — n=1500 All-Pairs

Step-by-step for running `scripts/run_n1500_allpairs.R` on Hellbender.
Do these steps manually after VPN into the cluster.

---

## 0. Before you leave your local machine

Make sure these files are ready locally:

```
output/sim_medsat_london_asthma_n1500_seed1_v2.rds   # n=1500 sim
output/sim_medsat_london_asthma_n500_seed1_v2.rds    # n=500 sim (for smoke test)
scripts/hb_smoke_test.R
scripts/run_n1500_allpairs.R
R/                                                    # all R/ source files
```

---

## 1. Transfer files to Hellbender

**File transfer tool: TBD** (fill in your preferred method — scp, rsync, FileZilla, etc.)

Suggested layout on Hellbender:

```
$HOME/medsat-pilot/
  R/
  scripts/
  output/
  logs/          # create this if it doesn't exist: mkdir -p logs
```

Example with rsync (adjust hostname and remote path):

```bash
# From your local machine:
rsync -avz \
  R/ scripts/ output/sim_medsat_london_asthma_n*.rds \
  <HELLBENDER_USER>@<HELLBENDER_HOST>:~/medsat-pilot/
```

---

## 2. Log into Hellbender

```bash
ssh <HELLBENDER_USER>@<HELLBENDER_HOST>
```

Navigate to the project root:

```bash
cd ~/medsat-pilot
```

---

## 3. Load R module

**R module command: TBD** — check available modules with:

```bash
module avail R
```

Then load (example, confirm the exact name):

```bash
module load R/<VERSION>
```

Verify:

```bash
Rscript --version
```

---

## 4. Set up personal R library (if needed)

**Personal R library path: TBD** — if the cluster doesn't have all packages
in the system library, install to a personal location:

```bash
mkdir -p ~/R/library
```

Then in R (or at the top of a script), set:

```r
.libPaths(c("~/R/library", .libPaths()))
```

Install required packages if missing (do this once in an interactive session,
not inside the production run):

```r
install.packages(c("Rcpp", "RcppEigen", "coda"), lib = "~/R/library")
```

---

## 5. Run the smoke test

This verifies the Hellbender R environment on the n=500 sim before
committing to the multi-hour production run.

```bash
cd ~/medsat-pilot
Rscript scripts/hb_smoke_test.R
```

Expected output: `SMOKE TEST PASSED` and wall time < 3 minutes.

If it fails, fix the environment (packages, C++ compilation) before proceeding.

---

## 6. Launch the production run

Create the log directory if needed:

```bash
mkdir -p logs
```

Launch with `nohup` so it survives logout, then `disown` so the shell doesn't
track it:

```bash
cd ~/medsat-pilot
nohup Rscript scripts/run_n1500_allpairs.R > logs/run_n1500_allpairs.log 2>&1 &
disown
```

Note the PID printed by the shell (e.g., `[1] 12345`) in case you need to
check or kill it later.

Expected settings:
- M=20, n_iter=3000, n_burn=1000, n_thin=1, n_draws=2000
- fit_interactions=TRUE (all-pairs: NO2×NDVI, NO2×IMD, NDVI×IMD)
- orthogonalize=TRUE, nu=1.0, log_rho_mu=log(0.2), log_rho_sd=1.0

---

## 7. Monitor progress

Tail the log in real time:

```bash
tail -f logs/run_n1500_allpairs.log
```

Check that the R process is still running:

```bash
ps aux | grep Rscript
```

Or by PID:

```bash
ps -p <PID>
```

The sampler prints iteration progress to stdout (captured in the log). You
should see lines like `[Gibbs] iter 100 / 3000 ...` at regular intervals.

---

## 8. Expected runtime

At n=1500, M=20, all-pairs interactions: roughly **2–5 hours** depending on
the node. The n=500 pilot (same settings) ran in ~30 min; n=1500 scales
roughly O(n^3) in the Cholesky step, so expect 3–10× longer.

Check timing in the log:
```
[fit] wall time: XXXX s  (XX.X min  |  XX.X h)
```

---

## 9. Retrieve results

Once the run completes (check the log for `PRODUCTION RUN complete`), the
output RDS is at:

```
~/medsat-pilot/output/fit_medsat_n1500_seed1_allpairs_hb.rds
```

Transfer back to your local machine (**file transfer tool: TBD**):

```bash
# From your local machine:
rsync -avz \
  <HELLBENDER_USER>@<HELLBENDER_HOST>:~/medsat-pilot/output/fit_medsat_n1500_seed1_allpairs_hb.rds \
  output/
```

Also retrieve the log for the record:

```bash
rsync -avz \
  <HELLBENDER_USER>@<HELLBENDER_HOST>:~/medsat-pilot/logs/run_n1500_allpairs.log \
  logs/
```

---

## 10. Placeholders to fill in before use

| Placeholder | What to fill in |
|---|---|
| `<HELLBENDER_USER>` | Your Hellbender username |
| `<HELLBENDER_HOST>` | Hellbender login node hostname |
| R module command | Output of `module avail R` |
| Personal R library path | e.g. `~/R/library` (confirm writable) |
| File transfer tool | scp / rsync / FileZilla / Globus / etc. |

---

## Notes

- The hostname guard in `run_n1500_allpairs.R` checks for "hellbender" in
  `Sys.info()[["nodename"]]`. If the node hostname does not match, the
  script stops immediately (before any compute). This prevents accidentally
  running the production script locally and overwriting nothing, but also
  means you need to confirm the actual node name.

- If you need to test the run script locally first (no compute committed),
  use the dry-run mode — it skips the hostname check and writes to a temp
  file:
  ```bash
  MEDSAT_DRY_RUN=1 Rscript scripts/run_n1500_allpairs.R
  ```

- The output RDS name `_hb.rds` is hardcoded (no host-tag logic). The
  production guard also prevents overwriting an existing fit — if the file
  already exists the script stops and asks you to rename or delete it first.
