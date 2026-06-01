# Bug log

A running list of bugs, gotchas, and the fixes. Add new entries at the top with date.

---

## 2026 — Centering of $f_1$ in interaction DGP

**Symptom:** [[Interaction 2x2 v2]] showed inflated RMSE for $f_{12}$ and degraded coverage everywhere.

**Cause:** $f_1 = 2\sin(\pi x)$ is not mean-zero on $[0,1]$ (integral is $4/\pi$), so the data-generating function violates the [[ANOVA identifiability]] constraint imposed during fitting.

**Fix:** Replace with $f_1 = 2(\sin(\pi x) - 2/\pi)$ in the DGP. After fix: 3–5× RMSE_$f_{12}$ improvement, coverage ≥90%.

**Lesson:** [[Centering for ANOVA identifiability]] — when simulating from an additive model with sum-to-zero constraints, the DGP components must also be mean-zero.

---

## 2026 — `vapply` vs `sapply` in `build_interaction_prior_precision`

**Symptom:** Type error in production simulation loop, only triggered when no interaction terms were active for a particular configuration.

**Cause:** `sapply` on an empty list returns a list (not a numeric of length 0), violating downstream type assumptions.

**Fix:** Use `vapply` with explicit return type and length.

```r
# Before (buggy)
result <- sapply(empty_list, FUN)

# After
result <- vapply(empty_list, FUN, numeric(1))
```

**Lesson:** Prefer `vapply` over `sapply` in any production loop where the list might be empty. The type signature documents the contract.

---

## 2026 — Off-by-one in `col_map_main`

**Symptom:** Misaligned columns in the design matrix when the interaction term was added; main-effect coefficients pulled wrong indices.

**Cause:** Off-by-one indexing in the `col_map_main` lookup that translates per-component coefficient indices to global design-matrix columns.

**Fix:** Corrected the index arithmetic. Verified by checking that fitted main-effect curves match the no-interaction baseline when the DGP has no interaction.

---

## DOS line endings on `.sub` files

**Symptom:** SBATCH submission fails with cryptic parse error.

**Cause:** Files edited on Windows have `\r\n` line endings; SBATCH expects `\n`.

**Fix:**
```bash
sed -i 's/\r$//' my_job.sub
```

Add this to your "before submitting" checklist.

---

## Template for new entries

```markdown
## YYYY-MM-DD — Short title

**Symptom:** What broke.

**Cause:** Root cause once you found it.

**Fix:** What you actually changed.

**Lesson:** Generalizable principle, if any. Link to a [[Decision]] note if it became a methodological insight.
```

## Connects to

- [[R file index]]
- [[Centering for ANOVA identifiability]]
- [[Hellbender HPC notes]]
