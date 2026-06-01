# Prior sensitivity is a strength

A reframe that turns a defensive analysis into a positive contribution.

## The standard frame (defensive)

"We did a prior sensitivity analysis to show our results don't depend on prior choices."

This sounds like covering your bases. It's correct but unmotivated.

## The better frame

"At $n=1000$, posterior inference is essentially identical across five prior bundles — including a deliberately misspecified $\rho$ prior (B5). This is a robustness result: the data dominate the prior in this regime, and reviewers can have confidence that downstream inferences (component estimates, coverage, model selection) are not artifacts of prior choice."

## Why this is dissertation-quality

1. **Negative results are rare and valuable.** Showing a sensitivity analysis where things *do* differ across priors would just be confusing. Showing they don't differ — and explaining *why* (data dominance) — is a finding.
2. **B5 is the active ingredient.** A deliberately misspecified prior is not a routine sanity check; it's a stress test. Robustness *to misspecification* is stronger than robustness within a reasonable range.
3. **It pairs with [[REML vs Bayes]].** REML doesn't have priors; Bayes does. Showing Bayes matches REML in point estimates AND is robust to prior choice is the case-closing argument that the Bayesian apparatus isn't smuggling in answers.

## How to present it

- Lead with the conclusion: "Inference is robust across B1–B5."
- Show the visual: posterior curves overlaid; they're indistinguishable.
- Point to B5 specifically: "Even under deliberate misspecification of the spatial range prior, the smooth components are recovered with the same accuracy and coverage."
- Finish with the caveat: "Robustness here reflects strong data ($n=1000$). At smaller $n$, prior choice should be revisited."

## Things to NOT say

- ❌ "Prior sensitivity analysis confirmed our priors were appropriate." (sounds defensive, suggests you tried priors until they worked)
- ❌ "We chose B1 because it gave the best fit." (this is unprincipled and reviewers will notice)

## Connects to

- [[Prior sensitivity B1-B5]]
- [[Prior sensitivity Sim4]]
- [[REML vs Bayes]]
