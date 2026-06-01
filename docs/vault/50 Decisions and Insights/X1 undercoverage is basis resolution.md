# X1 undercoverage is basis resolution

A diagnostic story. Important because it's a case where the temptation is to blame the model, but the right answer is "the basis can't represent the truth well enough at $M=6$."

## The observation

In several early simulations, the credible interval coverage for $f_1(x) = 2\sin(\pi x)$ (or its centered version, post-fix) was below nominal — around 80% instead of 95% — when knot count was $M=6$.

Coverage of $f_2, f_3$ was fine. Only $f_1$ was bad.

## The two hypotheses

**(H1) Model failure.** Something in the [[Collapsed Gibbs sampler]], [[RW2 prior]], or implementation is wrong. Coverage should improve with more data because the posterior concentrates around the *correct* answer.

**(H2) Basis-resolution / approximation bias.** With only $M=6$ knots, the LS basis can't represent $\sin(\pi x)$ exactly — there's a fixed approximation error that doesn't shrink with $n$. The posterior concentrates around the *best-the-basis-can-do* answer, not the truth, and intervals miss the truth systematically.

## The diagnostic that distinguished them

[[Sample size comparison]]: run the same model at $n=400$ and $n=1000$.

- Under H1: more data should improve coverage
- Under H2: more data shrinks intervals around the wrong center → coverage gets *worse*, not better

**Result:** coverage stayed the same (or got slightly worse) at $n=1000$. This rules out H1 and confirms H2.

## The fix

Increase $M$ to 20. With more knots, the LS approximation to $\sin(\pi x)$ is essentially exact, and coverage recovers.

## Why this is satisfying methodologically

It's a **principled diagnostic**, not just "vary settings until coverage looks ok." The two hypotheses make different predictions; we ran the experiment that distinguishes them; the data picked one.

## Why this matters for the writeup

- Don't sell $M=6$ results as a coverage failure
- Do mention $M$ as a tuning parameter with a known (and explainable) regime where approximation bias dominates
- Phase 2 simulations should use $M$ large enough to put approximation bias well below MCMC noise

## Connects to

- [[Sample size comparison]]
- [[LS basis]]
- [[RW2 prior]]
- [[Centering for ANOVA identifiability]] (the *other* coverage gotcha)
