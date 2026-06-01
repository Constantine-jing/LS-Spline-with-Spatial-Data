# Sample size comparison

Diagnostic study comparing $n=400$ vs $n=1000$. ✅ Complete. Used to confirm the [[X1 undercoverage is basis resolution]] insight.

## Purpose

Determine whether observed undercoverage of $f_1 = 2\sin(\pi x)$ at $M=6$ knots is:

(a) a model failure (prior, sampler, or implementation issue), or
(b) a basis-resolution / approximation-bias issue that disappears with more flexibility

Logic: if it's a sampler or prior issue, more data won't help. If it's basis resolution, undercoverage stays at $n=1000$ but disappears when you increase $M$.

## Setup

- Two sample sizes: $n = 400$, $n = 1000$
- Otherwise identical setup to [[Sim1 n1000]]

Script: `run_sample_size_comparison.R`

## Result

[Fill in specific numbers]

The pattern of $f_1$ undercoverage at $M=6$ persisted across both sample sizes, **ruling out** sampler / prior issues. Resolved by raising $M$ to 20 — confirming basis resolution is the cause.

## Why this is a satisfying diagnostic

It's a hypothesis-testing approach to a coverage problem, not "vary parameters and hope it works." See [[X1 undercoverage is basis resolution]] for the full reasoning.

## Connects to

- [[X1 undercoverage is basis resolution]]
- [[Sim1 n1000]]
