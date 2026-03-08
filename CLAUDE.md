# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git Workflow

Do not commit new features directly to `main`. Create a new branch for any new feature work and commit there instead.

## Overview

This is a [SageMath](https://www.sagemath.org/) research codebase for computing the **average polynomial-time algorithm** for the Hasse-Witt matrix (A_f) of hyperelliptic curves over finite fields. The work follows Harvey's lecture notes on computing zeta functions of hyperelliptic curves efficiently.

The primary mathematical object is `A_f`: a g×g matrix (where g = (d-1)/2, d = degree of f) that encodes the action of Frobenius on the de Rham cohomology of the curve y² = f(x).

## Running Code

All scripts are SageMath files (`.sage`). Run them with:

```bash
sage <filename>.sage
```

Scripts in `Harvey-Class/` use absolute load paths (e.g., `/home/janeshi/AWS26/...`), so they may need path adjustment when run locally. Scripts in `good-code/` use relative `load()` paths and are self-contained.

Tests must be run from the **repo root**:

```bash
sage tests/runtests.sage
```

## Repository Structure

- **`good-code/`** — Current, production-quality implementations:
  - `P7-utils.sage` — Core data structures and tree algorithms (load this first)
  - `P7-avg-poly.sage` — Main algorithm: average polynomial-time A_f computation across all primes up to N
  - `P6-Af.sage` — Benchmark/reference: naive O(p) computation of A_f via direct polynomial exponentiation
  - `P4-2-9-benchmark.sage` — Benchmarking utilities

- **`Harvey-Class/`** — Earlier iterations and exercises following Harvey's chapter-by-chapter problems:
  - `P7-prec-2/` — Precision-2 (mod p²) implementation
  - `P7-prec-3/` — Precision-3 (mod p³) implementation, extends ZP2 → ZP3
  - `Ch2Ex.sage`, `Ch3Ex.sage`, etc. — Chapter exercises

## Key Algorithms and Data Structures

### `ZP2` / `ZP3` (in `P7-utils.sage`)
Represent elements of R[P]/(P^μ) as tuples. `ZP2(x, y)` = x + y·P. The `realize(p)` method evaluates at P=p to get a concrete integer. These enable symbolic product trees that work for all primes simultaneously before specializing.

### `ZP2Matrix` / `ZP3Matrix`
Matrices over R[P]/(P^μ), stored as μ separate Sage matrices (the coefficient matrices of 1, P, P², ...). Arithmetic is done component-wise using Sage's native matrix routines for speed.

### Product Tree / Accumulating Remainder Tree
- `build_product_tree(list)` — builds a binary tree of products
- `remainder_tree_builder(value_tree, modulus_tree, identity)` — accumulating remainder tree that computes `value mod modulus` at each leaf simultaneously for all primes in O(N log² N) rather than O(N·p)

### Main Algorithm Flow (`compute_A_f_avg_poly`)
1. Build product trees of `T_bar` matrices (over ZP2) for i=0,1 (the two "shifts")
2. Use accumulating remainder tree to reduce mod p² for each prime p at once
3. Interpolate between the i=0 and i=1 reductions to get sprint matrices for arbitrary multiples of p
4. For each prime p > 5, iterate through l = 0..g-1 computing the accumulator vector via sprint steps and single T-matrix steps
5. Extract A_f[p] as last g entries of the accumulator

### Coefficient Convention
`F_coeffs` is the coefficient list of f(x) (as returned by `f.list()` in Sage), so `F_coeffs[0]` is the leading coefficient and `F_coeffs[d]` is the constant term.

### Performance Notes
- Average polynomial-time algorithm (good-code): ~59 seconds for all primes up to 50,000
- Naive sqrt-time algorithm (P6-Af.sage): ~30 minutes for the same range
- The test polynomial used throughout: `f = -(x^12 - x^10 + 6*x^9 - 7*x^8 + 5*x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1)` (genus 5)
