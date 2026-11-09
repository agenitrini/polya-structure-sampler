# polyagen — uniform random sampling of Pólya structures


A C99 program (no dependency other than `libm`) that draws uniformly at
random — up to floating-point precision — an unlabelled structure of a given
size `n`, in one of the following five families:

| family     | structures                                              | size                | OEIS    |
|------------|---------------------------------------------------------|---------------------|---------|
| `polya`    | unlabelled non-plane rooted trees                       | number of nodes     | A000081 |
| `hier`     | hierarchies (no unary node)                             | number of leaves    | A000669 |
| `identity` | identity trees (children pairwise non-isomorphic)       | number of nodes     | A004111 |
| `binary`   | non-plane binary trees (0 or 2 children per node), Wedderburn–Etherington numbers | number of leaves    | A001190 |
| `ubinary`  | unrooted binary trees (internal nodes of degree 3)      | number of leaves (`n ≥ 2`) | A129860 |

This code accompanies the paper

> Antoine Genitrini, *Smooth disturbances of tree functions: full asymptotic
> expansions and uniform random sampling of Pólya structures*.

An interactive visualisation of the samplers is available at
<https://agenitrini.github.io/dynamic_trees/polya.html>.

## Building and usage

```sh
make            # builds ./polyagen (cc -O3 -march=native)
make test       # self-tests (./polyagen --selftest)
make bench      # table of measurements (./bench.sh)
```

```
polyagen --family polya|hier|identity|binary|ubinary -n SIZE [--count K] [--seed S]
         [--format parent|brackets|dot|none] [--stats] [--kexact K0]
polyagen --selftest
```

* `-n` accepts the notation `1e6`;
* `--format parent` (default): one line per structure, the array of parents
  of the nodes numbered in preorder (the root is node 0, with parent `-1`);
  `--format brackets`: bracketed word, e.g. `((()())(()()))` (for
  `binary`, a leaf is written `()` and an internal node `(` child₁ child₂ `)`;
  the tree has `2n-1` nodes);
  `--format dot`: [Graphviz](https://graphviz.org/) format, one directed
  graph `digraph treeK { ... }` per structure (K = 0, 1, …), edges
  parent `->` child, nodes numbered in preorder (root 0), without
  labels; e.g. `./polyagen --family polya -n 100 --format dot | dot -Tpdf -o tree.pdf`
  (with `--count K > 1`, `dot` produces one page per structure);
  `--format none`: no output (for measurements);
* `--stats` (on `stderr`): for each structure, number of nodes, check of the
  size and of the class constraints, number of steps (examined pairs), ratio
  `steps/(n ln n)`, number of rejections (identity), time; then the means
  over the `--count` samples;
* `--kexact K0`: threshold between exact and asymptotic tables (default 3000).

Example:

```sh
$ ./polyagen --family polya -n 7 --count 2 --format brackets --seed 1
((()())(()()))
(((()()()())))
$ ./polyagen --family binary -n 5 --count 2 --format brackets --seed 1
(((()())())(()()))
(((()())())(()()))
```

## Method

### Recursive method, "Cayley-first" order

For `polya`, a tree of size `n` is a root together with a forest (multiset
of trees) of size `n-1`; for `hier`, a hierarchy of size `n ≥ 2` is a root
together with a multiset of at least two hierarchies of total size `n`. A
multiset of size `m` is drawn as follows (Nijenhuis–Wilf's RANRUT argument):
choose a pair `(d,j)` with probability

    P(d,j) = d · A_d · F_{m-jd} / (m · F_m),

recursively draw a uniform element `σ` of size `d` and a uniform multiset
`φ` of size `m - jd`, and return `φ` plus `j` copies of `σ`. The pairs are
scanned in the *Cayley-first* order: first all pairs with `j = 1` by
decreasing `d` (`d = m, m-1, …, 1`), then the pairs with `j ≥ 2` (here `d`
increasing then `j` increasing; the order of this tail is arbitrary). The
expected cost is about `(1/2) n ln n` examined pairs (theorem of the paper);
this number of "steps" is measured.
For the root of a hierarchy, the pair `(n,1)` is forbidden and the
normalisation becomes `n (F_n - T_n) = n T_n`.

For `identity`, a tree of size `n` is obtained by choosing `i` with
probability `i T_i V_i(n-i) / ((n-1) T_n)` (scan `i = n-1, …, 1`), drawing a
uniform identity tree `σ` of size `i`, then a tree `τ` of size `n-i` until
the root of `τ` has no child isomorphic to `σ` (rejection), and grafting `σ`
under the root of `τ`. The isomorphism test is exact and takes `O(1)`: each
completed subtree receives a canonical identifier by *hash-consing* (sorted
list of the children's identifiers, so no collision is possible). The steps
of rejected attempts are counted in the cost.

For `binary` (`T(z) = z + T(z)²/2 + T(z²)/2`, size = number of leaves), we
use the pointed form `z T' = z ζ'(z) + T · z T'` with `ζ(z) = z + T(z²)/2`,
that is, for `n ≥ 2`,

    n T_n = Σ_{d=1}^{n-1} d T_d T_{n-d} + [n even] (n/2) T_{n/2}.

A tree of size `n ≥ 2` is obtained either by a *Cayley step* `d`
(`1 ≤ d ≤ n-1`, probability `d T_d T_{n-d} / (n T_n)`): a root whose two
children are a uniform tree `σ` of size `d` and a uniform tree `φ` of size
`n-d`, drawn independently; or (if `n` is even) by the *decoration step*
(probability `T_{n/2} / (2 T_n)`): a root whose two children are two copies
of the same uniform tree `σ` of size `n/2`. A tree `{A,B}` with `A ≠ B` is
reached via `d = |A|` and `d = |B|`, with total weight `|A| + |B| = n`; a
tree `{A,A}` via `d = n/2` and via the decoration, with total weight `n`:
each tree therefore has probability `n/(n T_n) = 1/T_n`. *Cayley-first*
order: `d = n-1, n-2, …, 1`, then the decoration last; the number of examined
steps is measured as for the other families.

Implementation: explicit stack (no C recursion), tree stored in preorder as
the array `sz[p]` of subtree sizes; this representation is
translation-invariant, so the `j-1` copies of `σ` are plain `memcpy` calls
on a contiguous block (likewise for the second copy of `σ` in the decoration
step of `binary`). Pseudo-random generator xoshiro256** seeded by
splitmix64, uniform floats with 53 bits.

### Normalised tables

The integers `T_k` (which overflow) are never used; instead we use the
normalised weights `w_k = ρ^k T_k ~ k^{-3/2}` (and `f_k = ρ^k F_k`). The
probabilities can be rewritten without any exploding power of `ρ`:

    P(d,j) = d · w_d · f_{m-jd} · ρ^{(j-1)d} / (m · f_m),

with, for `polya`, `f_m = w_{m+1}` (since `F_m = T_{m+1}`), and for `hier`,
`f_0 = 1`, `f_1 = ρ`, `f_k = 2 w_k` (`k ≥ 2`). For `identity`,
`v_i(m) = ρ^m V_i(m) = Σ_j (-1)^j ρ^{ji} w_{m-ji}` is computed on the fly.
For `binary`, `P(d) = d w_d w_{n-d} / (n w_n)` and
`P(dec) = (n/2) w_{n/2} ρ^{n/2} / (n w_n)`.

* For `k ≤ K0` (default `K0 = min(n, 3000)`), `w_k` is computed by the
  normalised recurrence in Euler convolution form,
  `(n-1) w_n = Σ_k β_k w_{n-k}` with `β_k = Σ_{d|k} ε(k/d) d w_d ρ^{k-d}`
  (equivalent to the "(d,j)" recurrences of the paper, and in `O(K0²)`);
  for `binary`, `w_1 = ρ` and
  `2 w_n = Σ_{i=1}^{n-1} w_i w_{n-i} + [n even] ρ^{n/2} w_{n/2}`.
  The computations are carried out in **double-double** (~106 bits) rather
  than in `long double`: on Apple Silicon `long double` has only 53 bits.
* For `k > K0`, we use the asymptotic expansion
  `w_k ≈ (π k³)^{-1/2} Σ_{ℓ=0}^{15} τ_ℓ k^{-ℓ}`; `ρ` and the `τ_ℓ` (file
  `asymptotic_constants.json`, 25 digits) are embedded in the source.
  Agreement with the recurrence is checked at `k = K0` (warning if the
  relative gap exceeds `1e-15`). Measured: gap ≈ `2e-18` (polya),
  `3e-19` (hier), `8e-17` (identity), `6e-17` (binary) from `k = 120` on, and
  ≈ `2e-22`–`4e-22` at `k = 3000` (limited by the 25 digits of the constants).
  The constants for `binary` come from the computations of the paper
  (high-precision evaluation of the full expansion).

The precomputation costs `O(K0² + n)`: about 0.05 s for `n = 10^6`.

### Unrooted binary trees (`ubinary`)

An unrooted binary tree with `n` leaves is rooted at its *leaf centroid*:
either the unique internal node whose three branches have fewer than `n/2`
leaves each, or, for `n` even, the unique edge that splits the leaves into
two halves. The tree is thus a multiset of three rooted binary trees of
sizes `< n/2` and total size `n`, or a multiset of two rooted binary trees
of size `n/2`. These multisets are drawn with the cycle index of `S_3`
(resp. `S_2`): a pair (permutation, sequence of trees fixed by it) is drawn
uniformly, which makes the multiset uniform (Burnside's lemma), and the
trees themselves are drawn by the `binary` sampler. The output is the tree
rooted at the centroid node (or at an endpoint of the centroid edge), with
`2n-2` nodes.

The exact weights of the ordered triples of sizes `< n/2` require an
`O(n²)` precomputation, done once for all the samples of a given size:
about 0.7 s for `n = 10^5` and 6 s for `n = 3·10^5` on an Apple M4 Pro.
The generation itself costs the same as for `binary`.

### Uniformity: up to floating-point precision

The probabilities are handled as `double` (53-bit draw of `U`, successive
subtractions). The sampler is therefore uniform **up to rounding**: each
branching probability carries a relative error of order `10^-16`, and the
output distribution differs from the uniform distribution by an amount of
that order per decision. If, because of rounding, `U` exceeds the computed
sum of the probabilities, `U` is redrawn (an event of probability
~`10^-16`). An **exactly** uniform sampler would require certified bounds
and lazy refinement of the precision, in the spirit of Denise and Zimmermann
(*Uniform random generation of decomposable structures using floating-point
arithmetic*, TCS 1999); this is not implemented here.

## Tests (`make test`)

1. **(a)** exact recurrences (128-bit integers) = OEIS terms `T_1..T_20` of
   the four families; for `binary`, check of the pointed form
   `n T_n = Σ d T_d T_{n-d} + [n even] (n/2) T_{n/2}` for `n ≤ 30`;
   double-double tables = `ρ^n T_n` up to `~1e-31` (`n ≤ 30`); agreement
   between recurrence and asymptotics for `k` from 60 to 3000;
2. **(d)** the sum of the pair probabilities equals 1 (`< 1e-12`, measured
   `≤ 3e-15`) for sizes on either side of `K0`;
3. **(c)** validity of the structures (size, no unary node, no isomorphic
   siblings; for `binary`: `n` leaves, `2n-1` nodes, each node has 0 or 2
   children) for `n` from 1 to 50,000, including with `K0 = 300`
   (asymptotic tables);
4. **(b)** uniformity: 2·10^5 samples per size, `polya` and `hier` for
   `n = 5..9`, `identity` and `binary` for `n = 5..10`; we check that all
   isomorphism classes appear (e.g. 48 trees for polya `n = 7`, 33
   hierarchies for `n = 6`, 25 identity trees for `n = 8`, 52 for `n = 9`,
   23 binary trees for `n = 8`, 98 for `n = 10`) and run a chi-squared test
   (failure if `p < 1e-4`);
5. **(e)** root statistic at `n = 25` (exact probability of a unary root, of
   a binary root for hierarchies, and for `binary` of a root split into
   `{1, n-1}`, i.e. `T_24/T_25`): deviation < 4.5 standard deviations.

6. **(f)** `ubinary`: the centroid decomposition gives the OEIS terms of
   A129860 for `n = 2..20`; uniformity (2·10^5 samples, chi-squared test,
   canonical forms computed at the centre of the tree) for `n = 6..12`
   (e.g. 11 trees for `n = 10`, 37 for `n = 12`); validity (`n` leaves,
   `2n-2` nodes, internal degree 3) for `n` from 2 to 5000.

Result: **all tests pass** (chi-squared p-values between 0.08 and 0.92;
total running time ≈ 3 s).

## Measurements

Apple M4 Pro, Apple clang 21, `-O3 -march=native`, `--seed 2024`,
`--format none`. Mean time to generate one structure (excluding the
precomputation of the tables, which is reported separately and done once
for all samples); `steps` = number of examined pairs `(d,j)` (or sizes `i`).
Table produced by `./bench.sh` (`FAMILIES=binary ./bench.sh` for a single
family). For `binary`, a "step" is an examined Cayley step `d` or the
examined decoration step.

| family | n | samples | mean time / structure (s) | steps/(n ln n) | precomputation (s) |
|---|---:|---:|---:|---:|---:|
| polya | 1000 | 2000 | 0.000026 | 0.6177 | 0.0049 |
| polya | 2000 | 1000 | 0.000050 | 0.6101 | 0.0193 |
| polya | 10000 | 500 | 0.000253 | 0.5890 | 0.0434 |
| polya | 100000 | 100 | 0.002620 | 0.5582 | 0.0436 |
| polya | 1000000 | 20 | 0.028387 | 0.5791 | 0.0457 |
| polya | 10000000 | 3 | 0.307049 | 0.5896 | 0.0663 |
| hier | 1000 | 2000 | 0.000035 | 0.7696 | 0.0048 |
| hier | 2000 | 1000 | 0.000068 | 0.7459 | 0.0195 |
| hier | 10000 | 500 | 0.000344 | 0.7061 | 0.0431 |
| hier | 100000 | 100 | 0.003583 | 0.6637 | 0.0432 |
| hier | 1000000 | 20 | 0.037880 | 0.6707 | 0.0460 |
| hier | 10000000 | 3 | 0.396087 | 0.6476 | 0.0760 |
| identity | 1000 | 2000 | 0.000039 | 0.5932 | 0.0049 |
| identity | 2000 | 1000 | 0.000077 | 0.5876 | 0.0192 |
| identity | 10000 | 500 | 0.000387 | 0.5744 | 0.0430 |
| identity | 100000 | 100 | 0.004128 | 0.5545 | 0.0433 |
| identity | 1000000 | 20 | 0.043574 | 0.5533 | 0.0452 |
| identity | 10000000 | 3 | 0.566181 | 0.5689 | 0.0660 |
| binary | 1000 | 2000 | 0.000032 | 0.6031 | 0.0025 |
| binary | 2000 | 1000 | 0.000062 | 0.5912 | 0.0097 |
| binary | 10000 | 500 | 0.000312 | 0.5788 | 0.0219 |
| binary | 100000 | 100 | 0.003260 | 0.5629 | 0.0223 |
| binary | 1000000 | 20 | 0.034257 | 0.5518 | 0.0240 |
| binary | 10000000 | 3 | 0.363016 | 0.5541 | 0.0449 |

Remarks: the ratio `steps/(n ln n)` slowly decreases towards 1/2 (≈ 0.61 at
`n = 2000` for polya); for large sizes, few samples are drawn and the number
of steps has a high variance (it is dominated by the first levels of the
recursion), hence the non-monotone values. For `identity`, the mean number
of rejections is about `0.019 n` (≈ 19,000 for `n = 10^6`).
For `binary`, the precomputation is twice as fast (direct, symmetrised
quadratic recurrence) and over 30 samples at `n = 10^7` we measure
`steps/(n ln n) ≈ 0.547`.
Memory for `n = 10^7`: ≈ 120 MB (polya), 230 MB (hier), 530 MB (identity,
including the hash-consing table), 160 MB (binary, `2n-1` nodes).

## Files

* `polyagen.c` — the sampler and the self-tests;
* `Makefile` — targets `all`, `test`, `bench`, `clean`;
* `bench.sh` — measurement script (Markdown table);
* `LICENSE` — text of the GNU GPL v3 licence;
* `asymptotic_constants.json` — constants `ρ`, `τ_ℓ` (computed
  in high precision from the full asymptotic expansions of the paper).

## Licence

Copyright (C) 2026 Antoine Genitrini.

This program is free software, distributed under the terms of the
[GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html)
as published by the Free Software Foundation, either version 3 or (at your
option) any later version (`SPDX-License-Identifier: GPL-3.0-or-later`). It
is distributed without any warranty. The full text of the licence is in the
file [`LICENSE`](LICENSE).
