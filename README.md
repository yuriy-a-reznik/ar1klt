# AR1KLT

Version 1.1.

Exact fast algorithms for the Karhunen–Loève transform (KLT) of stationary
AR(1) sources, written in C.  The project includes two implementations:
exact short-length transforms for **N = 2, …, 8**, and a **general-N**
recursive implementation, with a test program verifying it up to
**N = 1000**.

The KLT of an AR(1) source with correlation coefficient `rho` is the
eigenvector basis of its covariance matrix `R_N(rho) = [rho^|i-j|]` — the
optimal decorrelating transform of classical transform coding.  This
library implements factorizations of the *exact* KLT transform, as
derived in:

> Y. Reznik, *"Direct Factorization of the Karhunen–Loève Transform of
> AR(1) Sources"*, submitted to the IEEE Transactions on Signal Processing,
> July 2026.

In the short-length modules, all algorithm constants are **algebraic
functions of `rho`**: rotation angles are closed-form arctangents, and the
remaining constants are roots of explicit cubic and quartic characteristic
equations, evaluated here by the Cardano–Ferrari radical formulas (with a
few Newton polishing iterations for floating-point hygiene).  Constants are
precomputed once per `(N, rho)` by `ar1klt_init()`; `ar1klt_apply()` then
runs the stage chain on data blocks.

As `rho -> 1` the modules converge to the corresponding fast DCT-II
modules (for `N = 4`, to Chen's factorization), and as `rho -> 0` to
DST-I modules.

## Layout

```
AR1KLT/
├── include/ar1klt.h          Public API (Doxygen-annotated)
├── src/ar1klt.c              Fixed short-length modules, N = 2..8
├── src/ar1klt_gen.c          General-N recursive factorization
├── test/test_ar1klt.c        Verification of the fixed modules
├── test/test_ar1klt_gen.c    Verification of the recursion (N up to 1000)
├── Makefile                  Build: library objects, test binaries, docs
├── Doxyfile                  Doxygen configuration
├── LICENSE                   MIT
└── README.md                 This file
```

## General-N recursive factorization

Beyond the fixed modules, the library implements the paper's recursive
factorizations for **arbitrary length** (`ar1klt_gen_init` /
`ar1klt_gen_apply` / `ar1klt_gen_free`):

* even orders split into a butterfly, **one shared half-length plan
  applied to both branches**, and two dense orthogonal correction
  stages;
* odd orders split into a butterfly, the shared half-length plan, a
  dense residual correction applied to both branches, and an arrowhead
  stage absorbing the center sample;
* outputs are interleaved by merging the two branch spectra in
  ascending generator-eigenvalue order; the merge is verified to
  alternate strictly between the branches at every recursion level.

Design-time constants (secular and arrowhead eigenvalues) are computed
by bracketed bisection on the interlacing intervals.  The correction
stages are stored densely in this reference implementation; their
Cauchy structure admits fast application (see the paper, Section on
complexity), which this package does not exploit.

```c
ar1klt_gen *plan;
double x[1000], y[1000];

if (ar1klt_gen_init(&plan, 1000, 0.95) != 0) { /* handle error */ }
ar1klt_gen_apply(plan, x, y);
ar1klt_gen_free(plan);
```

## Building and testing

Any C compiler will do; the code has no dependencies beyond `libc`/`libm`.

```sh
make            # builds the test binaries (strict flags: -std=c89 -pedantic -Wall -Wextra)
make test       # runs both verification suites
make docs       # generates API documentation with Doxygen (optional)
```

`test_ar1klt` verifies every fixed module against a **canonical KLT**
computed independently by numeric root finding: it builds `R_N(rho)`,
diagonalizes it with a cyclic Jacobi eigensolver, and sorts the
eigenvectors by decreasing eigenvalue.  For each `(N, rho)` case it
reports

* the maximum entrywise deviation between the fast module (materialized
  as a matrix) and the canonical KLT matrix, after aligning the free row
  signs; and
* the **MSE** between the two coefficient sets over a randomly generated
  sample set of Gaussian AR(1) vectors (5000 vectors per case).

Typical output (deviations at machine-precision level):

```
  N   rho      max|W_fast - W_ref|   sample MSE      result
  --  ----     -------------------   -------------   ------
   5  0.90     1.615e-15             5.305e-31       pass
   8  0.95     1.095e-13             5.015e-28       pass
...
RESULT: PASS
```

`test_ar1klt_gen` verifies the recursive factorization in three parts:
(1) against the fixed radical modules for N = 2..8, up to the free row
signs; (2) against a Jacobi eigensolver for every N = 2..64; and
(3) intrinsically for N up to 1000 (including odd and prime lengths),
checking orthogonality, diagonalization of the covariance, and the
decreasing variance order --- observed residues are at the 1e-13 to
1e-11 level throughout.

Both binaries exit with status 0 on success, so `make test` can serve
as a CI check.

## Usage

```c
#include "ar1klt.h"

ar1klt_ctx ctx;
double x[8], y[8];

if (ar1klt_init(&ctx, 8, 0.95) != 0) { /* handle bad arguments */ }
/* ... fill x ... */
ar1klt_apply(&ctx, x, y);   /* y[0..7] = KLT coefficients, decreasing variance */
```

Rows of the transform are defined up to sign, as usual for eigenvector
bases; `ar1klt_init()` fixes one deterministic convention.

## License

MIT — see [LICENSE](LICENSE).
