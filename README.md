# AR1KLT

Exact fast algorithms for the Karhunen–Loève transform (KLT) of stationary
AR(1) sources, at short lengths **N = 2, …, 8**, in strict ANSI C (C89).

The KLT of an AR(1) source with correlation coefficient `rho` is the
eigenvector basis of its covariance matrix `R_N(rho) = [rho^|i-j|]` — the
optimal decorrelating transform of classical transform coding.  This
library implements the *exact* KLT, at the actual `rho` of the source, as
short chains of butterflies and plane rotations followed by small dense
stages, as derived in:

> Y. A. Reznik, *"Direct Factorization of the Karhunen–Loève Transform of
> AR(1) Sources"*, Section 5 (exact short-length KLT factorizations).

All algorithm constants are **algebraic functions of `rho`**: rotation
angles are closed-form arctangents, and the remaining constants are roots
of explicit cubic and quartic characteristic equations, evaluated here by
the Cardano–Ferrari radical formulas (with a few Newton polishing
iterations for floating-point hygiene).  Constants are precomputed once
per `(N, rho)` by `ar1klt_init()`; `ar1klt_apply()` then runs the stage
chain on data blocks.

As `rho -> 1` the modules converge to the corresponding fast DCT-II
modules (for `N = 4`, to Chen's factorization), and as `rho -> 0` to
DST-I modules.

## Layout

```
AR1KLT/
├── include/ar1klt.h      Public API (Doxygen-annotated)
├── src/ar1klt.c          Library implementation
├── test/test_ar1klt.c    Verification & MSE measurement program
├── Makefile              Build: library object, test binary, docs
├── Doxyfile              Doxygen configuration
├── LICENSE               MIT
└── README.md             This file
```

## Building and testing

Any C89 compiler will do; the code has no dependencies beyond `libc`/`libm`.

```sh
make            # builds the test binary (strict flags: -std=c89 -pedantic -Wall -Wextra)
make test       # runs the verification suite
make docs       # generates API documentation with Doxygen (optional)
```

The test program verifies every module against a **canonical KLT**
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

The binary exits with status 0 on success, so `make test` can serve as a
CI check.

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
