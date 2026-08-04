/*!
 *  \file       test_ar1klt_gen.c
 *  \brief      Verification of the general-length recursive AR(1) KLT (ar1klt_gen_*) for N up to 1000.
 *
 *  \details
 *  The plan under test implements the recursive factorizations of Theorems III.3 and IV.4 of the paper
 *
 *      Y. Reznik, "Direct Factorization of the Karhunen-Loeve Transform of AR(1) Sources",
 *      submitted to the IEEE Transactions on Signal Processing, July 2026.
 *
 *  The verification has three parts:
 *
 *  Part 1, cross-check at N = 2..8: the recursive plan is materialized as a matrix and compared entrywise,
 *  up to the free row signs, against the fixed short-length modules of the library, which evaluate their
 *  constants by radical formulas.  Agreement confirms that two independent constant-evaluation paths meet.
 *
 *  Part 2, small lengths N = 2..64: the plan is compared against a canonical KLT computed by a cyclic
 *  Jacobi eigensolver (numeric root finding), with free row signs aligned.  The maximum entrywise
 *  deviation is reported, together with a mean squared error (MSE) over a random sample set.
 *
 *  Part 3, large lengths up to N = 1000: independent eigensolvers become the bottleneck, so correctness is
 *  established intrinsically.  The materialized transform must be orthogonal, must diagonalize the
 *  covariance, and must order the variances decreasingly; together these three properties characterize the
 *  KLT up to row signs.  The strict alternation of the two branch spectra is verified inside
 *  ar1klt_gen_init() at every recursion level.
 *
 *  The test exits with status 0 when all cases pass, so it can serve as a continuous-integration check.
 *
 *  \version    1.1
 *  \date       August 2026
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright  (C) 2026 Yuriy A. Reznik.
 *  \license    MIT License (see the LICENSE file).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ar1klt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*! \brief PASS threshold for entrywise matrix deviations. */
#define GTOL_DEV 1e-8
/*! \brief PASS threshold for orthogonality / diagonalization residues. */
#define GTOL_RES 1e-8
/*! \brief PASS threshold for sample MSE. */
#define GTOL_MSE 1e-16
/*! \brief Random sample vectors per MSE case. */
#define GNSAMP 2000

/* ********************************************************************************************** */
/*  Helpers                                                                                       */
/* ********************************************************************************************** */

/*!
 *  \brief      Materialize a recursive plan as an explicit n-by-n transform matrix.
 *
 *  \param[in]  plan   Plan built by ar1klt_gen_init().
 *  \param[in]  n      Transform length of the plan.
 *  \param[out] W      Transform matrix, stored row-major; obtained by applying the plan to unit vectors.
 *
 *  \return     0 on success; -1 when scratch memory could not be allocated.
 */
static int materialize(const ar1klt_gen *plan, int n, double *W)
{
    double *x = (double *)malloc((size_t)n * sizeof(double));
    double *y = (double *)malloc((size_t)n * sizeof(double));
    int i, j;
    if (x == NULL || y == NULL) { free(x); free(y); return -1; }
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) x[i] = 0.0;
        x[j] = 1.0;
        ar1klt_gen_apply(plan, x, y);
        for (i = 0; i < n; i++) W[i * n + j] = y[i];
    }
    free(x); free(y);
    return 0;
}

/*!
 *  \brief      Orthogonality residue of a materialized transform: the maximum entry of |W W^T - I|.
 *
 *  \param[in]  W   Transform matrix, stored row-major.
 *  \param[in]  n   Matrix order.
 *
 * \return     The maximum absolute deviation of W W^T from the identity matrix.
 */
static double orth_residue(const double *W, int n)
{
    int i, j, k;
    double worst = 0.0;
    for (i = 0; i < n; i++)
        for (j = i; j < n; j++) {
            double s = 0.0, e;
            for (k = 0; k < n; k++) s += W[i * n + k] * W[j * n + k];
            e = fabs(s - (i == j ? 1.0 : 0.0));
            if (e > worst) worst = e;
        }
    return worst;
}

/*!
 *  \brief      Diagonalization check of a materialized transform against the AR(1) covariance.
 *
 *  \details
 *  The routine forms E = W R W^T with R = [rho^|i-j|] and inspects it: a valid KLT must make E diagonal
 *  with a non-increasing diagonal, which lists the coefficient variances in decreasing order.
 *
 *  \param[in]  W      Transform matrix, stored row-major.
 *  \param[in]  n      Matrix order.
 *  \param[in]  rho    AR(1) correlation coefficient.
 *  \param[out] offd   Maximum absolute off-diagonal entry of E.
 *  \param[out] mono   1 when the diagonal of E is non-increasing (with tolerance 1e-9), 0 otherwise.
 *
 *  \return     0 on success; -1 when scratch memory could not be allocated.
 */
static int diag_residue(const double *W, int n, double rho, double *offd, int *mono)
{
    double *T = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    double *lam = (double *)malloc((size_t)n * sizeof(double));
    int i, j, k;
    if (T == NULL || lam == NULL) { free(T); free(lam); return -1; }
    /* T = W R, with R_kj = rho^|k-j| from a precomputed power table  */
    {
        double *rp = (double *)malloc((size_t)n * sizeof(double));
        if (rp == NULL) { free(T); free(lam); return -1; }
        rp[0] = 1.0;
        for (k = 1; k < n; k++) rp[k] = rho * rp[k - 1];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double s = 0.0;
                for (k = 0; k < n; k++)
                    s += W[i * n + k] * rp[k > j ? k - j : j - k];
                T[i * n + j] = s;
            }
        }
        free(rp);
    }
    *offd = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++) s += T[i * n + k] * W[j * n + k];
            if (i == j) lam[i] = s;
            else if (fabs(s) > *offd) *offd = fabs(s);
        }
    }
    *mono = 1;
    for (i = 0; i + 1 < n; i++)
        if (lam[i + 1] > lam[i] + 1e-9) *mono = 0;
    free(T); free(lam);
    return 0;
}

/*!
 *  \brief      Cyclic Jacobi eigensolver for a symmetric matrix held in heap storage.
 *
 *  \param[in,out] a   Symmetric input matrix of order n, row-major; destroyed during the iteration
 *                     (its diagonal converges to the eigenvalues).
 *  \param[out]    d   Eigenvalues, in the order left by the iteration (unsorted).
 *  \param[out]    v   Eigenvectors, stored as columns of the accumulated rotation product.
 *  \param[in]     n   Matrix order.
 */
static void jacobi_gen(double *a, double *d, double *v, int n)
{
    int i, j, k, sweep;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            v[i * n + j] = (i == j) ? 1.0 : 0.0;

    for (sweep = 0; sweep < 100; sweep++) {
        double off = 0.0;
        for (i = 0; i < n - 1; i++)
            for (j = i + 1; j < n; j++)
                off += a[i * n + j] * a[i * n + j];
        if (off < 1e-26) break;
        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n; j++) {
                double apq = a[i * n + j], th, t, c, s;
                if (fabs(apq) < 1e-300) continue;
                th = (a[j * n + j] - a[i * n + i]) / (2.0 * apq);
                t = (th >= 0.0)
                    ?  1.0 / (th + sqrt(1.0 + th * th))
                    : -1.0 / (-th + sqrt(1.0 + th * th));
                c = 1.0 / sqrt(1.0 + t * t);
                s = t * c;
                for (k = 0; k < n; k++) {
                    double aik = a[i * n + k], ajk = a[j * n + k];
                    a[i * n + k] = c * aik - s * ajk;
                    a[j * n + k] = s * aik + c * ajk;
                }
                for (k = 0; k < n; k++) {
                    double aki = a[k * n + i], akj = a[k * n + j];
                    a[k * n + i] = c * aki - s * akj;
                    a[k * n + j] = s * aki + c * akj;
                }
                for (k = 0; k < n; k++) {
                    double vki = v[k * n + i], vkj = v[k * n + j];
                    v[k * n + i] = c * vki - s * vkj;
                    v[k * n + j] = s * vki + c * vkj;
                }
            }
        }
    }

    for (i = 0; i < n; i++) d[i] = a[i * n + i];
}

/*!
 *  \brief      Uniform random variate in the open interval (0,1); never exactly zero or one.
 *  \return     One uniform variate.
 */
static double urand(void)
{
    return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

/*!
 *  \brief      Standard normal random variate, generated by the Box-Muller transformation.
 *  \return     One standard normal variate.
 */
static double grand(void)
{
    return sqrt(-2.0 * log(urand())) * cos(2.0 * M_PI * urand());
}

/* ********************************************************************************************** */
/*  Main program                                                                                  */
/* ********************************************************************************************** */

/*!
 *  \brief      Test driver: run all three verification parts and report a global result.
 */
int main(void)
{
    static const double rhos[] = { 0.10, 0.50, 0.90, 0.99 };
    static const int bigN[] = { 100, 127, 128, 255, 256, 257, 500, 512, 999, 1000 };
    static const double bigR[] = {0.1, 0.3, 0.50, 0.7, 0.90, 0.95, 0.99 };
    int fails = 0, n, r, i, j;

    srand(20260731);

    /* ---- Part 1: recursive plan vs fixed modules, N = 2..8 ------- */
    printf("Part 1: recursive plan vs fixed radical modules (N=2..8)\n");
    {
        double worst = 0.0;
        for (n = 2; n <= 8; n++) {
            for (r = 0; r < 4; r++) {
                ar1klt_gen *plan;
                ar1klt_ctx ctx;
                double Wg[8 * 8], x[8], y[8];
                if (ar1klt_gen_init(&plan, n, rhos[r]) != 0 ||
                    ar1klt_init(&ctx, n, rhos[r]) != 0) {
                    fails++; continue;
                }
                {
                    double Wm[8 * 8];
                    materialize(plan, n, Wg);
                    for (j = 0; j < n; j++) {
                        for (i = 0; i < n; i++) x[i] = 0.0;
                        x[j] = 1.0;
                        ar1klt_apply(&ctx, x, y);
                        for (i = 0; i < n; i++) Wm[i * n + j] = y[i];
                    }
                    /* rows are defined up to sign: align, then compare */
                    for (i = 0; i < n; i++) {
                        double dot = 0.0, sgn;
                        for (j = 0; j < n; j++)
                            dot += Wg[i * n + j] * Wm[i * n + j];
                        sgn = (dot < 0.0) ? -1.0 : 1.0;
                        for (j = 0; j < n; j++) {
                            double e = fabs(Wg[i * n + j]
                                            - sgn * Wm[i * n + j]);
                            if (e > worst) worst = e;
                        }
                    }
                }
                ar1klt_gen_free(plan);
            }
        }
        printf("  max |W_gen - W_module| = %.3e  %s\n\n", worst,
               worst < GTOL_DEV ? "pass" : "FAIL");
        if (!(worst < GTOL_DEV)) fails++;
    }

    /* ---- Part 2: vs Jacobi reference, N = 2..64 ------------------ */
    printf("Part 2: recursive plan vs Jacobi eigensolver (N=2..64)\n");
    {
        double worstdev = 0.0, worstmse = 0.0;
        for (n = 2; n <= 64; n++) {
            for (r = 0; r < 4; r++) {
                ar1klt_gen *plan;
                double *W, *A, *d, *V, *x, *yf, *yr;
                double dev = 0.0, acc = 0.0;
                int s, k;
                if (ar1klt_gen_init(&plan, n, rhos[r]) != 0) {
                    fails++; continue;
                }
                W = (double *)malloc((size_t)n * n * sizeof(double));
                A = (double *)malloc((size_t)n * n * sizeof(double));
                d = (double *)malloc((size_t)n * sizeof(double));
                V = (double *)malloc((size_t)n * n * sizeof(double));
                x = (double *)malloc((size_t)n * sizeof(double));
                yf = (double *)malloc((size_t)n * sizeof(double));
                yr = (double *)malloc((size_t)n * sizeof(double));
                materialize(plan, n, W);
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                        A[i * n + j] = pow(rhos[r],
                                           fabs((double)(i - j)));
                jacobi_gen(A, d, V, n);
                /* reference rows sorted by descending eigenvalue,
                 * signs aligned to the plan                          */
                for (i = 0; i < n; i++) {
                    int best = 0;
                    double dot;
                    for (j = 1; j < n; j++)
                        if (d[j] > d[best]) best = j;
                    d[best] = -1e300;
                    dot = 0.0;
                    for (j = 0; j < n; j++)
                        dot += V[j * n + best] * W[i * n + j];
                    for (j = 0; j < n; j++) {
                        double ref = (dot < 0.0 ? -1.0 : 1.0)
                                     * V[j * n + best];
                        double e = fabs(W[i * n + j] - ref);
                        if (e > dev) dev = e;
                        V[j * n + best] = ref; /* store aligned      */
                    }
                    /* also stash aligned row for the MSE test       */
                    for (j = 0; j < n; j++)
                        A[i * n + j] = V[j * n + best];
                }
                if (dev > worstdev) worstdev = dev;
                /* sample MSE                                        */
                for (s = 0; s < GNSAMP / 8; s++) {
                    x[0] = grand();
                    for (k = 1; k < n; k++)
                        x[k] = rhos[r] * x[k - 1]
                             + sqrt(1.0 - rhos[r] * rhos[r]) * grand();
                    ar1klt_gen_apply(plan, x, yf);
                    for (i = 0; i < n; i++) {
                        double v2 = 0.0;
                        for (j = 0; j < n; j++)
                            v2 += A[i * n + j] * x[j];
                        yr[i] = v2;
                    }
                    for (i = 0; i < n; i++)
                        acc += (yf[i] - yr[i]) * (yf[i] - yr[i]);
                }
                acc /= (double)((GNSAMP / 8) * n);
                if (acc > worstmse) worstmse = acc;
                free(W); free(A); free(d); free(V);
                free(x); free(yf); free(yr);
                ar1klt_gen_free(plan);
            }
        }
        printf("  max |W_gen - W_ref| = %.3e, max sample MSE = %.3e  %s\n\n",
               worstdev, worstmse,
               (worstdev < GTOL_DEV && worstmse < GTOL_MSE)
                   ? "pass" : "FAIL");
        if (!(worstdev < GTOL_DEV && worstmse < GTOL_MSE)) fails++;
    }

    /* ---- Part 3: intrinsic checks, large N ----------------------- */
    printf("Part 3: intrinsic checks at large N\n");
    printf("   N    rho    |WW'-I|      offdiag(WRW')  variances  result\n");
    for (i = 0; i < (int)(sizeof(bigN) / sizeof(bigN[0])); i++) {
        for (r = 0; r < (int)(sizeof(bigR) / sizeof(bigR[0])); r++) {
            ar1klt_gen *plan;
            double *W, orth, offd = 1.0;
            int mono = 0, bad;
            n = bigN[i];
            if (ar1klt_gen_init(&plan, n, bigR[r]) != 0) {
                printf("%5d  %.2f   init FAILED\n", n, bigR[r]);
                fails++; continue;
            }
            W = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
            if (W == NULL || materialize(plan, n, W) != 0) {
                printf("%5d  %.2f   alloc FAILED\n", n, bigR[r]);
                free(W); ar1klt_gen_free(plan); fails++; continue;
            }
            orth = orth_residue(W, n);
            diag_residue(W, n, bigR[r], &offd, &mono);
            bad = !(orth < GTOL_RES && offd < GTOL_RES && mono);
            fails += bad;
            printf("%5d  %.2f   %.3e    %.3e     %s   %s\n",
                   n, bigR[r], orth, offd,
                   mono ? "sorted" : "BAD   ", bad ? "FAIL" : "pass");
            free(W);
            ar1klt_gen_free(plan);
        }
    }

    printf("\n%s\n", fails ? "RESULT: FAIL" : "RESULT: PASS");
    return fails ? 1 : 0;
}
