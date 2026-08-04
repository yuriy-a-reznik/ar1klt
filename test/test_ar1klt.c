/*!
 *  \file test_ar1klt.c
 * 
 *  \brief Verification of the AR1KLTs fast modules against a canonical
 *         AR(1) KLT computed by numeric eigendecomposition.
 *
 *  \details
 *  For each supported length N = 2..8 and a sweep of correlation
 *  coefficients rho, the test
 *
 *   1. Builds the covariance matrix R_N(rho) = [rho^|i-j|] and computes
 *      its full eigendecomposition numerically with a cyclic Jacobi
 *      eigensolver (the "canonical" KLT: eigenvectors sorted by
 *      decreasing eigenvalue, i.e., numeric root finding of the
 *      characteristic equation to machine precision);
 * 
 *   2. Materializes the fast module of the library as a matrix (by
 *      applying it to unit vectors), aligns the free row signs, and
 *      reports the maximum entrywise deviation from the canonical KLT;
 * 
 *   3. Generates a random sample set of Gaussian AR(1) vectors, pushes
 *      each sample through both transforms, and reports the mean squared
 *      error (MSE) between the two coefficient sets.
 *
 *  The test exits with status 0 when every case passes the thresholds
 *  and 1 otherwise, so it can serve as a CI check.
 * 
 *  \author    Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright Copyright (C) 2026 Yuriy A. Reznik
 *  \license   MIT License (see LICENSE file).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ar1klt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*! \brief Number of random sample vectors per (N, rho) test case. */
#define NSAMPLES 5000

/*! \brief PASS threshold on the max entrywise matrix deviation. */
#define TOL_DEV 1e-9

/*! \brief PASS threshold on the sample-set MSE. */
#define TOL_MSE 1e-18

/*!
 * \brief Cyclic Jacobi eigensolver for a symmetric n x n matrix.
 *
 * \param[in,out] a  On input the symmetric matrix; destroyed on output (diagonal converges to the eigenvalues).
 * \param[out]    d  Eigenvalues.
 * \param[out]    v  Eigenvectors, stored as columns.
 * \param[in]     n  Matrix order (n <= AR1KLT_NMAX).
 */
static void jacobi_eig(double a[AR1KLT_NMAX][AR1KLT_NMAX], double *d, double v[AR1KLT_NMAX][AR1KLT_NMAX], int n)
{
    int i, j, k, sweep;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            v[i][j] = (i == j) ? 1.0 : 0.0;
    for (sweep = 0; sweep < 100; sweep++) {
        double off = 0.0;
        for (i = 0; i < n - 1; i++)
            for (j = i + 1; j < n; j++)
                off += a[i][j] * a[i][j];
        if (off < 1e-28) break;
        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n; j++) {
                double apq = a[i][j], theta, t, c, s;
                if (fabs(apq) < 1e-300) continue;
                theta = (a[j][j] - a[i][i]) / (2.0 * apq);
                t = (theta >= 0.0)
                    ?  1.0 / (theta + sqrt(1.0 + theta * theta))
                    : -1.0 / (-theta + sqrt(1.0 + theta * theta));
                c = 1.0 / sqrt(1.0 + t * t);
                s = t * c;
                for (k = 0; k < n; k++) {
                    double aik = a[i][k], ajk = a[j][k];
                    a[i][k] = c * aik - s * ajk;
                    a[j][k] = s * aik + c * ajk;
                }
                for (k = 0; k < n; k++) {
                    double aki = a[k][i], akj = a[k][j];
                    a[k][i] = c * aki - s * akj;
                    a[k][j] = s * aki + c * akj;
                }
                for (k = 0; k < n; k++) {
                    double vki = v[k][i], vkj = v[k][j];
                    v[k][i] = c * vki - s * vkj;
                    v[k][j] = s * vki + c * vkj;
                }
            }
        }
    }
    for (i = 0; i < n; i++) d[i] = a[i][i];
}

/*!
 * \brief Canonical AR(1) KLT: rows of W are the eigenvectors of R_N(rho) = [rho^|i-j|] sorted by decreasing eigenvalue.
 *
 * \param[in]  n    Transform length.
 * \param[in]  rho  Correlation coefficient.
 * \param[out] W    Canonical KLT matrix, rows in KLT order.
 */
static void klt_reference(int n, double rho, double W[AR1KLT_NMAX][AR1KLT_NMAX])
{
    double R[AR1KLT_NMAX][AR1KLT_NMAX];
    double V[AR1KLT_NMAX][AR1KLT_NMAX];
    double d[AR1KLT_NMAX];
    int idx[AR1KLT_NMAX];
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            R[i][j] = pow(rho, fabs((double)(i - j)));
    jacobi_eig(R, d, V, n);
    for (i = 0; i < n; i++) idx[i] = i;
    for (i = 0; i < n - 1; i++)          /* sort descending by d       */
        for (j = i + 1; j < n; j++)
            if (d[idx[j]] > d[idx[i]]) {
                k = idx[i]; idx[i] = idx[j]; idx[j] = k;
            }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            W[i][j] = V[j][idx[i]];      /* row i = eigenvector i      */
}


/*! \brief Uniform variate in (0,1), never exactly 0 or 1. */
static double urand(void)
{
    return ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
}

/*! \brief Standard normal variate (Box-Muller). */
static double grand(void)
{
    return sqrt(-2.0 * log(urand())) * cos(2.0 * M_PI * urand());
}

/*!
 * \brief Draw one Gaussian AR(1) sample vector of length n.
 *
 * \param[in]  n    Vector length.
 * \param[in]  rho  Correlation coefficient.
 * \param[out] x    Sample vector (unit-variance stationary AR(1)).
 */
static void ar1_sample(int n, double rho, double *x)
{
    int k;
    x[0] = grand();
    for (k = 1; k < n; k++)
        x[k] = rho * x[k - 1] + sqrt(1.0 - rho * rho) * grand();
}


/*!
 * \brief Run one (n, rho) verification case.
 *
 * \param[in]  n       Transform length.
 * \param[in]  rho     Correlation coefficient.
 * \param[out] maxdev  Max entrywise deviation of the module matrix
 *                     from the sign-aligned canonical KLT matrix.
 * \param[out] mse     Sample-set MSE between the two coefficient sets.
 * \return 0 on PASS, 1 on FAIL.
 */
static int run_case(int n, double rho, double *maxdev, double *mse)
{
    ar1klt_ctx ctx;
    double Wref[AR1KLT_NMAX][AR1KLT_NMAX];
    double Wmod[AR1KLT_NMAX][AR1KLT_NMAX];
    double x[AR1KLT_NMAX], y[AR1KLT_NMAX];
    double yf[AR1KLT_NMAX], yr[AR1KLT_NMAX];
    double dev, acc;
    int i, j, s;

    if (ar1klt_init(&ctx, n, rho) != 0) {
        *maxdev = -1.0; *mse = -1.0;
        return 1;
    }

    /* materialize the module as a matrix                              */
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) x[i] = 0.0;
        x[j] = 1.0;
        ar1klt_apply(&ctx, x, y);
        for (i = 0; i < n; i++) Wmod[i][j] = y[i];
    }

    /* canonical reference, with free row signs aligned to the module  */
    klt_reference(n, rho, Wref);
    for (i = 0; i < n; i++) {
        double dot = 0.0;
        for (j = 0; j < n; j++) dot += Wref[i][j] * Wmod[i][j];
        if (dot < 0.0)
            for (j = 0; j < n; j++) Wref[i][j] = -Wref[i][j];
    }

    dev = 0.0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double e = fabs(Wmod[i][j] - Wref[i][j]);
            if (e > dev) dev = e;
        }
    *maxdev = dev;

    /* sample-set MSE                                                  */
    acc = 0.0;
    for (s = 0; s < NSAMPLES; s++) {
        ar1_sample(n, rho, x);
        ar1klt_apply(&ctx, x, yf);
        for (i = 0; i < n; i++) {
            double v = 0.0;
            for (j = 0; j < n; j++) v += Wref[i][j] * x[j];
            yr[i] = v;
        }
        for (i = 0; i < n; i++)
            acc += (yf[i] - yr[i]) * (yf[i] - yr[i]);
    }
    *mse = acc / ((double)NSAMPLES * (double)n);

    return (dev < TOL_DEV && *mse < TOL_MSE) ? 0 : 1;
}


/*!
 * \brief Test driver: sweep all lengths and a set of rho values.
 * \return 0 if all cases pass; 1 otherwise.
 */
int main(void)
{
    static const double rhos[] = { 0.10, 0.30, 0.50, 0.70, 0.90, 0.95, 0.99 };
    const int nrho = (int)(sizeof(rhos) / sizeof(rhos[0]));
    int n, r, fails = 0;

    srand(20260730);

    printf("AR1KLTs verification against canonical Jacobi-based KLT\n");
    printf("(%d random AR(1) sample vectors per case)\n\n", NSAMPLES);
    printf("  N   rho      max|W_fast - W_ref|   sample MSE      result\n");
    printf("  --  ----     -------------------   -------------   ------\n");

    for (n = AR1KLT_NMIN; n <= AR1KLT_NMAX; n++) {
        for (r = 0; r < nrho; r++) {
            double dev, mse;
            int bad = run_case(n, rhos[r], &dev, &mse);
            fails += bad;
            printf("  %2d  %.2f     %.3e             %.3e       %s\n",
                   n, rhos[r], dev, mse, bad ? "FAIL" : "pass");
        }
    }

    printf("\n%s\n", fails ? "RESULT: FAIL" : "RESULT: PASS");
    return fails ? 1 : 0;
}
