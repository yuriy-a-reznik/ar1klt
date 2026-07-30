/*!
 *  \file ar1klt.c
 * 
 *  \brief Implementation of the exact short-length AR(1) KLT modules (N = 2..8) 
 * 
 *  \details
 *  This library implements the exact short-length KLT algorithms of
 *  Y. A. Reznik, "Direct Factorization of the Karhunen-Loeve Transform of AR(1) Sources", Section 5.
 * 
 *  Internal layout:
 * 
 *  Design-time constants are computed in ar1klt_init(): 
 *     rotation angles from closed-form arctangents; roots of the cubic and quartic characteristic 
 *     equations from the Cardano and Ferrari radical formulas, followed by a few Newton polishing iterations; 
 *     dense stage rows from the Cauchy forms of the rank-one (secular) and arrowhead eigenvector solutions, 
 *     normalized to unit length.
 *
 *  Run-time processing in ar1klt_apply():
 *     follows the printed stage chains literally.
 *
 *  All polynomial coefficient arrays are stored low order first, and all polynomials handled here are monic 
 *  with exclusively real roots (they are characteristic polynomials of small symmetric matrices).
 *
 *  \author    Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright Copyright (C) 2026 Yuriy A. Reznik
 *  \license   MIT License (see LICENSE file).
 */

#include <math.h>
#include <stddef.h>
#include "ar1klt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SQRT2 1.41421356237309504880

/********************
 *
 * Small helpers
 *
 */

/*! Plane rotation G(phi): (a,b) -> (a c + b s, -a s + b c). */
static void rot(double c, double s, double a, double b, double *u, double *v)
{
    *u =  c * a + s * b;
    *v = -s * a + c * b;
}

/*! Clamp x into [lo, hi]. */
static double clampd(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/*! Real cube root (sign-preserving). */
static double cbrt_r(double x)
{
    if (x >= 0.0) return pow(x, 1.0 / 3.0);
    return -pow(-x, 1.0 / 3.0);
}

/*********************
 * 
 * Polynomial utilities
 *
 */

/*! Polynomial product c = a * b; lc = la + lb - 1. */
static void pmul(const double *a, int la, const double *b, int lb, double *c)
{
    int i, j;
    for (i = 0; i < la + lb - 1; i++) c[i] = 0.0;
    for (i = 0; i < la; i++)
        for (j = 0; j < lb; j++)
            c[i + j] += a[i] * b[j];
}

/*! In-place scaled add: a += s * b (lb <= la). */
static void paxpy(double *a, double s, const double *b, int lb)
{
    int i;
    for (i = 0; i < lb; i++) a[i] += s * b[i];
}

/*! Evaluate polynomial of length len at t (Horner). */
static double peval(const double *c, int len, double t)
{
    int i;
    double v = c[len - 1];
    for (i = len - 2; i >= 0; i--) v = v * t + c[i];
    return v;
}

/*! Evaluate derivative of polynomial of length len at t. */
static double pdeval(const double *c, int len, double t)
{
    int i;
    double v = (len - 1) * c[len - 1];
    for (i = len - 2; i >= 1; i--) v = v * t + i * c[i];
    return v;
}

/*! Newton-polish nr roots of the polynomial (5 iterations each). */
static void polish(const double *c, int len, double *r, int nr)
{
    int i, it;
    for (i = 0; i < nr; i++) {
        for (it = 0; it < 5; it++) {
            double d = pdeval(c, len, r[i]);
            if (fabs(d) < 1e-300) break;
            r[i] -= peval(c, len, r[i]) / d;
        }
    }
}

/*! Sort n reals ascending (insertion sort). */
static void sort_asc(double *r, int n)
{
    int i, j;
    for (i = 1; i < n; i++) {
        double t = r[i];
        for (j = i - 1; j >= 0 && r[j] > t; j--) r[j + 1] = r[j];
        r[j + 1] = t;
    }
}

/***********************
 * 
 * Radical root solvers (all-real-root cases):
 * 
 **/

/*!
 *  \brief Roots of a monic cubic equation t^3 + c[2] t^2 + c[1] t + c[0] = 0 
 *         (coefficients low order first, c[3] = 1), via Cardano.
 *
 *  Handles both the three-real-root (trigonometric) and the one-real-root (cbrt) branch; returns the number of real roots
 *  written to r[] in ascending order.
 */
static int cubic_roots(const double *c, double r[3])
{
    double b = c[2], p, q, sh, disc;
    int n;
    sh = -b / 3.0;                       /* depression shift            */
    p  = c[1] - b * b / 3.0;
    q  = 2.0 * b * b * b / 27.0 - b * c[1] / 3.0 + c[0];
    disc = -4.0 * p * p * p - 27.0 * q * q;
    if (disc >= 0.0 && p < 0.0) {        /* three real roots            */
        double m  = 2.0 * sqrt(-p / 3.0);
        double ac = clampd(3.0 * q / (p * m), -1.0, 1.0);
        double th = acos(ac) / 3.0;
        r[0] = m * cos(th) + sh;
        r[1] = m * cos(th - 2.0 * M_PI / 3.0) + sh;
        r[2] = m * cos(th - 4.0 * M_PI / 3.0) + sh;
        n = 3;
    } else {                             /* one real root               */
        double d  = q * q / 4.0 + p * p * p / 27.0;
        double sd = sqrt(fabs(d));
        if (d >= 0.0)
            r[0] = cbrt_r(-q / 2.0 + sd) + cbrt_r(-q / 2.0 - sd) + sh;
        else {                           /* numerically borderline      */
            double m  = 2.0 * sqrt(fabs(p) / 3.0);
            double ac = clampd(3.0 * q / (p * m), -1.0, 1.0);
            r[0] = m * cos(acos(ac) / 3.0) + sh;
        }
        n = 1;
    }
	/* polish the roots (5 Newton iterations each) and sort ascending */
    polish(c, 4, r, n);
    sort_asc(r, n);
    return n;
}

/*!
 *  \brief Roots of a monic quartic equation t^4 + c[3] t^3 + ... + c[0] = 0,
 *         with four real roots, via Ferrari; ascending order in r[].
 */
static void quartic_roots(const double *c, double r[4])
{
    double b = c[3], sh, p, q, rr, res[4], m3[3], m, s, d1, d2;
    int nm;
    sh = -b / 4.0;                       /* depression shift            */
    p  = c[2] - 3.0 * b * b / 8.0;
    q  = c[1] - b * c[2] / 2.0 + b * b * b / 8.0;
    rr = c[0] - b * c[1] / 4.0 + b * b * c[2] / 16.0
         - 3.0 * b * b * b * b / 256.0;
    if (fabs(q) < 1e-14) {               /* biquadratic case            */
        double dd = sqrt(clampd(p * p - 4.0 * rr, 0.0, 1e300));
        double y1 = clampd((-p + dd) / 2.0, 0.0, 1e300);
        double y2 = clampd((-p - dd) / 2.0, 0.0, 1e300);
        r[0] =  sqrt(y1) + sh;  r[1] = -sqrt(y1) + sh;
        r[2] =  sqrt(y2) + sh;  r[3] = -sqrt(y2) + sh;
    } else {
        /* resolvent cubic: m^3 + p m^2 + ((p^2 - 4 r)/4) m - q^2/8 = 0 */
        res[0] = -q * q / 8.0;
        res[1] = (p * p - 4.0 * rr) / 4.0;
        res[2] = p;
        res[3] = 1.0;
        nm = cubic_roots(res, m3);
        m  = m3[nm - 1];                 /* largest root; m > 0 here    */
        m  = clampd(m, 1e-300, 1e300);
        s  = sqrt(2.0 * m);
        /* split into two quadratics, discriminants clamped to >= 0    */
        d1 = s * s - 4.0 * (p / 2.0 + m + q / (2.0 * s));
        d2 = s * s - 4.0 * (p / 2.0 + m - q / (2.0 * s));
        d1 = sqrt(clampd(d1, 0.0, 1e300));
        d2 = sqrt(clampd(d2, 0.0, 1e300));
        r[0] = ( s + d1) / 2.0 + sh;
        r[1] = ( s - d1) / 2.0 + sh;
        r[2] = (-s + d2) / 2.0 + sh;
        r[3] = (-s - d2) / 2.0 + sh;
    }
    /* polish the roots (5 Newton iterations each) and sort ascending */
    polish(c, 5, r, 4);
    sort_asc(r, 4);
}

/********************
 * 
 * Characteristic polynomials of the Section-5 modules
 * 
 */

/*!
 *  \brief Residual cubic chi(t) = (1-t)[(e-t)^2 - rho^2] - rho^2 (e-t),
 *         e = 1 + rho^2, returned monic (low order first, length 4).
 */
static void chi_residual(double rho, double *c)
{
    double e = 1.0 + rho * rho, r2 = rho * rho;
    double lin1[2], line[2], sq[3], t3[4];
    lin1[0] = 1.0; lin1[1] = -1.0;       /* (1 - t)  */
    line[0] = e;   line[1] = -1.0;       /* (e - t)  */
    pmul(line, 2, line, 2, sq);          /* (e-t)^2  */
    sq[0] -= r2;                         /* - rho^2  */
    pmul(lin1, 2, sq, 3, t3);            /* (1-t)[..] */
    paxpy(t3, -r2, line, 2);             /* - rho^2 (e-t) */
    c[0] = -t3[0]; c[1] = -t3[1]; c[2] = -t3[2]; c[3] = 1.0; /* monic  */
}

/*! \brief Parity quadratic (1-t)(tau-t) - rho^2, monic (length 3). */
static void chi_parity(double rho, double tau, double *c)
{
    c[0] = tau - rho * rho;
    c[1] = -(1.0 + tau);
    c[2] = 1.0;
}

/********************
 *
 * Dense stage assembly
 *
 */

/*!
 * \brief Rank-one secular rows: Q[i][j] ~ p[j] / (mu[j] - nu[i]),
 *        normalized; m poles, m roots.
 */
static void secular_stage(int m, const double *mu, const double *p, const double *nu, double Q[4][4])
{
    int i, j;
    for (i = 0; i < m; i++) {
        double s = 0.0;
        for (j = 0; j < m; j++) {
            Q[i][j] = p[j] / (mu[j] - nu[i]);
            s += Q[i][j] * Q[i][j];
        }
        s = 1.0 / sqrt(s);
        for (j = 0; j < m; j++) Q[i][j] *= s;
    }
}

/*!
 * \brief Arrowhead rows: A[i][j] ~ sqrt(2) rho q[j] / (dl[j] - nu[i]) for j < m, last entry 1; normalized; m poles, m+1 roots.
 */
static void arrow_stage(int m, double rho, const double *dl, const double *q, const double *nu, double A[4][5])
{
    int i, j;
    for (i = 0; i < m + 1; i++) {
        double s = 1.0;                  /* the trailing 1              */
        for (j = 0; j < m; j++) {
            A[i][j] = SQRT2 * rho * q[j] / (dl[j] - nu[i]);
            s += A[i][j] * A[i][j];
        }
        A[i][m] = 1.0;
        s = 1.0 / sqrt(s);
        for (j = 0; j <= m; j++) A[i][j] *= s;
    }
}

/*!
 *  \brief Init function.
 * 
 *  Offline processing of the AR(1) KLT module for a given length and correlation coefficient.
 */
int ar1klt_init(ar1klt_ctx *ctx, int n, double rho)
{
    double e, r2, ts, ta;
    if (ctx == NULL || n < AR1KLT_NMIN || n > AR1KLT_NMAX) return -1;
    if (!(rho > 0.0 && rho < 1.0)) return -1;

    ctx->n = n;
    ctx->rho = rho;
    e  = 1.0 + rho * rho;
    r2 = rho * rho;
    ts = 1.0 - rho + r2;
    ta = 1.0 + rho + r2;

    ctx->phi3  = 0.5 * atan(2.0 * SQRT2 / rho);
    ctx->phi_s = 0.5 * atan(-(1.0 - rho) / 2.0);
    ctx->phi_a = 0.5 * atan( (1.0 + rho) / 2.0);
    ctx->psi   = 0.5 * atan(rho / 2.0);

    if (n == 5) {
        /* delta_{1,2} = 1 + rho^2/2 -/+ sqrt(rho^2 + rho^4/4)          */
        double dl[2], q[2], nu[3], c4[4], lin1[2], line[2], t3[4];
        double cp = cos(ctx->psi), sp = sin(ctx->psi);
        dl[0] = 1.0 + r2 / 2.0 - sqrt(r2 + r2 * r2 / 4.0);
        dl[1] = 1.0 + r2 / 2.0 + sqrt(r2 + r2 * r2 / 4.0);
        q[0]  =  (cp - sp) / SQRT2;
        q[1]  = -(cp + sp) / SQRT2;
        /* cubic (1-t)(e-t)^2 - rho^2 [2(1-t) + (e-t)] = 0 (monic)      */
        lin1[0] = 1.0; lin1[1] = -1.0;
        line[0] = e;   line[1] = -1.0;
        pmul(line, 2, line, 2, t3);      /* (e-t)^2, length 3           */
        pmul(lin1, 2, t3, 3, c4);        /* (1-t)(e-t)^2, length 4      */
        paxpy(c4, -2.0 * r2, lin1, 2);
        paxpy(c4, -r2, line, 2);
        c4[0] = -c4[0]; c4[1] = -c4[1]; c4[2] = -c4[2]; c4[3] = 1.0;
        cubic_roots(c4, nu);
        arrow_stage(2, rho, dl, q, nu, ctx->Aar);
    }

    if (n == 6 || n == 7) {
        double mu[3], p[3];
        double c3 = cos(ctx->phi3), s3 = sin(ctx->phi3);
        mu[0] = 1.0 + r2 / 2.0 - sqrt(2.0 * r2 + r2 * r2 / 4.0);
        mu[1] = 1.0;
        mu[2] = 1.0 + r2 / 2.0 + sqrt(2.0 * r2 + r2 * r2 / 4.0);
        p[0] =  c3 / SQRT2;
        p[1] = -1.0 / SQRT2;
        p[2] = -s3 / SQRT2;
        if (n == 6) {
            /* branch cubics (1-t)[(e-t)(tau-t) - rho^2] - rho^2(tau-t) */
            double nu[3], c4[4], lin1[2], line[2], lint[2], t3[3], t4[4];
            int b;
            for (b = 0; b < 2; b++) {
                double tau = (b == 0) ? ts : ta;
                lin1[0] = 1.0; lin1[1] = -1.0;
                line[0] = e;   line[1] = -1.0;
                lint[0] = tau; lint[1] = -1.0;
                pmul(line, 2, lint, 2, t3);
                t3[0] -= r2;
                pmul(lin1, 2, t3, 3, t4);
                paxpy(t4, -r2, lint, 2);
                c4[0] = -t4[0]; c4[1] = -t4[1]; c4[2] = -t4[2]; c4[3] = 1.0;
                cubic_roots(c4, nu);
                secular_stage(3, mu, p, nu,
                              (b == 0) ? ctx->Qs : ctx->Qa);
            }
        } else {
            /* residual cubic and bordered quartic                      */
            double chi[4], dl[3], q[3], nu[4], line[2], c5[5], t5[5];
            double lin2[3];
            int i, j;
            chi_residual(rho, chi);
            cubic_roots(chi, dl);
            secular_stage(3, mu, p, dl, ctx->Qd);
            for (i = 0; i < 3; i++) {
                q[i] = 0.0;
                for (j = 0; j < 3; j++) q[i] += ctx->Qd[i][j] * p[j];
            }
            /* quartic chi(t)(e-t) - 2 rho^2 [(1-t)(e-t) - rho^2] = 0.
             * chi_residual() returns the monic cubic, which is the
             * NEGATION of the paper's chi (leading term -t^3); hence
             * the monic quartic is  -(chi_monic*(e-t)) - 2 rho^2 G.  */
            line[0] = e; line[1] = -1.0;
            pmul(chi, 4, line, 2, t5);
            lin2[0] = e - r2; lin2[1] = -(1.0 + e); lin2[2] = 1.0;
            c5[0] = -t5[0]; c5[1] = -t5[1]; c5[2] = -t5[2];
            c5[3] = -t5[3]; c5[4] = 1.0;
            paxpy(c5, -2.0 * r2, lin2, 3);
            quartic_roots(c5, nu);
            arrow_stage(3, rho, dl, q, nu, ctx->Aar);
        }
    }

    if (n == 8) {
        double mu[4], p[4], chs[3], cha[3], chi[4], t5[5], c5[5], nu[4];
        double cs = cos(ctx->phi_s), ss = sin(ctx->phi_s);
        double ca = cos(ctx->phi_a), sa = sin(ctx->phi_a);
        double sig_s = rho * (1.0 - rho), sig_a = rho * (1.0 + rho);
        int b;
        mu[0] = 1.0 - sig_s / 2.0 - (rho / 2.0) * sqrt((1.0 - rho) * (1.0 - rho) + 4.0);
        mu[1] = 1.0 + sig_a / 2.0 - (rho / 2.0) * sqrt((1.0 + rho) * (1.0 + rho) + 4.0);
        mu[2] = 1.0 - sig_s / 2.0 + (rho / 2.0) * sqrt((1.0 - rho) * (1.0 - rho) + 4.0);
        mu[3] = 1.0 + sig_a / 2.0 + (rho / 2.0) * sqrt((1.0 + rho) * (1.0 + rho) + 4.0);
        p[0] = 0.5 * (cs + ss);
        p[1] = -0.5 * (ca + sa);
        p[2] = 0.5 * (cs - ss);
        p[3] = 0.5 * (sa - ca);
        chi_parity(rho, ts, chs);
        chi_parity(rho, ta, cha);
        chi_residual(rho, chi);
        for (b = 0; b < 2; b++) {
            /* paper: chi_s chi_a + sig * chi_paper with sig = -sig_s
             * (sum) or +sig_a (diff); chi_residual() returns the
             * monic cubic = -chi_paper, so the coefficient flips.    */
            double sig = (b == 0) ? sig_s : -sig_a;
            pmul(chs, 3, cha, 3, t5);    /* chi_s * chi_a, length 5     */
            paxpy(t5, sig, chi, 4);      /* + sig * (-chi_paper)        */
            c5[0] = t5[0]; c5[1] = t5[1]; c5[2] = t5[2];
            c5[3] = t5[3]; c5[4] = 1.0;  /* already monic               */
            quartic_roots(c5, nu);
            secular_stage(4, mu, p, nu, (b == 0) ? ctx->Qs : ctx->Qa);
        }
    }

    return 0;
}

/*************************
 *
 *  Run-time processing
 *
 */

/*! Inner 3-point KLT stage (N = 3, 6, 7): x[3] -> g[3]. */
static void w3_apply(double phi3, const double *x, double *g)
{
    double c = cos(phi3), s = sin(phi3);
    double u = (x[0] + x[2]) / SQRT2;
    g[1] = (x[0] - x[2]) / SQRT2;
    rot(c, s, u, x[1], &g[0], &g[2]);
}

/*! Inner 4-point KLT stage (N = 4, 8): x[4] -> g[4]. */
static void w4_apply(double phi_s, double phi_a, const double *x, double *g)
{
    double cs = cos(phi_s), ss = sin(phi_s);
    double ca = cos(phi_a), sa = sin(phi_a);
    double u0 = (x[0] + x[3]) / SQRT2, u1 = (x[1] + x[2]) / SQRT2;
    double v0 = (x[0] - x[3]) / SQRT2, v1 = (x[1] - x[2]) / SQRT2;
    double s0 = (u0 + u1) / SQRT2, s1 = (u0 - u1) / SQRT2;
    double d0 = (v0 + v1) / SQRT2, d1 = (v0 - v1) / SQRT2;
    rot(cs, ss, s0, s1, &g[0], &g[2]);
    rot(ca, sa, d0, d1, &g[1], &g[3]);
}

/*! Dense m x m stage: y = Q x. */
static void dense_apply(int m, const double Q[4][4], const double *x, double *y)
{
    int i, j;
    for (i = 0; i < m; i++) {
        double s = 0.0;
        for (j = 0; j < m; j++) s += Q[i][j] * x[j];
        y[i] = s;
    }
}

/*! Arrowhead (m+1) x (m+1) stage: y = A (x, c). */
static void arrow_apply(int m, const double A[4][5], const double *x, double c, double *y)
{
    int i, j;
    for (i = 0; i < m + 1; i++) {
        double s = A[i][m] * c;
        for (j = 0; j < m; j++) s += A[i][j] * x[j];
        y[i] = s;
    }
}

/*! 
 *  \brief Apply the AR(1) KLT module to a vector x of length n, producing y. 
 * 
 *  Online processing of the AR(1) KLT module for a given length and correlation coefficient.
 */
void ar1klt_apply(const ar1klt_ctx *ctx, const double *x, double *y)
{
    switch (ctx->n) {

    case 2: {
        y[0] = (x[0] + x[1]) / SQRT2;
        y[1] = (x[0] - x[1]) / SQRT2;
        break;
    }

    case 3: {
        double g[3];
        w3_apply(ctx->phi3, x, g);
        y[0] = g[0]; y[1] = g[1]; y[2] = g[2];
        break;
    }

    case 4: {
        double g[4];
        w4_apply(ctx->phi_s, ctx->phi_a, x, g);
        y[0] = g[0]; y[1] = g[1]; y[2] = g[2]; y[3] = g[3];
        break;
    }

    case 5: {
        double u[2], v[2], a0, a1, b0, b1, r[2], t[2], z[3];
        double cp = cos(ctx->psi), sp = sin(ctx->psi);
        int i;
        for (i = 0; i < 2; i++) {
            u[i] = (x[i] + x[4 - i]) / SQRT2;
            v[i] = (x[i] - x[4 - i]) / SQRT2;
        }
        a0 = (u[0] + u[1]) / SQRT2; a1 = (u[0] - u[1]) / SQRT2;
        b0 = (v[0] + v[1]) / SQRT2; b1 = (v[0] - v[1]) / SQRT2;
        rot(cp, sp, a0, a1, &r[0], &r[1]);
        rot(cp, sp, b0, b1, &t[0], &t[1]);
        arrow_apply(2, ctx->Aar, r, x[2], z);
        y[0] = z[0]; y[2] = z[1]; y[4] = z[2];
        y[1] = t[0]; y[3] = t[1];
        break;
    }

    case 6: {
        double u[3], v[3], g[3], h[3], zs[3], za[3];
        int i;
        for (i = 0; i < 3; i++) {
            u[i] = (x[i] + x[5 - i]) / SQRT2;
            v[i] = (x[i] - x[5 - i]) / SQRT2;
        }
        w3_apply(ctx->phi3, u, g);
        w3_apply(ctx->phi3, v, h);
        dense_apply(3, ctx->Qs, g, zs);
        dense_apply(3, ctx->Qa, h, za);
        y[0] = zs[0]; y[2] = zs[1]; y[4] = zs[2];
        y[1] = za[0]; y[3] = za[1]; y[5] = za[2];
        break;
    }

    case 7: {
        double u[3], v[3], g[3], h[3], r[3], t[3], z[4];
        int i;
        for (i = 0; i < 3; i++) {
            u[i] = (x[i] + x[6 - i]) / SQRT2;
            v[i] = (x[i] - x[6 - i]) / SQRT2;
        }
        w3_apply(ctx->phi3, u, g);
        w3_apply(ctx->phi3, v, h);
        dense_apply(3, ctx->Qd, g, r);
        dense_apply(3, ctx->Qd, h, t);
        arrow_apply(3, ctx->Aar, r, x[3], z);
        y[0] = z[0]; y[2] = z[1]; y[4] = z[2]; y[6] = z[3];
        y[1] = t[0]; y[3] = t[1]; y[5] = t[2];
        break;
    }

    case 8: {
        double u[4], v[4], g[4], h[4], zs[4], za[4];
        int i;
        for (i = 0; i < 4; i++) {
            u[i] = (x[i] + x[7 - i]) / SQRT2;
            v[i] = (x[i] - x[7 - i]) / SQRT2;
        }
        w4_apply(ctx->phi_s, ctx->phi_a, u, g);
        w4_apply(ctx->phi_s, ctx->phi_a, v, h);
        dense_apply(4, ctx->Qs, g, zs);
        dense_apply(4, ctx->Qa, h, za);
        y[0] = zs[0]; y[2] = zs[1]; y[4] = zs[2]; y[6] = zs[3];
        y[1] = za[0]; y[3] = za[1]; y[5] = za[2]; y[7] = za[3];
        break;
    }

    default:
        break;
    }
}
