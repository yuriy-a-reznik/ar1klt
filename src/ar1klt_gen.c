/*!
 *  \file       ar1klt_gen.c
 *  \brief      General-length recursive AR(1) KLT: the self-similar factorizations for arbitrary N >= 1.
 *
 *  \details
 *  This file implements the recursive factorizations of Theorems III.3 and IV.4 of the paper
 *
 *      Y. Reznik, "Direct Factorization of the Karhunen-Loeve Transform of AR(1) Sources",
 *      submitted to the IEEE Transactions on Signal Processing, July 2026.
 *
 *  Design time (ar1klt_gen_init) proceeds bottom-up.  Each recursion node stores a pointer to ONE shared
 *  half-length plan, used by both the sum branch and the difference branch; the dense correction stages;
 *  the output interleaving; the generator eigenvalues in output order; and the vector of last components of
 *  the transform rows, which is the fold vector consumed by the parent level.  Secular and arrowhead
 *  eigenvalues are found by bracketed bisection on the interlacing intervals.  Eigenvector rows follow the
 *  Cauchy forms of the paper and are normalized to unit length.
 *
 *  Run time (ar1klt_gen_apply) executes the stage chain literally: butterfly fold, one recursive call per
 *  branch on the shared half-length plan, dense correction stage or stages, the arrowhead stage at odd
 *  lengths, and the interleaving permutation.
 *
 *  \version    1.1
 *  \date       August 2026
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright  (C) 2026 Yuriy A. Reznik.
 *  \license    MIT License (see the LICENSE file).
 */

#include <math.h>
#include <stdlib.h>
#include <stddef.h>

#include "ar1klt.h"

#ifndef SQRT2
#define SQRT2 1.41421356237309504880
#endif

/* ********************************************************************************************** */
/*  Plan node                                                                                     */
/* ********************************************************************************************** */

struct ar1klt_gen {
    int n;                    /* transform length at this node         */
    struct ar1klt_gen *half;  /* shared half-length plan (M = n/2)     */
    double *Qs;               /* even: M*M sum-branch stage            */
    double *Qa;               /* even: M*M difference-branch stage     */
    double *Qd;               /* odd:  M*M residual stage              */
    double *Aa;               /* odd:  (M+1)*(M+1) arrowhead stage     */
    int    *obr;              /* output row -> branch (0 sum, 1 diff)  */
    int    *oix;              /* output row -> index within branch     */
    double *geneig;           /* n generator eigenvalues, output order */
    double *pl;               /* n last components of rows (fold vec)  */
    double *buf;              /* scratch: 6*(M+1) doubles              */
    int     alt_ok;           /* branch spectra alternate strictly     */
};

/* ********************************************************************************************** */
/*  Secular solvers (design time)                                                                 */
/* ********************************************************************************************** */

/*!
 *  \brief      Evaluate the rank-one secular function f(nu) = 1 + sigma * sum_i p2_i / (d_i - nu).
 * 
 *  \param[in]  m       Number of poles.
 *  \param[in]  d       Poles d_1 < ... < d_m.
 *  \param[in]  p2      Squared components of the perturbation vector.
 *  \param[in]  sigma   Perturbation strength (of either sign).
 *  \param[in]  nu      Evaluation point, distinct from every pole.
 * 
 *  \return     The value of the secular function at nu.
 */
static double secfun(int m, const double *d, const double *p2, double sigma, double nu)
{
    int i;
    double s = 0.0;
    for (i = 0; i < m; i++) s += p2[i] / (d[i] - nu);
    return 1.0 + sigma * s;
}

/*!
 *  \brief      All roots of the rank-one secular equation f(nu) = 0 for the matrix diag(d) + sigma * p p^T.
 *
 *  \details
 *  The poles d must be given in ascending order, and every component of the perturbation vector must be
 *  nonzero, as the paper guarantees for the fold vector.  Each root is isolated in its interlacing
 *  interval, which lies one position above the corresponding pole when sigma is positive and one position
 *  below it when sigma is negative; the extreme interval is bounded by the total spectral shift
 *  sigma * ||p||^2.  Within each interval the secular function is strictly monotone, so plain bisection is
 *  safe; two hundred bisection steps drive the bracket to the resolution of double precision.
 *
 *  \param[in]  m       Number of poles (and of roots).
 *  \param[in]  d       Poles d_1 < ... < d_m, in ascending order.
 *  \param[in]  p2      Squared components of the perturbation vector.
 *  \param[in]  sigma   Perturbation strength (of either sign, nonzero).
 *  \param[out] r       The m roots, in ascending order, strictly interlacing the poles.
 */
static void secular_roots(int m, const double *d, const double *p2, double sigma, double *r)
{
    int i, it;
    double shift = 0.0;
    for (i = 0; i < m; i++) shift += p2[i];
    shift *= sigma;                       /* total spectral shift      */
    for (i = 0; i < m; i++) {
        double lo, hi, a, b, fs;
        if (sigma > 0.0) {                /* root i in (d_i, d_{i+1})  */
            lo = d[i];
            hi = (i < m - 1) ? d[i + 1] : d[m - 1] + shift;
            fs = -1.0;                    /* sign of f at lo+          */
        } else {                          /* root i in (d_{i-1}, d_i)  */
            hi = d[i];
            lo = (i > 0) ? d[i - 1] : d[0] + shift;
            fs = 1.0;
        }
        a = lo; b = hi;
        for (it = 0; it < 200; it++) {
            double mid = 0.5 * (a + b), fm;
            if (!(mid > a && mid < b)) break;
            fm = secfun(m, d, p2, sigma, mid);
            if ((fm >= 0.0 ? 1.0 : -1.0) == fs) a = mid; else b = mid;
        }
        r[i] = 0.5 * (a + b);
    }
}

/*!
 *  \brief      Evaluate the arrowhead secular function g(nu) = alpha - nu - sum_j b2_j / (dl_j - nu).
 * 
 *  \param[in]  m       Number of shaft entries.
 *  \param[in]  dl      Shaft diagonal dl_1 < ... < dl_m.
 *  \param[in]  b2      Squared border entries.
 *  \param[in]  alpha   Corner entry of the arrowhead matrix.
 *  \param[in]  nu      Evaluation point, distinct from every shaft entry.
 * 
 *  \return     The value of the arrowhead secular function at nu.
 */
static double arrfun(int m, const double *dl, const double *b2,
                     double alpha, double nu)
{
    int j;
    double s = 0.0;
    for (j = 0; j < m; j++) s += b2[j] / (dl[j] - nu);
    return alpha - nu - s;
}

/*!
 *  \brief      All eigenvalues of the arrowhead matrix [[diag(dl), b], [b^T, alpha]].
 *
 *  The m + 1 eigenvalues are the roots of the arrowhead secular equation g(nu) = 0 and strictly interlace
 *  the shaft diagonal from outside: nu_1 < dl_1 < nu_2 < ... < dl_m < nu_{m+1}.  The two extreme roots are
 *  bracketed by the norm of the border vector; the interior roots are bracketed by consecutive shaft
 *  entries.  The secular function is strictly decreasing within every bracket, so plain bisection is safe.
 *
 *  \param[in]  m       Number of shaft entries.
 *  \param[in]  dl      Shaft diagonal dl_1 < ... < dl_m, in ascending order.
 *  \param[in]  b2      Squared border entries; all must be nonzero.
 *  \param[in]  alpha   Corner entry of the arrowhead matrix.
 *  \param[out] r       The m + 1 eigenvalues, in ascending order.
 */
static void arrow_roots(int m, const double *dl, const double *b2, double alpha, double *r)
{
    int i, it;
    double nrm = 0.0;
    for (i = 0; i < m; i++) nrm += b2[i];
    nrm = sqrt(nrm) + 1e-12;
    for (i = 0; i <= m; i++) {
        double lo, hi, a, b;
        if (i == 0) {
            lo = (dl[0] < alpha ? dl[0] : alpha) - nrm;
            hi = dl[0];
        } else if (i == m) {
            lo = dl[m - 1];
            hi = (dl[m - 1] > alpha ? dl[m - 1] : alpha) + nrm;
        } else {
            lo = dl[i - 1];
            hi = dl[i];
        }
        a = lo; b = hi;                   /* g(lo+) > 0 > g(hi-)       */
        for (it = 0; it < 200; it++) {
            double mid = 0.5 * (a + b);
            if (!(mid > a && mid < b)) break;
            if (arrfun(m, dl, b2, alpha, mid) > 0.0) a = mid;
            else b = mid;
        }
        r[i] = 0.5 * (a + b);
    }
}

/* ********************************************************************************************** */
/*  Stage assembly (design time)                                                                  */
/* ********************************************************************************************** */

/*!
 *  \brief      Assemble a correction stage in its row-normalized Cauchy form: Q[i][j] ~ p_j / (mu_j - nu_i).
 *  \param[in]  m    Stage order.
 *  \param[in]  p    Fold-vector components of the half-length transform; all nonzero.
 *  \param[in]  mu   Generator eigenvalues of the half-length transform, in output (ascending) order.
 *  \param[in]  nu   Secular roots of the branch, in ascending order.
 *  \param[out] Q    Assembled orthogonal stage of order m, stored row-major.
 */
static void cauchy_stage(int m, const double *p, const double *mu,
                         const double *nu, double *Q)
{
    int i, j;
    for (i = 0; i < m; i++) {
        double s = 0.0;
        for (j = 0; j < m; j++) {
            double v = p[j] / (mu[j] - nu[i]);
            Q[i * m + j] = v;
            s += v * v;
        }
        s = 1.0 / sqrt(s);
        for (j = 0; j < m; j++) Q[i * m + j] *= s;
    }
}

/*!
 *  \brief      Merge the two ascending branch spectra into the global output order of the transform.
 *
 *  \details
 *  The merged order is ascending in the generator eigenvalue, which is the KLT order of decreasing
 *  variance.  For every output row the routine records which branch it comes from and its index within
 *  that branch, and it verifies that consecutive output rows always come from different branches.  This
 *  strict alternation is the parity-alternation property that justifies the interleaving permutation of
 *  the paper; it is checked here rather than assumed.
 *
 *  \param[in]  ns    Number of sum-branch eigenvalues.
 *  \param[in]  es    Sum-branch eigenvalues, in ascending order.
 *  \param[in]  na    Number of difference-branch eigenvalues.
 *  \param[in]  ea    Difference-branch eigenvalues, in ascending order.
 *  \param[out] obr   Per output row: 0 for the sum branch, 1 for the difference branch.
 *  \param[out] oix   Per output row: index of the row within its branch.
 *  \param[out] eig   Merged eigenvalues, in ascending order.
 *  \return     1 if the merged order alternates strictly between the branches, 0 otherwise.
 */
static int merge_branches(int ns, const double *es, int na,
                          const double *ea, int *obr, int *oix,
                          double *eig)
{
    int i = 0, j = 0, k = 0, alt = 1;
    while (i < ns || j < na) {
        int takes = (j >= na) || (i < ns && es[i] < ea[j]);
        if (takes) { obr[k] = 0; oix[k] = i; eig[k] = es[i]; i++; }
        else       { obr[k] = 1; oix[k] = j; eig[k] = ea[j]; j++; }
        if (k > 0 && obr[k] == obr[k - 1]) alt = 0;
        k++;
    }
    return alt;
}

/* ********************************************************************************************** */
/*  Construction                                                                                  */
/* ********************************************************************************************** */

/*!
 *  \brief      Release one recursion node together with its entire chain of half-length sub-plans.
 *  \param[in]  g   Node to release; NULL is allowed and ignored.
 */
static void gen_free_node(struct ar1klt_gen *g)
{
    if (g == NULL) return;
    gen_free_node(g->half);
    free(g->Qs); free(g->Qa); free(g->Qd); free(g->Aa);
    free(g->obr); free(g->oix); free(g->geneig); free(g->pl);
    free(g->buf);
    free(g);
}

/*!
 *  \brief      Release a plan created by ar1klt_gen_init(), including all of its half-length sub-plans.
 *  \param[in]  plan   Plan to release; NULL is allowed and ignored.
 */
void ar1klt_gen_free(ar1klt_gen *plan)
{
    gen_free_node(plan);
}

/*!
 *  \brief      Allocate one recursion node with all pointers cleared.
 *  \param[in]  n   Transform length stored at the node.
 *  \return     The newly allocated node, or NULL when memory is exhausted.
 */
static struct ar1klt_gen *node_alloc(int n)
{
    struct ar1klt_gen *g =
        (struct ar1klt_gen *)malloc(sizeof(struct ar1klt_gen));
    if (g == NULL) return NULL;
    g->n = n; g->half = NULL;
    g->Qs = g->Qa = g->Qd = g->Aa = NULL;
    g->obr = g->oix = NULL;
    g->geneig = g->pl = g->buf = NULL;
    g->alt_ok = 1;
    return g;
}

static int build(struct ar1klt_gen **out, int n, double rho);

/*!
 *  \brief      Recursively build the plan node for one transform length.
 *
 *  \details
 *  The routine first builds the shared half-length plan, then constructs the correction stages of the
 *  current level from the eigenvalue data and the fold vector of the half-length plan, merges the branch
 *  spectra into the output order, and finally computes its own fold vector by applying the freshly built
 *  transform to the last unit coordinate vector.  The lengths one and two are explicit base cases.
 *
 *  \param[out] out   Receives the completed node on success.
 *  \param[in]  n     Transform length of the node.
 *  \param[in]  rho   AR(1) correlation coefficient.
 *  \return     0 on success; -1 on allocation failure; -2 if the branch spectra failed to alternate.
 */
static int build(struct ar1klt_gen **out, int n, double rho)
{
    struct ar1klt_gen *g = node_alloc(n);
    int m, i, j, rc;
    double *mu, *p, *p2, *e1, *e2, *q, *b2;

    *out = NULL;
    if (g == NULL) return -1;
    g->geneig = (double *)malloc((size_t)n * sizeof(double));
    g->pl     = (double *)malloc((size_t)n * sizeof(double));
    if (g->geneig == NULL || g->pl == NULL) { gen_free_node(g); return -1; }

    if (n == 1) {
        g->geneig[0] = 1.0 - rho * rho;
        g->pl[0] = 1.0;
        *out = g;
        return 0;
    }
    if (n == 2) {
        g->geneig[0] = 1.0 - rho;
        g->geneig[1] = 1.0 + rho;
        g->pl[0] =  1.0 / SQRT2;
        g->pl[1] = -1.0 / SQRT2;
        *out = g;
        return 0;
    }

    m = n / 2;
    rc = build(&g->half, m, rho);
    if (rc != 0) { gen_free_node(g); return rc; }
    mu = g->half->geneig;
    p  = g->half->pl;

    g->obr = (int *)malloc((size_t)n * sizeof(int));
    g->oix = (int *)malloc((size_t)n * sizeof(int));
    g->buf = (double *)malloc((size_t)(6 * (m + 1)) * sizeof(double));
    p2 = (double *)malloc((size_t)m * sizeof(double));
    e1 = (double *)malloc((size_t)(m + 1) * sizeof(double));
    e2 = (double *)malloc((size_t)(m + 1) * sizeof(double));
    if (g->obr == NULL || g->oix == NULL || g->buf == NULL ||
        p2 == NULL || e1 == NULL || e2 == NULL) {
        free(p2); free(e1); free(e2); gen_free_node(g); return -1;
    }
    for (i = 0; i < m; i++) p2[i] = p[i] * p[i];

    if (n % 2 == 0) {
        g->Qs = (double *)malloc((size_t)(m * m) * sizeof(double));
        g->Qa = (double *)malloc((size_t)(m * m) * sizeof(double));
        if (g->Qs == NULL || g->Qa == NULL) {
            free(p2); free(e1); free(e2); gen_free_node(g); return -1;
        }
        secular_roots(m, mu, p2, -rho * (1.0 - rho), e1);
        secular_roots(m, mu, p2,  rho * (1.0 + rho), e2);
        cauchy_stage(m, p, mu, e1, g->Qs);
        cauchy_stage(m, p, mu, e2, g->Qa);
        g->alt_ok = merge_branches(m, e1, m, e2,
                                   g->obr, g->oix, g->geneig);
    } else {
        g->Qd = (double *)malloc((size_t)(m * m) * sizeof(double));
        g->Aa = (double *)malloc((size_t)((m + 1) * (m + 1))
                                 * sizeof(double));
        q  = (double *)malloc((size_t)m * sizeof(double));
        b2 = (double *)malloc((size_t)m * sizeof(double));
        if (g->Qd == NULL || g->Aa == NULL || q == NULL || b2 == NULL) {
            free(p2); free(e1); free(e2); free(q); free(b2);
            gen_free_node(g); return -1;
        }
        secular_roots(m, mu, p2, rho * rho, e2);      /* deltas       */
        cauchy_stage(m, p, mu, e2, g->Qd);
        for (i = 0; i < m; i++) {                     /* q = Qd p     */
            double s = 0.0;
            for (j = 0; j < m; j++) s += g->Qd[i * m + j] * p[j];
            q[i] = s;
            b2[i] = 2.0 * rho * rho * s * s;          /* (sqrt2 rho q)^2 */
        }
        arrow_roots(m, e2, b2, 1.0 + rho * rho, e1);  /* nus          */
        for (i = 0; i <= m; i++) {                    /* arrowhead rows */
            double s = 1.0;
            for (j = 0; j < m; j++) {
                double v = SQRT2 * rho * q[j] / (e2[j] - e1[i]);
                g->Aa[i * (m + 1) + j] = v;
                s += v * v;
            }
            g->Aa[i * (m + 1) + m] = 1.0;
            s = 1.0 / sqrt(s);
            for (j = 0; j <= m; j++) g->Aa[i * (m + 1) + j] *= s;
        }
        g->alt_ok = merge_branches(m + 1, e1, m, e2,
                                   g->obr, g->oix, g->geneig);
        free(q); free(b2);
    }
    free(p2); free(e1); free(e2);

    /* fold vector for the parent: pl = W_n e_{n-1}                   */
    {
        double *ein = (double *)malloc((size_t)n * sizeof(double));
        if (ein == NULL) { gen_free_node(g); return -1; }
        for (i = 0; i < n; i++) ein[i] = 0.0;
        ein[n - 1] = 1.0;
        ar1klt_gen_apply(g, ein, g->pl);
        free(ein);
    }

    if (!g->alt_ok) { gen_free_node(g); return -2; }
    *out = g;
    return 0;
}

/*!
 *  \brief      Build the recursive factorization plan for the length-n AR(1) KLT at a given rho.
 *  \param[out] plan   Receives the newly allocated plan; release it with ar1klt_gen_free().
 *  \param[in]  n      Transform length, n >= 1.
 *  \param[in]  rho    AR(1) correlation coefficient, 0 < rho < 1.
 *  \return     0 on success; -1 on invalid arguments or allocation failure; -2 if the branch spectra
 *              failed to alternate strictly (never observed; it would indicate a numerical breakdown).
 */
int ar1klt_gen_init(ar1klt_gen **plan, int n, double rho)
{
    if (plan == NULL || n < 1) return -1;
    if (!(rho > 0.0 && rho < 1.0)) return -1;
    return build(plan, n, rho);
}

/* ********************************************************************************************** */
/*  Application                                                                                   */
/* ********************************************************************************************** */

/*!
 *  \brief      Apply the recursive AR(1) KLT plan to one input block.
 *
 *  \details
 *  The routine folds the input by the butterfly, applies the shared half-length plan once to each branch,
 *  applies the dense correction stage or stages and, at odd lengths, the arrowhead stage that absorbs the
 *  center sample, and finally interleaves the branch outputs into the global KLT order.
 *
 *  \param[in]  plan   Plan built by ar1klt_gen_init().
 *  \param[in]  x      Input block of plan length; must not alias y.
 *  \param[out] y      KLT coefficients, ordered by decreasing variance.
 */
void ar1klt_gen_apply(const ar1klt_gen *plan, const double *x, double *y)
{
    int n = plan->n, m, i, j;
    double *u, *v, *gb, *hb, *zb, *tb;

    if (n == 1) { y[0] = x[0]; return; }
    if (n == 2) {
        double a = x[0], b = x[1];
        y[0] = (a + b) / SQRT2;
        y[1] = (a - b) / SQRT2;
        return;
    }
    m = n / 2;
    u  = plan->buf;
    v  = u + (m + 1);
    gb = v + (m + 1);
    hb = gb + (m + 1);
    zb = hb + (m + 1);
    tb = zb + (m + 1);

    for (i = 0; i < m; i++) {
        u[i] = (x[i] + x[n - 1 - i]) / SQRT2;
        v[i] = (x[i] - x[n - 1 - i]) / SQRT2;
    }
    ar1klt_gen_apply(plan->half, u, gb);   /* shared half-plan, twice */
    ar1klt_gen_apply(plan->half, v, hb);

    if (n % 2 == 0) {
        for (i = 0; i < m; i++) {          /* zb = Qs gb, tb = Qa hb  */
            double s = 0.0, t = 0.0;
            for (j = 0; j < m; j++) {
                s += plan->Qs[i * m + j] * gb[j];
                t += plan->Qa[i * m + j] * hb[j];
            }
            zb[i] = s; tb[i] = t;
        }
    } else {
        double c = x[m];
        for (i = 0; i < m; i++) {          /* residual stage, both    */
            double s = 0.0, t = 0.0;
            for (j = 0; j < m; j++) {
                s += plan->Qd[i * m + j] * gb[j];
                t += plan->Qd[i * m + j] * hb[j];
            }
            u[i] = s; tb[i] = t;           /* reuse u for r            */
        }
        u[m] = c;
        for (i = 0; i <= m; i++) {         /* arrowhead: zb = Aa (r,c) */
            double s = 0.0;
            for (j = 0; j <= m; j++)
                s += plan->Aa[i * (m + 1) + j] * u[j];
            zb[i] = s;
        }
    }
    for (i = 0; i < n; i++)                /* interleave               */
        y[i] = (plan->obr[i] == 0) ? zb[plan->oix[i]] : tb[plan->oix[i]];
}
