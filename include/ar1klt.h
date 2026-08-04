/*!
 *  \file       ar1klt.h
 *  \brief      AR1KLT: exact fast algorithms for the Karhunen-Loeve transform (KLT) of AR(1) sources.
 *
 *  \details
 *  This library implements the AR(1) KLT factorizations of the paper
 *
 *      Y. Reznik, "Direct Factorization of the Karhunen-Loeve Transform of AR(1) Sources",
 *      submitted to the IEEE Transactions on Signal Processing, July 2026.
 *
 *  The transforms diagonalize the covariance matrix R_N(rho) = [rho^|i-j|] of a stationary first-order
 *  autoregressive (AR(1)) source with correlation coefficient rho in (0,1). 
 *
 *  Two interfaces are provided:
 *
 *  The first interface (ar1klt_init / ar1klt_apply) offers fixed short-length modules for N = 2..8,
 *  corresponding to Section V of the paper.  Each module is a short chain of butterflies and plane rotations
 *  followed by small dense stages.  All constants are algebraic functions of rho: rotation angles are
 *  closed-form arctangents, and the remaining constants are roots of cubic and quartic characteristic
 *  equations, evaluated by the Cardano and Ferrari radical formulas and polished by Newton iterations.
 *
 *  The second interface (ar1klt_gen_init / ar1klt_gen_apply / ar1klt_gen_free) implements the recursive
 *  self-similar factorizations of Theorems III.3 and IV.4 of the paper for arbitrary length N >= 1.
 *  At even N the transform is a butterfly, two applications of the half-length transform, and two dense
 *  correction stages; at odd N it is a butterfly, two applications of the half-length transform followed by
 *  a shared residual correction stage, and one arrowhead stage that absorbs the center sample.
 *  One half-length plan is built once and shared by both branches, so the recursion is genuine.
 *
 *  In both interfaces the output coefficients y_0, ..., y_{N-1} appear in KLT order (decreasing variance),
 *  and the rows of the transform are defined up to sign, as is usual for eigenvector bases.
 *
 *  The code is strict ANSI C (C89) and has no dependencies beyond the standard C library and libm.
 *
 *  \version    1.1
 *  \date       August 2026
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright  (C) 2026 Yuriy A. Reznik.
 *  \license    MIT License (see the LICENSE file).
 */

#ifndef _AR1KLT_H_
#define _AR1KLT_H_

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Smallest supported transform length. */
#define AR1KLT_NMIN 2
/*! \brief Largest supported transform length. */
#define AR1KLT_NMAX 8

/*!
 *  \brief Precomputed constants for one AR(1) KLT of a given length and correlation coefficient.
 */
typedef struct ar1klt_ctx {
    int    n;          /*!< Transform length, AR1KLT_NMIN..AR1KLT_NMAX.   */
    double rho;        /*!< AR(1) correlation coefficient, 0 < rho < 1.   */
    double phi3;       /*!< 3-point rotation angle (N = 3, 6, 7).         */
    double phi_s;      /*!< 4-point sum-branch angle (N = 4, 8).          */
    double phi_a;      /*!< 4-point difference-branch angle (N = 4, 8).   */
    double psi;        /*!< Residual rotation angle (N = 5).              */
    double Qs[4][4];   /*!< Sum-branch secular stage (N = 6: 3x3; 8: 4x4).*/
    double Qa[4][4];   /*!< Diff-branch secular stage (N = 6: 3x3; 8: 4x4).*/
    double Qd[4][4];   /*!< Residual secular stage (N = 7: 3x3).          */
    double Aar[4][5];  /*!< Arrowhead stage (N = 5: 3x3; N = 7: 4x4).     */
} ar1klt_ctx;

/*!
 *  \brief Precompute all constants of the length-n AR(1) KLT at a given rho.
 *
 *  \details
 *  Evaluates the closed-form angles and the roots of the cubic/quartic characteristic equations
 *  of the Section-5 module for the requested length, and assembles the dense secular and arrowhead 
 *  stage matrices. The call is intended for design time; its cost is negligible and
 *  independent of any data.
 *
 *  \param[out] ctx  Plan to fill; must be non-NULL.
 *  \param[in]  n    Transform length, AR1KLT_NMIN <= n <= AR1KLT_NMAX.
 *  \param[in]  rho  AR(1) correlation coefficient, 0 < rho < 1.
 *  \return 0 on success; -1 if ctx is NULL or an argument is out of range.
 */
int ar1klt_init(ar1klt_ctx *ctx, int n, double rho);

/*!
 *  \brief Apply the precomputed AR(1) KLT to one input block.
 *
 *  \details
 *  Executes the stage chain of the corresponding Section-5 module (butterflies, plane rotations,
 *  and the small dense secular stages) exactly as printed in the paper. 
 *  The input and output arrays must each hold ctx->n doubles and may not alias.
 *
 *  \param[in]  ctx  Plan previously filled by ar1klt_init().
 *  \param[in]  x    Input block x_0, ..., x_{n-1}.
 *  \param[out] y    Output KLT coefficients y_0, ..., y_{n-1}, ordered by decreasing variance.
 */
void ar1klt_apply(const ar1klt_ctx *ctx, const double *x, double *y);

/* ********************************************************************************************** */
/*  General-length recursive factorization (any N >= 1)                                           */
/* ********************************************************************************************** */

/*!
 *  \brief Opaque recursive constants set for the AR(1) KLT of arbitrary length.
 */
typedef struct ar1klt_gen ar1klt_gen;

/*!
 *  \brief Build the recursive constants set for the length-n AR(1) KLT at rho.
 *
 *  \param[out] plan  Receives the newly allocated plan; free with ar1klt_gen_free().
 *  \param[in]  n     Transform length, n >= 1.
 *  \param[in]  rho   AR(1) correlation coefficient, 0 < rho < 1.
 * 
 *  \return 0 on success; 
 *          -1 on bad arguments or allocation failure;
 *          -2 if the branch spectra failed to alternate (indicates a numerical breakdown).
 */
int ar1klt_gen_init(ar1klt_gen **plan, int n, double rho);

/*!
 *  \brief Apply the recursive AR(1) KLT plan to one input block.
 *
 *  \param[in]  plan  Plan built by ar1klt_gen_init().
 *  \param[in]  x     Input block of plan length.
 *  \param[out] y     KLT coefficients, ordered by decreasing variance.
 */
void ar1klt_gen_apply(const ar1klt_gen *plan, const double *x, double *y);

/*!
 *  \brief Release a plan created by ar1klt_gen_init().
 * 
 *  \param[in] plan  Plan to free (NULL is allowed).
 */
void ar1klt_gen_free(ar1klt_gen *plan);

#ifdef __cplusplus
}
#endif

#endif /* AR1KLT_H */
