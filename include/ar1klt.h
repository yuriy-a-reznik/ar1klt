/*!
 *  \file ar1klt.h
 *  \brief AR1KLT: exact fast factorizations of the AR(1) Karhunen-Loeve transform (KLT) for short lengths N = 2..8.
 *
 *  \details
 *  This library implements the exact short-length KLT algorithms of
 *
 *   Y. A. Reznik, "Direct Factorization of the Karhunen-Loeve Transform of AR(1) Sources", Section 5.
 *
 *  The transforms diagonalize the covariance matrix R_N(rho) = [rho^|i-j|] of a stationary AR(1) source 
 *  with correlation coefficient rho in (0,1).
 *
 *  Each transform is realized as a short chain of butterflies and plane rotations followed by small 
 *  dense "secular" stages.  All constants are algebraic functions of rho (elementary or radical expressions; 
 *  roots of cubic and quartic characteristic equations are obtained by the Cardano-Ferrari formulas).  
 *  Constants are precomputed once by ar1klt_init() and reused by ar1klt_apply().
 *
 *  Output coefficients y_0, ..., y_{N-1} appear in KLT order (decreasing variance). 
 *  Rows of the transform are defined up to sign, as usual for eigenvector bases.
 *
 *  The code is strict ANSI C (C89) and has no dependencies beyond libc/libm.
 *
 *  \author    Yuriy A. Reznik <yreznik@mit.edu>
 *  \copyright Copyright (C) 2026 Yuriy A. Reznik
 *  \license   MIT License (see LICENSE file).
 */

#ifndef AR1KLT_H
#define AR1KLT_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Smallest supported transform length. */
#define AR1KLT_NMIN 2
/*! \brief Largest supported transform length. */
#define AR1KLT_NMAX 8

/*!
 * \brief Precomputed constants ("plan") for one AR(1) KLT of a given
 *        length and correlation coefficient.
 *
 * \details
 * The structure is filled by ar1klt_init() and must be treated as opaque
 * by callers; it holds the rotation angles and the small dense stage
 * matrices of the corresponding Section-5 module.
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
 * \brief Precompute all constants of the length-n AR(1) KLT at a given rho.
 *
 * \details
 * Evaluates the closed-form angles and the roots of the cubic/quartic
 * characteristic equations of the Section-5 module for the requested
 * length, and assembles the dense secular and arrowhead stage matrices.
 * The call is intended for design time; its cost is negligible and
 * independent of any data.
 *
 * \param[out] ctx  Plan to fill; must be non-NULL.
 * \param[in]  n    Transform length, AR1KLT_NMIN <= n <= AR1KLT_NMAX.
 * \param[in]  rho  AR(1) correlation coefficient, 0 < rho < 1.
 * \return 0 on success; -1 if ctx is NULL or an argument is out of range.
 */
int ar1klt_init(ar1klt_ctx *ctx, int n, double rho);

/*!
 * \brief Apply the precomputed AR(1) KLT to one input block.
 *
 * \details
 * Executes the stage chain of the corresponding Section-5 module
 * (butterflies, plane rotations, and the small dense secular stages)
 * exactly as printed in the paper.  The input and output arrays must
 * each hold ctx->n doubles and may not alias.
 *
 * \param[in]  ctx  Plan previously filled by ar1klt_init().
 * \param[in]  x    Input block x_0, ..., x_{n-1}.
 * \param[out] y    Output KLT coefficients y_0, ..., y_{n-1},
 *                  ordered by decreasing variance.
 */
void ar1klt_apply(const ar1klt_ctx *ctx, const double *x, double *y);

#ifdef __cplusplus
}
#endif

#endif /* AR1KLT_H */
