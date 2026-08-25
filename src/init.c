#include <R_ext/RS.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* Fortran MCMC kernels.
 *   qrb_al_mcmc  (38 args) : adaptive-lasso binary quantile regression
 *                           (34 + the four derivation-version parameters
 *                            alpha_om, beta_om, b0_mean and b0_prec, which
 *                            select the V3 / V4 / V5 hierarchy; see the
 *                            declaration block in QRb_AL_mcmc.f95)
 *   qrb_l_mcmc   (35 args) : lasso binary quantile regression
 *                           (31 + alpha_t, beta_t, b0_mean, b0_prec: the
 *                            proper hyperprior on tau_h and the proper
 *                            intercept prior)
 *   qrb_bqr_mcmc (19 args) : unpenalised binary quantile regression
 *                           (17 + b0_mean, b0_prec)
 *
 * Continuous-response counterparts (derivations/Mathematical_Derivations_V5C.tex).
 * Each drops the y* draw and the five identification-anchor switches, always
 * samples sigma, and gains q_exp (the exponent on sigma in the local-scale
 * prior) and sigma_init:
 *   qrc_al_mcmc  (34 args) : adaptive-lasso continuous quantile regression
 *                           (NCLIP is 16 here, not 18: C_YSTAR and C_V4SD
 *                            guarded the deleted y* draw)
 *   qrc_l_mcmc   (32 args) : lasso continuous quantile regression
 *   qrc_bqr_mcmc (20 args) : unpenalised continuous quantile regression
 *                           (GAINS a, b, q_exp and sigma_init over its binary
 *                            sibling, which never sampled sigma)
 * Argument counts are verified against the subroutine declarations in
 * QRb_AL_mcmc.f95, QRb_L_mcmc.f95 and QRb_BQR_mcmc.f95 respectively. A
 * mismatch here corrupts memory rather than erroring, so keep them in step.
 */

#define V void *

extern void F77_NAME(qrb_al_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V, V, V, V);

extern void F77_NAME(qrb_l_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V);

extern void F77_NAME(qrb_bqr_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V);

extern void F77_NAME(qrc_al_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V);

extern void F77_NAME(qrc_l_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V);

extern void F77_NAME(qrc_bqr_mcmc)(
    V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V, V,
    V, V, V);

#undef V

static const R_FortranMethodDef FortranEntries[] = {
    {"qrb_al_mcmc",  (DL_FUNC) &F77_NAME(qrb_al_mcmc),  38},
    {"qrb_l_mcmc",   (DL_FUNC) &F77_NAME(qrb_l_mcmc),   35},
    {"qrb_bqr_mcmc", (DL_FUNC) &F77_NAME(qrb_bqr_mcmc), 19},
    {"qrc_al_mcmc",  (DL_FUNC) &F77_NAME(qrc_al_mcmc),  34},
    {"qrc_l_mcmc",   (DL_FUNC) &F77_NAME(qrc_l_mcmc),   32},
    {"qrc_bqr_mcmc", (DL_FUNC) &F77_NAME(qrc_bqr_mcmc), 20},
    {NULL, NULL, 0}
};

void R_init_bbqr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
