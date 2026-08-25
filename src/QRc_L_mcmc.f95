! =====================================================================
! Continuous quantile regression with the (standard) LASSO
! =====================================================================
! Continuous-response counterpart of QRb_L_mcmc.f95; implements the lasso
! section of derivations/Mathematical_Derivations_V5C.tex.
!
! THE LASSO LAYER IS A DIFFERENT HIERARCHY FROM THE ADAPTIVE ONE, not the
! adaptive one with k = 1.  It descends from Benoit, Alhamzawi & Yu (2013):
! a single global lambda^2 parameterised as a RATE with a Gamma prior and a
! second-level hyperparameter tau_h, where QRc_AL_mcmc uses per-coefficient
! lambda_j^2 as a SCALE with an inverse-Gamma prior and the hyperparameter
! omega.  Reading one off the other gives a model neither appendix describes.
!
!   s_j     | sigma, lambda^2 ~ Exp(rate = sigma^q * lambda^2 / 2)
!   lambda^2| delta, tau_h    ~ Gamma(delta, tau_h)
!   tau_h                     ~ Gamma(alpha_t, beta_t)   [V5 correction]
!   delta                     ~ Gamma(alpha_d, beta_d)   [V5 correction]
!   sigma                     ~ Gamma(a, b)
!   beta0                     ~ N(b0_mean, 1/b0_prec)    [V5 correction]
!
! alpha_t = beta_t = 0 recovers the published p(tau_h) prop 1/tau_h, whose
! joint posterior does not exist -- the same wedge that kills p(omega) prop
! 1/omega in the adaptive layer.
!
! *** DEFAULT q IS 0 FOR THIS LAYER, NOT 1. ***
! Benoit et al.'s penalty is lambda*|beta_j|, which is sigma-free, i.e. q = 0;
! Alhamzawi, Yu & Benoit's adaptive penalty is (sigma^(1/2)/lambda_j)*|beta_j|,
! i.e. q = 1.  The binary kernels reproduce each faithfully, so the two layers
! of this package have always sat at different exponents -- invisibly, because
! the binary model anchors sigma.  Each continuous kernel therefore defaults to
! its own binary sibling's exponent, so that a paired binary/continuous
! comparison isolates the response and nothing else.  See Remark 5C.
!
! At q = 0 the s_j prior is free of sigma, so it contributes nothing to the
! sigma conditional: the shape is a + 3n/2 with no +k and the rate is built
! from the likelihood alone, exactly as in QRb_L_mcmc.
!
! Differences from QRb_L_mcmc, all of them model differences:
!   - y is observed and real; the y* truncated-normal draw and rtnorm_geweke
!     are DELETED.
!   - sigma is identified and is drawn every sweep; fix_sigma, sigma_fixed,
!     constrain_beta_norm, fix_beta1 and norm_slopes_only are absent.
!   - q_exp and sigma_init are new arguments.
subroutine QRc_L_mcmc(n, k, r, keep, &
     y, tau, X, &
     a, b, q_exp, sigma_init, &
     delta_init, tau_init, mh_sd0, target_acc, burn_in, use_boundaries, &
     alpha_d, beta_d, &
     alpha_t, beta_t, b0_mean, b0_prec, &
     betadraw, beta0draw, sigmadraw, lambda2draw, tauhyperdraw, deltadraw, accdraw, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- inputs --------
  integer,  intent(in) :: n, k, r, keep
  real(dp), intent(in) :: y(n)              ! OBSERVED continuous response
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  ! Gamma(shape = alpha_t, rate = beta_t) hyperprior on tau_h.  The published
  ! specification is p(tau_h) prop 1/tau_h = the point alpha_t = beta_t = 0,
  ! which reproduces the paper's Gamma(delta, lambda^2) conditional exactly and
  ! leaves the joint posterior improper.
  real(dp), intent(in) :: alpha_t, beta_t
  real(dp), intent(in) :: b0_mean, b0_prec
  real(dp), intent(in) :: a, b              ! Gamma(shape=a, rate=b) for sigma
  ! Exponent q >= 0 on sigma in the s_j rate.  q = 0 is Benoit et al. (2013) and
  ! the default; q = 2 is response-scale equivariant (Proposition 5C).
  real(dp), intent(in) :: q_exp
  real(dp), intent(in) :: sigma_init
  real(dp), intent(in) :: delta_init, tau_init
  real(dp), intent(in) :: mh_sd0, target_acc
  integer,  intent(in) :: burn_in
  logical,  intent(in) :: use_boundaries
  real(dp), intent(in) :: alpha_d, beta_d
  real(dp), intent(in) :: beta_init(k), beta0_init

  ! -------- outputs (thinned by `keep`) --------
  real(dp), intent(out) :: betadraw(r/keep, k)
  real(dp), intent(out) :: beta0draw(r/keep)
  real(dp), intent(out) :: sigmadraw(r/keep)
  real(dp), intent(out) :: lambda2draw(r/keep)
  real(dp), intent(out) :: tauhyperdraw(r/keep)
  real(dp), intent(out) :: deltadraw(r/keep)
  real(dp), intent(out) :: accdraw(r/keep)

  ! -------- locals --------
  integer :: i, j, it, outix, rk
  real(dp) :: theta, phisq, denom, num, m0, s0sq, prec0
  real(dp) :: sigma, lambda2, tau_h, delta, mh_sd, acceptance_rate
  integer  :: accepted
  real(dp) :: zstd, prop_std, log_curr, log_prop, delta_prop
  real(dp) :: lp_curr, lp_prop, log_acc, u
  real(dp) :: gamma_m, scale_g, sum_s
  real(dp) :: beta0, beta_old_j, var_j, mu_j
  real(dp) :: chi, psi, muIG
  real(dp) :: psi_all, chi_i
  real(dp) :: sigq, A1, A2, shape_sig
  real(dp), allocatable :: beta(:), s(:)
  real(dp), allocatable :: z(:), Xbeta(:), resid(:), resid2(:)
  real(dp), parameter :: one=1.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp
  real(dp), parameter :: s_min=1.0e-4_dp, s_max=1.0e4_dp
  real(dp), parameter :: mh_min=1.0e-4_dp, mh_max=2.0_dp

  if (n <= 0 .or. k <= 0) return
  if (r <= 0 .or. keep <= 0) return
  rk = r/keep; if (rk*keep /= r) return
  if (q_exp < 0.0_dp) return

  call rndstart()

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))

  allocate(beta(k), s(k))
  allocate(z(n), Xbeta(n), resid(n), resid2(n))

  beta  = beta_init
  beta0 = beta0_init
  z = 1.0_dp; s = 1.0_dp
  lambda2 = 1.0_dp
  sigma = sigma_init
  if (sigma <= 0.0_dp) sigma = one
  delta = delta_init; tau_h = tau_init

  mh_sd = mh_sd0; accepted = 0; acceptance_rate = 0.0_dp

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! ---- the y* draw of the binary sweep IS DELETED. ----

    ! (2) beta0 | .
    denom = 0.0_dp; num = 0.0_dp
    do i=1,n
      denom = denom + one / z(i)
      num   = num   + (y(i) - Xbeta(i) - theta*z(i)) / z(i)
    end do
    if (b0_prec > 0.0_dp) then
      prec0 = b0_prec + (sigma / phisq) * denom
      s0sq  = one / prec0
      m0    = s0sq * ( b0_prec*b0_mean + (sigma / phisq) * num )
    else
      m0   = num / denom
      s0sq = (phisq / sigma) / denom
    end if
    call rnorm(zstd); beta0 = m0 + sqrt(s0sq) * zstd

    ! (3) z_i | .
    psi_all = sigma * (theta*theta + 2.0_dp*phisq) / phisq
    do i=1,n
      resid(i) = y(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      call rgig_half(psi_all, chi_i, z(i))
      if (use_boundaries) z(i) = max(eps_small, min(z(i), big))
    end do

    ! (4) beta_j | .
    do j=1,k
      beta_old_j = beta(j)
      muIG = 0.0_dp
      do i=1,n
        muIG = muIG + (X(i,j)*X(i,j)) / z(i)
      end do
      denom = muIG + (phisq / sigma) / s(j)
      num = 0.0_dp
      do i=1,n
        num = num + X(i,j) * ( y(i) - beta0 - theta*z(i) - (Xbeta(i) - X(i,j)*beta_old_j) ) / z(i)
      end do
      mu_j  = num / denom
      var_j = (phisq / sigma) / denom
      call rnorm(zstd); beta(j) = mu_j + sqrt(var_j) * zstd
      do i=1,n
        Xbeta(i) = Xbeta(i) + X(i,j) * (beta(j) - beta_old_j)
      end do
    end do

    ! ---- the norm-projection block of the binary sweep IS DELETED. ----

    ! (5) s_j | . ~ GIG(1/2, psi = sigma^q * lambda^2, chi = beta_j^2)
    !     *** q ENTERS (1 of 3) ***  At q = 0 this is the published sigma-free
    !     psi = lambda^2 of QRb_L_mcmc.
    sigq = sigma**q_exp
    do j = 1, k
      chi   = beta(j)*beta(j)
      psi   = sigq * lambda2
      call rgig_half(psi, chi, s(j))
      if (use_boundaries) s(j) = max(s_min, min(s(j), s_max))
    end do

    ! (6) sigma | . -- ALWAYS EXECUTED.
    !     *** q ENTERS (2 of 3) ***
    !     shape = a + q*k + 3n/2 ; rate carries A2 = lambda^2 * sum(s_j) / 2
    !     multiplying sigma^q.  At q = 0 both extra terms vanish and this is
    !     exactly QRb_L_mcmc's Gamma(a + 3n/2, likelihood-only + b).
    do i = 1, n
      resid2(i) = y(i) - (beta0 + Xbeta(i)) - theta*z(i)
    end do
    A1 = b
    do i = 1, n
      A1 = A1 + ( (resid2(i)*resid2(i)) / (2.0_dp*phisq*z(i)) + z(i) )
    end do
    sum_s = 0.0_dp
    do j = 1, k
      sum_s = sum_s + s(j)
    end do
    A2 = 0.5_dp * lambda2 * sum_s
    shape_sig = a + q_exp*real(k,dp) + 1.5_dp*real(n,dp)
    if (q_exp <= 0.0_dp) then
      call rgamma(shape_sig, 1.0_dp/A1, sigma)
    else if (abs(q_exp - one) < 1.0e-12_dp) then
      call rgamma(shape_sig, 1.0_dp/(A1 + A2), sigma)
    else
      call slice_sigma(sigma, shape_sig, A1, A2, q_exp)
    end if
    if (use_boundaries) sigma = max(eps_small, min(sigma, big))

    ! (7) lambda^2 | . ~ Gamma(k + delta, sigma^q * sum(s_j)/2 + tau_h)
    !     *** q ENTERS (3 of 3) ***
    sigq = sigma**q_exp
    sum_s = 0.0_dp
    do j = 1, k
      sum_s = sum_s + s(j)
    end do
    muIG = real(k,dp) + delta
    scale_g = 0.5_dp*sigq*sum_s + tau_h
    call rgamma(muIG, 1.0_dp/scale_g, lambda2)
    if (use_boundaries) lambda2 = max(s_min, min(lambda2, s_max))

    ! (8) tau_h | . ~ Gamma(alpha_t + delta, beta_t + lambda^2).  Free of q.
    call rgamma(alpha_t + delta, 1.0_dp/(beta_t + lambda2), tau_h)
    if (use_boundaries) tau_h = max(1.0e-3_dp, min(tau_h, 1.0e3_dp))

    ! (9) delta | .  RW-MH on log delta.  Free of q.
    log_curr = log(delta)
    call rnorm(prop_std)
    log_prop  = log_curr + mh_sd * prop_std
    delta_prop = exp(log_prop)
    lp_curr = logpost_delta_L(delta,      tau_h, lambda2, alpha_d, beta_d)
    lp_prop = logpost_delta_L(delta_prop, tau_h, lambda2, alpha_d, beta_d)
    log_acc = (lp_prop + log(delta_prop)) - (lp_curr + log(delta))
    call runif(u)
    if (log(u) < log_acc) then
      delta = delta_prop
      if (use_boundaries) delta = max(1.0e-3_dp, min(delta, 1.0e3_dp))
      accepted = accepted + 1
    end if
    if (it <= burn_in) then
      gamma_m = 1.0_dp / sqrt( real(it,dp) )
      mh_sd = exp( log(mh_sd) + gamma_m * ( merge(1.0_dp, 0.0_dp, log(u) < log_acc ) - target_acc ) )
      if (use_boundaries) mh_sd = max(mh_min, min(mh_sd, mh_max))
    end if
    acceptance_rate = real(accepted,dp) / real(it,dp)

    if (mod(it, keep) == 0) then
      outix = outix + 1
      betadraw(outix, :)  = beta
      beta0draw(outix)    = beta0
      sigmadraw(outix)    = sigma
      lambda2draw(outix)  = lambda2
      tauhyperdraw(outix) = tau_h
      deltadraw(outix)    = delta
      accdraw(outix)      = acceptance_rate
    end if

    if (mod(it, 100000) == 0) then
      call intpr('Current iteration :', -1, it, 1)
    end if

  end do

  call rndend()

  deallocate(beta, s, z, Xbeta, resid, resid2)

contains

  ! Log of the unnormalised sigma conditional up to an additive constant.
  pure real(dp) function lt_sigma(x, sh, c1, c2, qq) result(lv)
    real(dp), intent(in) :: x, sh, c1, c2, qq
    if (x <= 0.0_dp) then
      lv = -1.0d300
    else
      lv = (sh - 1.0_dp)*log(x) - c1*x - c2*(x**qq)
    end if
  end function lt_sigma

  ! Stepping-out / shrinkage slice sampler; see QRc_AL_mcmc.f95 for the
  ! log-concavity argument that makes this exact and tuning-free.
  subroutine slice_sigma(x, sh, c1, c2, qq)
    real(dp), intent(inout) :: x
    real(dp), intent(in)    :: sh, c1, c2, qq
    real(dp) :: lvy, w, lo, hi, uu, xnew
    integer  :: istep
    ! MAXSTEP is a runaway backstop, not a tuning constant.  TRUNCATED
    ! stepping-out is reversible only with Neal's (2003) randomised allocation
    ! of steps between the two directions; unbounded stepping-out is reversible
    ! as it stands.  With w set to the target's own standard deviation below,
    ! 300 steps reaches 300 sd from the mode, so the cap cannot bind on a
    ! log-concave target and the procedure is the unbounded one in practice.
    integer, parameter :: MAXSTEP = 300, MAXSHRINK = 200
    real(dp), parameter :: tiny_sig = 1.0e-12_dp

    if (x <= 0.0_dp) x = 1.0_dp
    call runif(uu)
    if (uu <= 0.0_dp) return
    lvy = lt_sigma(x, sh, c1, c2, qq) + log(uu)

    ! Initial interval width.  This MUST NOT depend on the current state x:
    ! the stepping-out / shrinkage procedure is reversible only for a width
    ! chosen independently of where the chain currently sits (Neal 2003).  An
    ! earlier version used w = max(x, eps_small), which violates that and was
    ! replaced on the argument alone -- the bias it induces was NOT detectable
    ! here (see the validation note below), but the argument does not depend on
    ! detecting it.  sqrt(sh)/c1 is the standard deviation of the
    ! Gamma(sh, rate = c1) the target reduces to at c2 = 0, so it is correctly
    ! scaled and depends only on the conditioning quantities, which are fixed
    ! for this draw.
    !
    ! VALIDATION.  At q = 1 the conjugate branch above and this slice branch
    ! target the same law, so passing q = 1 + 1e-9 forces the slice path on a
    ! target the Gamma path samples exactly.  Five chains per branch,
    ! n = 300, k = 4, 40k sweeps: Gamma mean 1.66340 (sd across seeds 0.00156),
    ! slice mean 1.66304 (0.00118); Welch t = 0.41, p = 0.69.  A single pair of
    ! chains is NOT enough to check this -- one such pair differed by 0.28% with
    ! a Kolmogorov-Smirnov p of 0.001, which replication showed to be seed
    ! variation.  Re-run the five-seed control, not a single pair, after any
    ! change to this routine.
    w = sqrt(sh) / max(c1, eps_small)
    if (.not. (w > 0.0_dp)) w = eps_small
    call runif(uu)
    lo = max(x - w*uu, tiny_sig)
    hi = lo + w

    istep = 0
    do while (lt_sigma(lo, sh, c1, c2, qq) > lvy .and. lo > tiny_sig &
              .and. istep < MAXSTEP)
      lo = max(lo - w, tiny_sig)
      istep = istep + 1
    end do
    istep = 0
    do while (lt_sigma(hi, sh, c1, c2, qq) > lvy .and. istep < MAXSTEP)
      hi = hi + w
      istep = istep + 1
    end do

    do istep = 1, MAXSHRINK
      call runif(uu)
      xnew = lo + uu*(hi - lo)
      ! Shrink, never expand.  An earlier version assigned lo = xnew here,
      ! which WIDENS the interval when xnew is below the floor and is not a
      ! shrinkage step at all.
      if (xnew <= tiny_sig) then
        lo = tiny_sig
        cycle
      end if
      if (lt_sigma(xnew, sh, c1, c2, qq) >= lvy) then
        x = xnew
        return
      end if
      if (xnew < x) then
        lo = xnew
      else
        hi = xnew
      end if
    end do
  end subroutine slice_sigma

  ! ---------------------------------------------------------------------------
  ! COPIED VERBATIM from QRb_L_mcmc.f95 (lines 358-413 and 474-509), minus
  ! rtnorm_geweke, whose only caller was the y* draw.  Keep byte-identical.
  ! ---------------------------------------------------------------------------

  pure real(dp) function logpost_delta_L(d, tau_h, lam2, alpha_d, beta_d) result(lp)
    implicit none
    real(dp), intent(in) :: d, tau_h, lam2, alpha_d, beta_d
    if (.not.(d > 0.0_dp .and. tau_h > 0.0_dp .and. lam2 > 0.0_dp)) then
      lp = -1.0d300; return
    end if
    lp = d*log(tau_h*lam2) - log_gamma(d) + (alpha_d - 1.0_dp)*log(d) - beta_d*d
  end function logpost_delta_L

  subroutine runif(fn_val)
    real(dp), intent(out) :: fn_val
    real(dp) :: unifrnd
    fn_val = unifrnd()
  end subroutine

  subroutine rnorm(fn_val)
    real(dp), intent(out) :: fn_val
    real(dp) :: normrnd
    fn_val = normrnd()
  end subroutine

  ! Draw GIG(1/2, psi, chi). For chi > 0, its reciprocal is
  ! inverse-Gaussian with mean sqrt(psi/chi) and shape psi. The expression for
  ! q is rationalized to avoid the cancellation present in the usual formula.
  ! At chi = 0 the limiting law is Gamma(1/2, rate=psi/2).
  subroutine rgig_half(psi, chi, fn_val)
    implicit none
    real(dp), intent(in)  :: psi, chi
    real(dp), intent(out) :: fn_val
    real(dp) :: mu, lambda, nu, q, u, t
    if (chi == 0.0_dp) then
      call rgamma(0.5_dp, 2.0_dp/psi, fn_val)
      return
    end if
    mu = sqrt(psi/chi)
    lambda = psi
    call rnorm(nu)
    nu = nu*nu
    if (nu == 0.0_dp) then
      q = mu
    else
      t = mu*nu/lambda
      if (t < sqrt(huge(one))) then
        q = 2.0_dp*mu / (2.0_dp + t + sqrt(t)*sqrt(4.0_dp + t))
      else
        q = (2.0_dp*lambda/nu) / &
            (1.0_dp + sqrt(1.0_dp + 4.0_dp/t) + 2.0_dp/t)
      end if
    end if
    call runif(u)
    if (u <= 1.0_dp/(1.0_dp + q/mu)) then
      fn_val = 1.0_dp/q
    else
      fn_val = (q/mu)/mu
    end if
  end subroutine rgig_half

  subroutine rgamma(shape, scale, fn_val)
    implicit none
    real(dp), intent(in)  :: shape, scale
    real(dp), intent(out) :: fn_val
    real(dp) :: a, d, c, x, v, u
    logical :: flag
    if (shape < 1.0_dp) then
      a = shape + 1.0_dp
    else
      a = shape
    end if
    d = a - 1.0_dp/3.0_dp
    c = 1.0_dp / sqrt(9.0_dp*d)
    flag = .true.
    do while (flag)
      v = 0.0_dp
      do while (v <= 0.0_dp)
        call rnorm(x)
        v = (1.0_dp + c*x)**3
      end do
      call runif(u)
      if (u < (1.0_dp - 0.0331_dp*(x**4))) then
        fn_val = d*v
        flag = .false.
      else if ( log(u) < (0.5_dp*x*x + d*(1.0_dp - v + log(v))) ) then
        fn_val = d*v
        flag = .false.
      end if
    end do
    if (shape < 1.0_dp) then
      call runif(u)
      fn_val = (fn_val * (u**(1.0_dp/shape))) * scale
    else
      fn_val = fn_val * scale
    end if
  end subroutine rgamma

end subroutine QRc_L_mcmc
