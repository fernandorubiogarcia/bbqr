! =====================================================================
! Binary quantile regression with the (standard) LASSO
! =====================================================================
! Faithful implementation of:
!   Benoit, D.F., Alhamzawi, R. & Yu, K. (2013). "Bayesian lasso binary
!   quantile regression." Computational Statistics, 28(6): 2861-2873.
!
! This is the STANDARD lasso (a single, shared penalty for all slopes),
! NOT the adaptive lasso.  It differs from this package's QRb_AL_mcmc in
! exactly the penalty layer:
!
!   QRb_AL_mcmc (adaptive):   one local scale lambda_j^2 per coefficient,
!                             with an omega hyperparameter and an
!                             Inverse-Gamma update for each lambda_j^2.
!
!   QRb_L_mcmc  (this file):  a SINGLE global scale lambda^2 (the paper's
!                             v^2) shared by all coefficients, with the
!                             paper's full hyper-hierarchy
!                                 s_j     ~ GIG(1/2, beta_j^2, lambda^2)
!                                 lambda^2~ Gamma(k + delta, sum(s_j)/2 + tau)
!                                 tau     ~ Gamma(delta, lambda^2)          [eq. after (13)]
!                                 delta   ~ Metropolis, target (13):
!                                            p(delta|tau,lambda^2)
!                                              proportional to (tau*lambda^2)^delta / Gamma(delta)
!
! The hierarchical model and full conditionals are those of Section 3 of
! Benoit et al. (2013).  As in QRb_AL_mcmc we write the ALD likelihood via
! the Kozumi & Kobayashi (2011) normal scale mixture,
!       y*_i | . ~ N( beta0 + x_i'beta + theta z_i , phisq z_i / sigma ),
! with theta = (1-2p)/(p(1-p)) and phisq = 2/(p(1-p)); this is an
! algebraically equivalent augmentation to the (xi, zeta) mixture written
! in the paper and keeps the code identical in style to the rest of the
! package.  In this parameterisation the standard-lasso prior on a slope
! is beta_j ~ N(0, s_j) with a *sigma-free* local variance s_j, so:
!   * the s_j full conditional uses psi = lambda^2 (no sigma), and
!   * the sigma full conditional has shape a + 3n/2 (no +k) and a rate
!     built only from the likelihood -- matching the paper's sigma update.
!
! sigma is sampled by default (as in the paper).  Because a binary
! threshold model only identifies beta up to a positive scale, the same
! identification switches as QRb_AL_mcmc are provided: set fix_sigma to
! hold sigma at sigma_fixed.
!
! Input arguments:
!   - n, k, r, keep  : #obs, #slopes, #iterations, thinning
!   - y              : binary response (0/1), length n
!   - tau            : quantile of interest p in (0,1)
!   - X              : (n x k) regressor matrix (NO intercept column)
!   - a, b           : Gamma(shape=a, rate=b) prior for sigma (paper: 0.1, 0.1)
!   - delta_init     : starting value for delta
!   - tau_init       : starting value for the hyperparameter tau
!   - mh_sd0         : initial RW-MH log-scale sd for delta
!   - target_acc     : target acceptance rate for the delta MH step
!   - burn_in        : iterations over which mh_sd is adapted
!   - use_boundaries : logical, clamp latent quantities for numerical safety
!   - fix_sigma      : logical, hold sigma fixed at sigma_fixed
!   - sigma_fixed    : value sigma is held at when fix_sigma is .true.
!   - alpha_d, beta_d: Gamma(shape=alpha_d, rate=beta_d) prior for delta
!   - beta_init(k), beta0_init : starting values
!
! Output arguments (thinned by `keep`):
!   - betadraw       : (r/keep x k) slope draws
!   - beta0draw      : (r/keep)     intercept draws
!   - sigmadraw      : (r/keep)     sigma draws (= sigma_fixed if fix_sigma)
!   - lambda2draw    : (r/keep)     global penalty lambda^2 (= paper's v^2)
!   - tauhyperdraw   : (r/keep)     hyperparameter tau
!   - deltadraw      : (r/keep)     hyperparameter delta
!   - accdraw        : (r/keep)     running MH acceptance rate for delta
! =====================================================================
subroutine QRb_L_mcmc(n, k, r, keep, &
     y, tau, X, &
     a, b, delta_init, tau_init, mh_sd0, target_acc, burn_in, use_boundaries, &
     fix_sigma, sigma_fixed, constrain_beta_norm, fix_beta1, norm_slopes_only, &
     alpha_d, beta_d, &
     alpha_t, beta_t, b0_mean, b0_prec, &   ! V5: proper tau_h and intercept priors
     betadraw, beta0draw, sigmadraw, lambda2draw, tauhyperdraw, deltadraw, accdraw, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- inputs --------
  integer,  intent(in) :: n, k, r, keep
  integer,  intent(in) :: y(n)
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  ! -------- V5 corrections --------
  ! Gamma(shape = alpha_t, rate = beta_t) hyperprior on the scale
  ! hyperparameter tau_h.  The published specification is p(tau_h) prop
  ! 1/tau_h, which is the point alpha_t = beta_t = 0 and reproduces the
  ! paper's Gamma(delta, lambda^2) conditional exactly.  That prior leaves the
  ! joint posterior improper: with g = lambda^2*tau_h ~ Gamma(delta,1),
  ! u_j = lambda^2*s_j/2 ~ Exp(1) and q_j = beta_j/sqrt(s_j) ~ N(0,1) the prior
  ! measure factorises as (d tau_h / tau_h) times proper factors, and along
  ! tau_h -> 0 one has lambda^2 -> infinity, s_j -> 0 and beta_j = O(sqrt(tau_h)),
  ! so the likelihood is bounded below while the measure is not integrable.
  ! This is the same defect the adaptive-lasso kernel carried in omega.
  real(dp), intent(in) :: alpha_t, beta_t
  ! N(b0_mean, 1/b0_prec) prior on the intercept, carried as a PRECISION so
  ! that the flat prior is the exact interior point b0_prec = 0.
  real(dp), intent(in) :: b0_mean, b0_prec
  real(dp), intent(in) :: a, b                 ! Gamma(shape=a, rate=b) for sigma
  real(dp), intent(in) :: delta_init, tau_init
  real(dp), intent(in) :: mh_sd0, target_acc
  integer,  intent(in) :: burn_in
  logical,  intent(in) :: use_boundaries
  logical,  intent(in) :: fix_sigma
  real(dp), intent(in) :: sigma_fixed
  logical,  intent(in) :: constrain_beta_norm, fix_beta1
  logical,  intent(in) :: norm_slopes_only   ! TRUE: ||beta||=1 over SLOPES only;
                                             ! FALSE: joint (beta0,beta) norm
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

  ! -------- locals / work arrays / constants --------
  integer  :: i, j, it, outix, rk
  integer  :: accepted
  real(dp) :: theta, phisq, denom, num, m0, s0sq, prec0
  real(dp) :: sigma, lambda2, tau_h, delta, mh_sd, acceptance_rate
  real(dp) :: zstd, prop_std, log_curr, log_prop, delta_prop
  real(dp) :: lp_curr, lp_prop, log_acc, u
  real(dp) :: gamma_m, scale_g, sum_s
  real(dp) :: beta0, beta_old_j, var_j, mu_j
  real(dp) :: mu_tr, sd_tr, chi, psi, muIG, lamIG, tmpIG
  real(dp) :: chi_i, psi_all
  real(dp) :: bnorm
  real(dp), allocatable :: beta(:), s(:)
  real(dp), allocatable :: z(:), ystar(:), Xbeta(:), resid(:), resid2(:)
  real(dp), parameter :: one=1.0_dp, zero=0.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp
  real(dp), parameter :: s_min=1.0e-4_dp, s_max=1.0e4_dp
  real(dp), parameter :: mh_min=1.0e-4_dp, mh_max=2.0_dp

  if (n <= 0 .or. k <= 0) return          ! guarded R-side; never terminate R
  if (r <= 0 .or. keep <= 0) return       ! guarded R-side; never terminate R
  rk = r/keep; if (rk*keep /= r) return   ! guarded R-side; never terminate R

  call rndstart()

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))

  allocate(beta(k), s(k))
  allocate(z(n), ystar(n), Xbeta(n), resid(n), resid2(n))

  beta  = beta_init
  beta0 = beta0_init
  if (fix_beta1) beta(1) = 1.0_dp   ! anchor: hold first slope at 1, never sampled

  z = 1.0_dp; s = 1.0_dp
  lambda2 = 1.0_dp
  sigma   = merge(sigma_fixed, 1.0_dp, fix_sigma)
  delta   = delta_init
  tau_h   = tau_init

  ! initialise y*
  Xbeta = matmul(X, beta)
  do i = 1, n
    mu_tr = beta0 + Xbeta(i) + theta*z(i)
    sd_tr = sqrt(phisq * z(i) / sigma)
    if (y(i) == 1) then
      call rtnorm_geweke(0.0_dp, .true.,  mu_tr, sd_tr, ystar(i))
    else
      call rtnorm_geweke(0.0_dp, .false., mu_tr, sd_tr, ystar(i))
    end if
  end do

  mh_sd = mh_sd0; accepted = 0; acceptance_rate = 0.0_dp

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! (1) y* | .
    do i = 1, n
      mu_tr = beta0 + Xbeta(i) + theta*z(i)
      sd_tr = sqrt(phisq * z(i) / sigma)
      if (y(i) == 1) then
        call rtnorm_geweke(0.0_dp, .true.,  mu_tr, sd_tr, ystar(i))
      else
        call rtnorm_geweke(0.0_dp, .false., mu_tr, sd_tr, ystar(i))
      end if
    end do

    ! (2) beta0 | .   ~ N(m0*, V0*) with prior precision b0_prec.
    !   1/V0* = b0_prec + (sigma/phisq)*C0,  C0 = sum 1/z_i
    !   m0*   = V0* * ( b0_prec*b0_mean + (sigma/phisq)*S0 )
    ! b0_prec = 0 is the flat prior, and its limb below is the historical
    ! expression evaluated in its original order: the two forms agree
    ! algebraically but not bitwise, and a last-bit difference would make a
    ! flat-intercept run diverge from the earlier build over a few thousand
    ! sweeps.
    denom = 0.0_dp; num = 0.0_dp
    do i = 1, n
      denom = denom + one / z(i)
      num   = num   + (ystar(i) - Xbeta(i) - theta*z(i)) / z(i)
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

    ! (3) z_i | .   ~ GIG(1/2, chi_i, psi_all)
    psi_all = sigma * (theta*theta + 2.0_dp*phisq) / phisq
    do i = 1, n
      resid(i) = ystar(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      call rgig_half(psi_all, chi_i, z(i))
      if (use_boundaries) z(i) = max(eps_small, min(z(i), big))
    end do

    ! (4) beta_j | .   ~ Normal, prior beta_j ~ N(0, s_j)  (prior precision 1/s_j)
    do j = 1, k
      if (fix_beta1 .and. j == 1) cycle   ! anchor: skip beta_1 update
      beta_old_j = beta(j)
      muIG = 0.0_dp
      do i = 1, n
        muIG = muIG + (X(i,j)*X(i,j)) / z(i)
      end do
      denom = muIG + (phisq / sigma) / s(j)
      num = 0.0_dp
      do i = 1, n
        num = num + X(i,j) * ( ystar(i) - beta0 - theta*z(i) - (Xbeta(i) - X(i,j)*beta_old_j) ) / z(i)
      end do
      mu_j  = num / denom
      var_j = (phisq / sigma) / denom
      call rnorm(zstd); beta(j) = mu_j + sqrt(var_j) * zstd
      do i = 1, n
        Xbeta(i) = Xbeta(i) + X(i,j) * (beta(j) - beta_old_j)
      end do
    end do


    ! (4b) optional identification restriction: rescale (beta0, beta) and adjust
    ! the ALD inverse-scale. The observed-data likelihood is invariant under
    ! (beta0, beta, sigma) -> (c*beta0, c*beta, sigma/c), so this step leaves
    ! the LIKELIHOOD untouched. It does NOT leave the penalised posterior
    ! untouched: the shrinkage layer is stated on a fixed scale, so renormalising
    ! the coefficients moves them relative to the penalty. That asymmetry is the
    ! object of study, not a defect -- see the Identification section of ?bbqr.
    if (constrain_beta_norm) then
      if (norm_slopes_only) then
        bnorm = sqrt(sum(beta*beta))                     ! slope vector only
      else
        bnorm = sqrt(beta0*beta0 + sum(beta*beta))       ! intercept included
      end if
      if (bnorm > eps_small) then
        beta0 = beta0 / bnorm
        beta  = beta  / bnorm
        Xbeta = Xbeta / bnorm
        sigma = sigma * bnorm
        if (use_boundaries) sigma = max(eps_small, min(sigma, big))
      end if
    end if

    ! (5) s_j | .   ~ GIG(1/2, chi=beta_j^2, psi=lambda^2)   [single global lambda^2]
    do j = 1, k
      chi   = beta(j)*beta(j)
      psi   = lambda2
      call rgig_half(psi, chi, s(j))
      if (use_boundaries) s(j) = max(s_min, min(s(j), s_max))
    end do

    ! (6) sigma | .   ~ Gamma(shape = a + 3n/2, rate = likelihood-only + b)
    !     (s_j is sigma-free in the standard lasso, so no +k in the shape
    !      and no sum(s)/lambda term in the rate)
    if (.not. fix_sigma) then
      do i = 1, n
        resid2(i) = ystar(i) - (beta0 + Xbeta(i)) - theta*z(i)
      end do
      scale_g = 0.0_dp
      do i = 1, n
        scale_g = scale_g + ( (resid2(i)*resid2(i)) / (2.0_dp*phisq*z(i)) + z(i) )
      end do
      scale_g = scale_g + b
      muIG = a + 1.5_dp*real(n,dp)
      call rgamma(muIG, 1.0_dp/scale_g, sigma)
      if (use_boundaries) sigma = max(eps_small, min(sigma, big))
    end if

    ! (7) lambda^2 | .   ~ Gamma(shape = k + delta, rate = sum(s_j)/2 + tau_h)
    sum_s = 0.0_dp
    do j = 1, k
      sum_s = sum_s + s(j)
    end do
    muIG = real(k,dp) + delta
    scale_g = 0.5_dp*sum_s + tau_h
    call rgamma(muIG, 1.0_dp/scale_g, lambda2)
    if (use_boundaries) lambda2 = max(s_min, min(lambda2, s_max))

    ! (8) tau (hyper) | . ~ Gamma(alpha_t + delta, beta_t + lambda^2).
    ! alpha_t = beta_t = 0 recovers the published Gamma(delta, lambda^2)
    ! exactly; positive values are the V5 correction.
    call rgamma(alpha_t + delta, 1.0_dp/(beta_t + lambda2), tau_h)
    if (use_boundaries) tau_h = max(1.0e-3_dp, min(tau_h, 1.0e3_dp))

    ! (9) delta | .   RW-MH on log-scale, target eq. (13) times Gamma prior
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

    ! -------- store --------
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

    ! progress ping
    if (mod(it, 100000) == 0) then
      call intpr('Current iteration :', -1, it, 1)
    end if
  end do

  call rndend()

  deallocate(beta, s, z, ystar, Xbeta, resid, resid2)

contains

  ! log target for delta (up to an additive constant): from Benoit et al.
  ! (2013), eq. (13), p(delta | tau, v^2) proportional to (tau*v^2)^delta / Gamma(delta),
  ! multiplied by a Gamma(shape=alpha_d, rate=beta_d) prior on delta.
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

  subroutine rtnorm_geweke(a, lb, mu, sigma, fn_val)
    implicit none
    logical,  intent(in)  :: lb
    real(dp), intent(in)  :: a, mu, sigma
    real(dp), intent(out) :: fn_val
    real(dp) :: z, az, c, u1, u2, phiz
    az = (a - mu)/sigma
    if (lb) then
      c = az
    else
      c = -az
    end if
    ! A NaN truncation point means an upstream latent has already diverged, and
    ! both rejection loops below would spin forever: `u2 < NaN` is always false.
    ! Propagate the non-finite value instead, so the R layer records a failed
    ! fit rather than burning the wall clock.  A large finite c is not a hazard
    ! -- the exponential-proposal branch accepts almost immediately.  The
    ! approximate mu fallback stays confined to boundary-protection mode.
    if (c /= c) then
      if (use_boundaries) then
        fn_val = mu
      else
        fn_val = mu + sigma*c
      end if
      return
    end if
    if (use_boundaries .and. abs(c) > 1.0e12_dp) then
      fn_val = mu
      return
    end if
    if (c < 0.45_dp) then
      do
        call rnorm(z)
        if (z > c) exit
      end do
      if (lb) then
        fn_val = mu + sigma*z
      else
        fn_val = mu - sigma*z
      end if
    else
      do
        call runif(u1)
        z = -log(u1)/c
        phiz = exp(-0.5_dp*z*z)
        call runif(u2)
        if (u2 < phiz) exit
      end do
      z = z + c
      if (lb) then
        fn_val = mu + sigma*z
      else
        fn_val = mu - sigma*z
      end if
    end if
  end subroutine rtnorm_geweke

  ! Marsaglia & Tsang (2000) gamma generator -- same routine as in
  ! QRb_AL_mcmc.  Returns a Gamma(shape, scale) draw (scale = 1/rate).
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

end subroutine QRb_L_mcmc
