! =====================================================================
! Binary quantile regression -- Bayesian, NO penalty
! =====================================================================
! Faithful implementation of:
!   Benoit, D.F. & Van den Poel, D. (2012). "Binary Quantile Regression:
!   A Bayesian Approach Based on the Asymmetric Laplace Distribution."
!   Journal of Applied Econometrics, 27(7): 1174-1188.
!
! The 2012 paper augments the binary response with latent y* ~ ALD and
! places a *vague* (non-informative) Normal prior on the regression
! parameters -- there is no shrinkage/variable-selection layer.  The
! scale of the ALD is fixed to sigma = 1 for identification (a binary
! threshold model identifies beta only up to a positive scale; see the
! paper's Section 3.3).
!
! This routine keeps the *exact same code style* as QRb_AL_mcmc (the
! adaptive-lasso kernel of this package): the ALD likelihood is written
! as the Kozumi & Kobayashi (2011) location-scale mixture of normals,
!       y*_i | . ~ N( beta0 + x_i'beta + theta z_i , phisq z_i / sigma )
!       z_i        latent mixing variable (GIG / inverse-Gaussian),
! with theta = (1-2p)/(p(1-p)) and phisq = 2/(p(1-p)).  Sampling y* from
! this augmentation with a Gibbs step is exactly the "improved (Gibbs
! instead of Metropolis-Hastings)" sampler used in Benoit's own bayesQR
! package for QRb, and targets the same posterior as the 2012 paper.
!
! Difference vs. QRb_AL_mcmc:
!   * sigma is held at 1 (never sampled).
!   * beta_j has a fixed Normal prior N(0, beta_var); i.e. the per-
!     coefficient local variance s_j is a *constant* (= beta_var) and is
!     never updated.  There is no s_j / lambda^2 / omega / delta layer.
!
! Input arguments:
!   - n              : number of observations
!   - k              : number of slope regressors (intercept handled separately)
!   - r              : number of MCMC iterations
!   - keep           : thinning parameter
!   - y              : binary response (0/1), length n
!   - tau            : quantile of interest p in (0,1)
!   - X              : (n x k) matrix of regressors (NO intercept column)
!   - beta_var       : prior variance of each slope (Benoit 2012 uses 100)
!   - use_boundaries : logical, clamp latent quantities for numerical safety
!   - beta_init      : (k) starting values for the slopes
!   - beta0_init     : starting value for the intercept
!
! Output arguments (thinned by `keep`):
!   - betadraw       : (r/keep x k) slope draws
!   - beta0draw      : (r/keep)     intercept draws
! =====================================================================
subroutine QRb_BQR_mcmc(n, k, r, keep, &
     y, tau, X, &
     beta_var, use_boundaries, &
     constrain_beta_norm, fix_beta1, norm_slopes_only, &
     b0_mean, b0_prec, &                    ! V5: proper intercept prior
     betadraw, beta0draw, sigmadraw, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- inputs --------
  integer,  intent(in) :: n, k, r, keep
  integer,  intent(in) :: y(n)
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  real(dp), intent(in) :: beta_var
  ! N(b0_mean, 1/b0_prec) prior on the intercept, as a PRECISION so that the
  ! flat prior of Benoit and Van den Poel (2012) is the exact point
  ! b0_prec = 0.  With sigma held at 1 there is no unidentified scale ray, so
  ! the flat intercept here is admissible whenever the two outcome classes
  ! overlap in the convex hull of the slope covariates; a proper prior removes
  ! that condition and is the only improper component this kernel had.
  real(dp), intent(in) :: b0_mean, b0_prec
  logical,  intent(in) :: use_boundaries
  logical,  intent(in) :: constrain_beta_norm, fix_beta1
  logical,  intent(in) :: norm_slopes_only   ! TRUE: ||beta||=1 over SLOPES only;
                                             ! FALSE: joint (beta0,beta) norm
  real(dp), intent(in) :: beta_init(k), beta0_init

  ! -------- outputs (thinned by `keep`) --------
  real(dp), intent(out) :: betadraw(r/keep, k)
  real(dp), intent(out) :: beta0draw(r/keep)
  real(dp), intent(out) :: sigmadraw(r/keep)

  ! -------- locals / work arrays / constants --------
  integer  :: i, j, it, outix, rk
  real(dp) :: theta, phisq, denom, num, m0, s0sq, prec0
  real(dp) :: sigma, prior_prec
  real(dp) :: zstd, beta0, beta_old_j, var_j, mu_j
  real(dp) :: mu_tr, sd_tr, muIG, lamIG, tmpIG
  real(dp) :: chi_i, psi_all
  real(dp) :: bnorm
  real(dp), allocatable :: beta(:)
  real(dp), allocatable :: z(:), ystar(:), Xbeta(:), resid(:)
  real(dp), parameter :: one=1.0_dp, zero=0.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp

  if (n <= 0 .or. k <= 0) return          ! guarded R-side; never terminate R
  if (r <= 0 .or. keep <= 0) return       ! guarded R-side; never terminate R
  rk = r/keep; if (rk*keep /= r) return   ! guarded R-side; never terminate R

  call rndstart()

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))

  allocate(beta(k))
  allocate(z(n), ystar(n), Xbeta(n), resid(n))

  beta  = beta_init
  beta0 = beta0_init
  if (fix_beta1) beta(1) = 1.0_dp   ! anchor: hold first slope at 1, never sampled

  z     = 1.0_dp
  sigma = 1.0_dp                       ! fixed for identification (Benoit 2012)
  prior_prec = one / max(beta_var, eps_small)   ! 1 / prior variance of each slope

  ! initialise y* consistent with observed 0/1
  Xbeta = matmul(X, beta)
  do i = 1, n
    mu_tr = beta0 + Xbeta(i) + theta*z(i)
    sd_tr = sqrt(phisq * z(i) / sigma)
    if (y(i) == 1) then
      call rtnorm_geweke(0.0_dp, .true.,  mu_tr, sd_tr, ystar(i))   ! y* > 0
    else
      call rtnorm_geweke(0.0_dp, .false., mu_tr, sd_tr, ystar(i))   ! y* <= 0
    end if
  end do

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! (1) y* | .   ~ truncated normal (Geweke, 1991)
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
    ! b0_prec = 0 takes the flat-prior limb, which is the historical
    ! expression in its original evaluation order so that a flat run stays
    ! bitwise identical to the earlier build.
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

    ! (3) z_i | .   ~ GIG(1/2, chi_i, psi_all)  (sampled as 1/InvGaussian)
    psi_all = sigma * (theta*theta + 2.0_dp*phisq) / phisq
    do i = 1, n
      resid(i) = ystar(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      call rgig_half(psi_all, chi_i, z(i))
      if (use_boundaries) z(i) = max(eps_small, min(z(i), big))
    end do

    ! (4) beta_j | .   ~ Normal, with fixed vague prior N(0, beta_var)
    do j = 1, k
      if (fix_beta1 .and. j == 1) cycle   ! anchor: skip beta_1 update
      beta_old_j = beta(j)
      muIG = 0.0_dp
      do i = 1, n
        muIG = muIG + (X(i,j)*X(i,j)) / z(i)
      end do
      denom = muIG + (phisq / sigma) * prior_prec
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

    ! -------- store --------
    if (mod(it, keep) == 0) then
      outix = outix + 1
      betadraw(outix, :) = beta
      beta0draw(outix)   = beta0
      sigmadraw(outix)   = sigma
    end if

    ! progress ping
    if (mod(it, 100000) == 0) then
      call intpr('Current iteration :', -1, it, 1)
    end if
  end do

  call rndend()

  deallocate(beta, z, ystar, Xbeta, resid)

contains

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

  ! One draw from the inverse-Gaussian(mu, lambda) using Michael, Schucany
  ! & Haas (1976) -- same routine as in QRb_AL_mcmc.
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

  ! One draw from a truncated normal via Geweke (1991) -- same routine as
  ! in QRb_AL_mcmc.  lb=.true. => support (a, +Inf); lb=.false. => (-Inf, a).
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

end subroutine QRb_BQR_mcmc
