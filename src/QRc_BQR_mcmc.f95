! =====================================================================
! Continuous quantile regression, UNPENALISED
! =====================================================================
! Continuous-response counterpart of QRb_BQR_mcmc.f95; implements the
! unpenalised section of derivations/Mathematical_Derivations_V5C.tex.
!
! *** THIS IS THE ONE LAYER THAT GAINS CODE RATHER THAN LOSING IT. ***
! QRb_BQR_mcmc holds sigma at 1 and never samples it -- it has no sigma draw
! at all, and carries no (a, b) arguments, because in a binary threshold model
! sigma is exactly unidentified and sigma = 1 IS the identification anchor.
! Once the response is observed, sigma is identified by the data and must be
! sampled, so this kernel adds:
!
!   * the Gamma(shape = a, rate = b) prior arguments a and b;
!   * a starting value sigma_init;
!   * step (6), sigma | . ~ Gamma(a + 3n/2, A1), V5C eq. (fc_sig_none5c);
!   * sigmadraw now varies rather than being a column of ones.
!
! There is no penalty layer, so the exponent q of the other two continuous
! kernels does not arise: the sigma conditional has shape a + 3n/2 for every q.
! For symmetry of the R interface the argument is accepted and ignored, and a
! comment marks the single place it would enter if a penalty were added.
!
! Hierarchy:
!   beta_j            ~ N(0, beta_var)          [fixed, never updated]
!   beta0             ~ N(b0_mean, 1/b0_prec)
!   sigma             ~ Gamma(a, b)
!
! Propriety: Proposition 1C with k = 0.  The bound
! Z <= (tau(1-tau))^n * Gamma(a+n) / (Gamma(a) b^n) holds for every design,
! every n and every beta_var < infinity.  With a FLAT intercept (b0_prec = 0)
! Proposition 4C gives propriety whenever a + n - 1 > 0 -- no overlap condition
! and no restriction on a.  The binary kernel's admissibility argument, which
! needs the two outcome classes to overlap, does NOT transfer and must not be
! quoted here.
!
! Differences from QRb_BQR_mcmc, all of them model differences:
!   - y is observed and real; the y* truncated-normal draw and rtnorm_geweke
!     are DELETED.
!   - sigma is sampled (see above); constrain_beta_norm, fix_beta1 and
!     norm_slopes_only are absent, there being no scale ray to anchor.
subroutine QRc_BQR_mcmc(n, k, r, keep, &
     y, tau, X, &
     beta_var, a, b, q_exp, sigma_init, use_boundaries, &
     b0_mean, b0_prec, &
     betadraw, beta0draw, sigmadraw, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- inputs --------
  integer,  intent(in) :: n, k, r, keep
  real(dp), intent(in) :: y(n)              ! OBSERVED continuous response
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  real(dp), intent(in) :: beta_var          ! fixed slope prior variance
  real(dp), intent(in) :: a, b              ! Gamma(shape=a, rate=b) for sigma
  ! Accepted for interface symmetry with QRc_AL_mcmc and QRc_L_mcmc and
  ! deliberately unused: with no penalty layer there is no s_j prior for sigma^q
  ! to appear in.  Declared so the R layer can pass one prior object to all
  ! three kernels without special-casing.
  real(dp), intent(in) :: q_exp
  real(dp), intent(in) :: sigma_init
  ! N(b0_mean, 1/b0_prec) prior on the intercept, as a PRECISION so the flat
  ! prior is the exact point b0_prec = 0.  Admissible here by Proposition 4C
  ! whenever a + n - 1 > 0; the binary kernel's overlap condition does not
  ! apply.
  real(dp), intent(in) :: b0_mean, b0_prec
  logical,  intent(in) :: use_boundaries
  real(dp), intent(in) :: beta_init(k), beta0_init

  ! -------- outputs (thinned by `keep`) --------
  real(dp), intent(out) :: betadraw(r/keep, k)
  real(dp), intent(out) :: beta0draw(r/keep)
  real(dp), intent(out) :: sigmadraw(r/keep)

  ! -------- locals --------
  integer :: i, j, it, outix, rk
  real(dp) :: theta, phisq, denom, num, m0, s0sq, prec0
  real(dp) :: sigma, prior_prec
  real(dp) :: zstd
  real(dp) :: beta0, beta_old_j, var_j, mu_j
  real(dp) :: muIG, psi_all, chi_i, A1
  real(dp), allocatable :: beta(:)
  real(dp), allocatable :: z(:), Xbeta(:), resid(:), resid2(:)
  real(dp), parameter :: one=1.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp

  if (n <= 0 .or. k <= 0) return
  if (r <= 0 .or. keep <= 0) return
  rk = r/keep; if (rk*keep /= r) return
  if (beta_var <= 0.0_dp) return

  call rndstart()

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))
  prior_prec = one / beta_var

  allocate(beta(k))
  allocate(z(n), Xbeta(n), resid(n), resid2(n))

  beta  = beta_init
  beta0 = beta0_init
  z = 1.0_dp
  sigma = sigma_init
  if (sigma <= 0.0_dp) sigma = one

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! ---- the y* draw of the binary sweep IS DELETED. ----

    ! (2) beta0 | .
    denom = 0.0_dp; num = 0.0_dp
    do i = 1, n
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
    do i = 1, n
      resid(i) = y(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      call rgig_half(psi_all, chi_i, z(i))
      if (use_boundaries) z(i) = max(eps_small, min(z(i), big))
    end do

    ! (4) beta_j | . ~ Normal, with fixed vague prior N(0, beta_var)
    do j = 1, k
      beta_old_j = beta(j)
      muIG = 0.0_dp
      do i = 1, n
        muIG = muIG + (X(i,j)*X(i,j)) / z(i)
      end do
      denom = muIG + (phisq / sigma) * prior_prec
      num = 0.0_dp
      do i = 1, n
        num = num + X(i,j) * ( y(i) - beta0 - theta*z(i) - (Xbeta(i) - X(i,j)*beta_old_j) ) / z(i)
      end do
      mu_j  = num / denom
      var_j = (phisq / sigma) / denom
      call rnorm(zstd); beta(j) = mu_j + sqrt(var_j) * zstd
      do i = 1, n
        Xbeta(i) = Xbeta(i) + X(i,j) * (beta(j) - beta_old_j)
      end do
    end do

    ! ---- the norm-projection block of the binary sweep IS DELETED. ----

    ! (6) sigma | . ~ Gamma(a + 3n/2, A1).   *** NEW RELATIVE TO THE BINARY
    ! KERNEL, WHICH HAS NO SIGMA DRAW AT ALL. ***  V5C eq. (fc_sig_none5c).
    ! With a penalty layer this shape would be a + q*k + 3n/2 and the rate would
    ! gain the A2 term; with no penalty layer both are absent for every q.
    do i = 1, n
      resid2(i) = y(i) - (beta0 + Xbeta(i)) - theta*z(i)
    end do
    A1 = b
    do i = 1, n
      A1 = A1 + ( (resid2(i)*resid2(i)) / (2.0_dp*phisq*z(i)) + z(i) )
    end do
    call rgamma(a + 1.5_dp*real(n,dp), 1.0_dp/A1, sigma)
    if (use_boundaries) sigma = max(eps_small, min(sigma, big))

    if (mod(it, keep) == 0) then
      outix = outix + 1
      betadraw(outix, :) = beta
      beta0draw(outix)   = beta0
      sigmadraw(outix)   = sigma
    end if

    if (mod(it, 100000) == 0) then
      call intpr('Current iteration :', -1, it, 1)
    end if

  end do

  call rndend()

  deallocate(beta, z, Xbeta, resid, resid2)

contains

  ! ---------------------------------------------------------------------------
  ! COPIED VERBATIM from QRb_BQR_mcmc.f95 (lines 235-283 and 285-320), minus
  ! rtnorm_geweke, whose only caller was the y* draw.  Keep byte-identical.
  ! ---------------------------------------------------------------------------

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

end subroutine QRc_BQR_mcmc
