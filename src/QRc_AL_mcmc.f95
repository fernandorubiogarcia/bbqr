! Continuous quantile regression with adaptive lasso (ALD working likelihood).
!
! This is the continuous-response counterpart of QRb_AL_mcmc.f95 and implements
! derivations/Mathematical_Derivations_V5C.tex.  Read that appendix before
! changing anything here; the differences from the binary kernel are model
! differences, not refactorings.
!
! Three structural differences from QRb_AL_mcmc:
!
!   (1) y is OBSERVED and real.  The truncated-normal draw of y* -- step (1) of
!       the binary sweep -- and the rtnorm_geweke routine that served it are
!       DELETED, not evaluated at y.  V5C section 1 explains why substituting y
!       into a truncated-normal conditional is not a limiting case of anything.
!
!   (2) sigma is IDENTIFIED and is sampled every sweep.  The binary kernel gates
!       its sigma draw behind `.not. fix_sigma` and the package default
!       (anchor = "sigma1") never executes it.  There is no anchor here: the
!       continuous likelihood is not invariant along the scale ray, so
!       fix_sigma, constrain_beta_norm, fix_beta1 and norm_slopes_only are all
!       absent rather than defaulted off.  Two clip sites that are unreachable
!       in the binary default path (C_SIGMA, C_V4RSIG) are therefore live here.
!
!   (3) The exponent q on sigma in the local-scale prior is an ARGUMENT.  V5's
!       hierarchy is s_j | sigma, lambda_j^2 ~ Exp(sigma / (2 lambda_j^2)),
!       which is q = 1; the response-scale-equivariant form is q = 2.  See
!       Proposition 5C.  q enters in exactly three places, each marked "q ENTERS"
!       below: the psi argument of the s_j draw, the scale of the lambda_j^2
!       draw, and the shape and second exponential term of the sigma draw.
!
! At q = 1 every full conditional reduces algebraically to the binary kernel's
! with y* replaced by y, and the sigma draw is conjugate Gamma.  At q /= 1 the
! sigma conditional leaves the Gamma family and is sampled by a stepping-out
! slice sampler on a log-concave target -- exact, no Metropolis step, no tuning.
!
! NOTE ON THE FILE NAME: R/packages/bbalqr/ carries a vendored copy of
! bayesQR's own QRc_AL_mcmc.f95.  That is a different package, a different
! hierarchy (per-coefficient Gamma priors, no omega/delta layer, a clamped
! sigma rate) and a different DLL.  The two must not be confused when reading
! diffs.
subroutine QRc_AL_mcmc(n, k, r, keep, &
     y, tau, X, &
     a, b, q_exp, sigma_init, &
     delta_init, omega_init, mh_sd0, target_acc, burn_in, clip, zs_sampler, &
     alpha_d, beta_d, &                              ! Gamma(shape, rate) hyperprior on delta
     alpha_om, beta_om, b0_mean, b0_prec, &          ! omega hyperprior; intercept prior (precision form)
     betadraw, beta0draw, sigmadraw, lambdasqdraw, omegadraw, deltadraw, accdraw, &
     cliphits, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- per-site clamp switches --------
  ! Same contract as the binary kernel: one switch per numerical clamp, in the
  ! order the sites are reached inside one sweep, shared with the R layer.
  !
  ! TWO SITES ARE GONE relative to QRb_AL_mcmc's eighteen.  C_YSTAR guarded the
  ! truncated-normal draw and C_V4SD floored its standard deviation in the
  ! initialisation sweep; neither draw exists here.  The remaining sixteen keep
  ! their meanings.  The R layer's .BBQR_CLIP_SITES must carry a SEPARATE
  ! sixteen-name vector for the continuous kernels -- reusing the binary
  ! eighteen-name vector would silently misalign every site.
  integer, parameter :: NCLIP = 16
  integer, parameter :: C_Z = 1, C_S = 2, C_SIGMA = 3, &
                        C_LAMBDASQ = 4, C_OMEGA = 5, C_DELTA = 6, C_MHSD = 7
  integer, parameter :: C_V4ZPSI = 8, C_V4ZCHI = 9, C_V4ZBND = 10, &
                        C_V4SLAM = 11, C_V4SPSI = 12, C_V4SCHI = 13, &
                        C_V4SBND = 14, &
                        C_V4RSIG = 15, C_V4ROM = 16

  ! -------- which routine draws the latents z and s --------
  ! Identical contract to the binary kernel, including the RNG-stream alignment
  ! property: both consume one normal and one uniform per draw.
  integer, parameter :: ZS_GIG = 0, ZS_V4INVGAUS = 1

  ! -------- inputs --------
  integer, intent(in) :: n, k, r, keep
  real(dp), intent(in) :: y(n)          ! OBSERVED continuous response
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  real(dp), intent(in) :: a, b          ! Gamma(shape=a, rate=b) prior for inverse-scale sigma
  ! Exponent q >= 0 in the local-scale prior rate sigma^q / (2 lambda_j^2).
  !   q = 1 : V5's hierarchy; conjugate Gamma draw for sigma.
  !   q = 2 : response-scale equivariant in the penalty block (Proposition 5C);
  !           slice draw for sigma.
  !   q = 0 : Park & Casella's unconditioned hierarchy; the s_j prior drops out
  !           of the sigma conditional entirely.
  ! Values in [0,1) make the sigma log-target non-concave and are not
  ! recommended; the slice sampler still runs but its reliability is not argued.
  real(dp), intent(in) :: q_exp
  real(dp), intent(in) :: sigma_init    ! starting value; sigma is always sampled
  real(dp), intent(in) :: delta_init, omega_init
  real(dp), intent(in) :: mh_sd0, target_acc
  integer,  intent(in) :: burn_in
  logical,  intent(in) :: clip(NCLIP)
  integer,  intent(in) :: zs_sampler    ! ZS_GIG or ZS_V4INVGAUS
  real(dp), intent(in) :: alpha_d, beta_d
  real(dp), intent(in) :: alpha_om, beta_om
  ! N(b0_mean, 1/b0_prec) prior on the intercept, carried as a PRECISION so the
  ! flat prior is the exact interior point b0_prec = 0.
  !
  ! Unlike the binary model, the flat intercept costs NOTHING here: Proposition
  ! 4C gives propriety whenever a + n - 1 > 0, with no overlap condition and no
  ! restriction on a.  The binary requirement a > 1 must NOT be ported.
  real(dp), intent(in) :: b0_mean, b0_prec
  real(dp), intent(in) :: beta_init(k), beta0_init

  ! -------- outputs (thinned by `keep`) --------
  real(dp), intent(out) :: betadraw(r/keep, k)
  real(dp), intent(out) :: beta0draw(r/keep)
  real(dp), intent(out) :: sigmadraw(r/keep)
  real(dp), intent(out) :: lambdasqdraw(r/keep, k)
  real(dp), intent(out) :: omegadraw(r/keep)
  real(dp), intent(out) :: deltadraw(r/keep)
  real(dp), intent(out) :: accdraw(r/keep)
  ! Denominators, as in the binary kernel but with two changes: C_SIGMA and
  ! C_V4RSIG are now out of r unconditionally (sigma is always sampled), and
  ! there is no site with denominator n.
  real(dp), intent(out) :: cliphits(NCLIP)

  ! -------- locals / work arrays / constants --------
  integer :: i, j, it, outix, rk
  real(dp) :: theta, phisq, denom, num, m0, s0sq, prec0
  real(dp) :: sigma, omega, delta, mh_sd, acceptance_rate
  integer  :: accepted
  real(dp) :: zstd, prop_std, log_curr, log_prop, delta_prop
  real(dp) :: lp_curr, lp_prop, log_acc, u
  real(dp) :: gamma_m, scale_g
  real(dp) :: beta0, beta_old_j, var_j, mu_j
  real(dp) :: chi, psi, muIG, lamIG, tmpIG
  real(dp) :: chi_i, psi_all
  real(dp) :: psi_g, chi_g
  real(dp) :: sigq                    ! sigma**q_exp, recomputed at each use site
  real(dp) :: A1, A2, shape_sig
  real(dp), allocatable :: beta(:), s(:), lambdasq(:)
  real(dp), allocatable :: z(:), Xbeta(:), resid(:), resid2(:)
  real(dp), parameter :: one=1.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp
  real(dp), parameter :: s_min=1.0e-4_dp, s_max=1.0e4_dp
  real(dp), parameter :: mh_min=1.0e-4_dp, mh_max=2.0_dp

  if (n <= 0 .or. k <= 0) return          ! guarded R-side; never terminate R
  if (r <= 0 .or. keep <= 0) return       ! guarded R-side; never terminate R
  rk = r/keep; if (rk*keep /= r) return   ! guarded R-side; never terminate R
  if (q_exp < 0.0_dp) return              ! guarded R-side; never terminate R

  call rndstart()

  cliphits = 0.0_dp

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))

  allocate(beta(k), s(k), lambdasq(k))
  allocate(z(n), Xbeta(n), resid(n), resid2(n))

  beta  = beta_init
  beta0 = beta0_init

  z = 1.0_dp; s = 1.0_dp; lambdasq = 1.0_dp
  sigma = sigma_init
  if (sigma <= 0.0_dp) sigma = one
  delta = delta_init; omega = omega_init

  ! No initialisation sweep.  The binary kernel needs one because ystar has no
  ! value until it is drawn; y is data.

  mh_sd = mh_sd0; accepted = 0; acceptance_rate = 0.0_dp

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! ---- step (1) of the binary sweep, y* | . , IS DELETED. ----

    ! (2) beta0 | . -- V5C eq. (fc_b05c).  Identical to the binary update with
    ! y in place of ystar; the intercept enters the augmented likelihood the
    ! same way in both models.  The b0_prec = 0 limb is kept in its historical
    ! floating-point form for the same reason as the binary kernel.
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

    ! (3) z_i | . -- V5C eq. (fc_z5c).  Free of q.
    psi_all = sigma * (theta*theta + 2.0_dp*phisq) / phisq
    do i=1,n
      resid(i) = y(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      psi_g = floor_rec(psi_all, eps_small, clip(C_V4ZPSI), cliphits(C_V4ZPSI))
      chi_g = floor_rec(chi_i,   eps_small, clip(C_V4ZCHI), cliphits(C_V4ZCHI))
      if (zs_sampler == ZS_V4INVGAUS) then
        muIG  = sqrt(psi_g / chi_g)
        lamIG = psi_g
        call rinvgaus(muIG, lamIG, tmpIG)
        z(i) = one / tmpIG
      else
        call rgig_half(psi_g, chi_g, z(i))
      end if
      call clip_rec(z(i), 1.0e-8_dp, 1.0e8_dp, clip(C_V4ZBND), cliphits(C_V4ZBND))
      call clip_rec(z(i), eps_small, big, clip(C_Z), cliphits(C_Z))
    end do

    ! (4) beta_j | . -- V5C eq. (fc_bj5c).  Free of q.
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

    ! ---- step (4b) of the binary sweep, the norm projection, IS DELETED. ----
    ! There is no identification ray: the continuous likelihood is not invariant
    ! under (beta0, beta, sigma) -> (c*beta0, c*beta, sigma/c).

    ! (5) s_j | . -- V5C eq. (fc_sj5c).  *** q ENTERS (1 of 3) ***
    !   psi = sigma^q / lambda_j^2   (binary kernel: sigma / lambda_j^2)
    sigq = sigma**q_exp
    do j=1,k
      chi   = beta(j)*beta(j)
      psi   = sigq / floor_rec(lambdasq(j), eps_small, clip(C_V4SLAM), &
                               cliphits(C_V4SLAM))
      psi_g = floor_rec(psi, eps_small, clip(C_V4SPSI), cliphits(C_V4SPSI))
      chi_g = floor_rec(chi, eps_small, clip(C_V4SCHI), cliphits(C_V4SCHI))
      if (zs_sampler == ZS_V4INVGAUS) then
        muIG  = sqrt(psi_g / chi_g)
        lamIG = psi_g
        call rinvgaus(muIG, lamIG, tmpIG)
        s(j) = one / tmpIG
      else
        call rgig_half(psi_g, chi_g, s(j))
      end if
      call clip_rec(s(j), 1.0e-8_dp, 1.0e8_dp, clip(C_V4SBND), cliphits(C_V4SBND))
      call clip_rec(s(j), s_min, s_max, clip(C_S), cliphits(C_S))
    end do

    ! (6) sigma | . -- V5C eq. (fc_sig5c).  *** q ENTERS (2 of 3) ***
    ! ALWAYS EXECUTED.  The binary kernel wraps this in `if (.not. fix_sigma)`
    ! and the package default never reaches it.
    !
    !   p(sigma | .) prop sigma^(a + q*k + 3n/2 - 1) exp{-A1*sigma - A2*sigma^q}
    !
    ! with A1 the likelihood-plus-z-prior-plus-b term and A2 the penalty term.
    ! At q = 1 the two exponentials merge and the draw is Gamma(a + k + 3n/2,
    ! A1 + A2), which is exactly the binary kernel's (a_1, a_2).
    do i=1,n
      resid2(i) = y(i) - (beta0 + Xbeta(i)) - theta*z(i)
    end do
    A1 = b
    do i=1,n
      A1 = A1 + ( (resid2(i)*resid2(i)) / (2.0_dp*phisq*z(i)) + z(i) )
    end do
    A2 = 0.0_dp
    do j=1,k
      A2 = A2 + s(j) / (2.0_dp*floor_rec(lambdasq(j), eps_small, &
                                         clip(C_V4RSIG), cliphits(C_V4RSIG)))
    end do
    shape_sig = a + q_exp*real(k,dp) + 1.5_dp*real(n,dp)
    if (abs(q_exp - one) < 1.0e-12_dp) then
      ! conjugate limb: identical arithmetic to the binary kernel
      scale_g = A1 + A2
      call rgamma(shape_sig, 1.0_dp/floor_rec(scale_g, eps_small, &
                                              clip(C_V4RSIG), cliphits(C_V4RSIG)), sigma)
    else if (A2 <= 0.0_dp .or. k == 0) then
      ! the penalty term is absent (q = 0 collapses it, k = 0 removes it)
      call rgamma(shape_sig, 1.0_dp/floor_rec(A1, eps_small, &
                                              clip(C_V4RSIG), cliphits(C_V4RSIG)), sigma)
    else
      call slice_sigma(sigma, shape_sig, A1, A2, q_exp)
    end if
    call clip_rec(sigma, eps_small, big, clip(C_SIGMA), cliphits(C_SIGMA))

    ! (7) lambda_j^2 | . -- V5C eq. (fc_lam5c).  *** q ENTERS (3 of 3) ***
    !   scale = sigma^q * s_j / 2 + omega   (binary kernel: sigma * s_j / 2 + omega)
    ! The kappa_lam barrier of the binary kernel is NOT carried: V4 and V5 both
    ! refuse it, and V5C's appendix excludes it.
    sigq = sigma**q_exp
    do j=1,k
      muIG = delta + 1.0_dp
      scale_g = 1.0_dp / ( 0.5_dp*sigq*s(j) + omega )
      call rgamma(muIG, scale_g, tmpIG)
      lambdasq(j) = 1.0_dp / tmpIG
      call clip_rec(lambdasq(j), s_min, s_max, clip(C_LAMBDASQ), cliphits(C_LAMBDASQ))
    end do

    ! (8) omega | . -- V5C eq. (fc_om5c).  Free of sigma and therefore of q;
    ! identical to the binary kernel.
    muIG = alpha_om + real(k,dp) * delta
    scale_g = beta_om
    do j=1,k
      scale_g = scale_g + 1.0_dp / floor_rec(lambdasq(j), eps_small, &
                                             clip(C_V4ROM), cliphits(C_V4ROM))
    end do
    call rgamma(muIG, 1.0_dp/floor_rec(scale_g, eps_small, &
                                      clip(C_V4ROM), cliphits(C_V4ROM)), omega)
    call clip_rec(omega, 1.0e-3_dp, 1.0e3_dp, clip(C_OMEGA), cliphits(C_OMEGA))

    ! (9) delta | . -- random-walk Metropolis on log delta.  Free of q.
    log_curr = log(delta)
    call rnorm(prop_std)
    log_prop  = log_curr + mh_sd * prop_std
    delta_prop = exp(log_prop)
    lp_curr = logpost_delta(delta,      k, omega, lambdasq, alpha_d, beta_d)
    lp_prop = logpost_delta(delta_prop, k, omega, lambdasq, alpha_d, beta_d)
    log_acc = (lp_prop + log(delta_prop)) - (lp_curr + log(delta))
    call runif(u)
    if (log(u) < log_acc) then
      delta = delta_prop
      call clip_rec(delta, 1.0e-3_dp, 1.0e3_dp, clip(C_DELTA), cliphits(C_DELTA))
      accepted = accepted + 1
    end if
    if (it <= burn_in) then
      gamma_m = 1.0_dp / sqrt( real(it,dp) )
      mh_sd = exp( log(mh_sd) + gamma_m * ( merge(1.0_dp, 0.0_dp, log(u) < log_acc ) - target_acc ) )
      call clip_rec(mh_sd, mh_min, mh_max, clip(C_MHSD), cliphits(C_MHSD))
    end if
    acceptance_rate = real(accepted,dp) / real(it,dp)

    if (mod(it, keep) == 0) then
      outix = outix + 1
      betadraw(outix, :)    = beta
      beta0draw(outix)      = beta0
      sigmadraw(outix)      = sigma
      lambdasqdraw(outix,:) = lambdasq
      omegadraw(outix)      = omega
      deltadraw(outix)      = delta
      accdraw(outix)        = acceptance_rate
    end if

    if (mod(it, 100000) == 0) then
        call intpr('Current iteration :', -1, it, 1)
    end if

  end do

  call rndend()

  deallocate(beta, s, lambdasq, z, Xbeta, resid, resid2)

contains

  ! ---------------------------------------------------------------------------
  ! Log of the unnormalised sigma conditional, up to an additive constant:
  !   (sh - 1) log(x) - c1 x - c2 x^qq
  ! ---------------------------------------------------------------------------
  pure real(dp) function lt_sigma(x, sh, c1, c2, qq) result(lv)
    real(dp), intent(in) :: x, sh, c1, c2, qq
    if (x <= 0.0_dp) then
      lv = -huge(1.0_dp)
    else
      lv = (sh - one)*log(x) - c1*x - c2*(x**qq)
    end if
  end function lt_sigma

  ! ---------------------------------------------------------------------------
  ! Stepping-out / shrinkage slice sampler for the sigma conditional at q /= 1.
  !
  ! The target is log-concave whenever sh >= 1 and qq >= 1 (V5C, section on the
  ! q /= 1 case: l''(x) = -(sh-1)/x^2 - qq(qq-1) c2 x^(qq-2) < 0), so a slice
  ! sampler is exact and needs no tuning.  A Metropolis step would need a
  ! proposal scale and would put a second acceptance rate in the output; this
  ! does not.
  !
  ! Both loops are bounded.  If the bounds are hit the current value is
  ! returned unchanged, which is a valid -- if lazy -- Markov transition, and
  ! the caller's C_SIGMA clamp still applies.
  ! ---------------------------------------------------------------------------
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

    if (x <= 0.0_dp) x = one

    ! vertical level: log(y) = l(x) - Exp(1), drawn as l(x) + log(U)
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

    ! stepping out
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

    ! shrinkage
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
    ! bounds exhausted: leave x unchanged
  end subroutine slice_sigma

  ! ---------------------------------------------------------------------------
  ! Everything below is COPIED VERBATIM from QRb_AL_mcmc.f95 (lines 447-576 and
  ! 639-674 at the time of writing), with exactly one deletion: rtnorm_geweke,
  ! whose only caller was the y* draw.  Do not "tidy" these -- keeping them
  ! byte-identical is what makes a binary/continuous paired comparison
  ! attributable to the model rather than to the arithmetic, and rgig_half's
  ! rationalised q and rinvgaus's deliberately un-rationalised one are both
  ! load-bearing (see their own comments).
  ! ---------------------------------------------------------------------------

  subroutine clip_rec(v, lo, hi, on, hits)
    implicit none
    real(dp), intent(inout) :: v
    real(dp), intent(in)    :: lo, hi
    logical,  intent(in)    :: on
    real(dp), intent(inout) :: hits
    real(dp) :: w
    w = max(lo, min(v, hi))
    if (w /= v) hits = hits + 1.0_dp
    if (on) v = w
  end subroutine clip_rec

  ! One-sided version of clip_rec, for the v4 guards that floor a quantity
  ! feeding into a sampler rather than clamping sampled state.  Same counting
  ! convention: the hit is recorded whether or not `on` is set, and the
  ! arithmetic is v4's own max(v, lo).
  real(dp) function floor_rec(v, lo, on, hits) result(w)
    implicit none
    real(dp), intent(in)    :: v, lo
    logical,  intent(in)    :: on
    real(dp), intent(inout) :: hits
    real(dp) :: m
    m = max(v, lo)
    if (m /= v) hits = hits + 1.0_dp
    if (on) then
      w = m
    else
      w = v
    end if
  end function floor_rec

  pure real(dp) function logpost_delta(d, k, omega, lambdasq, alpha_d, beta_d) result(lp)
    implicit none
    real(dp), intent(in) :: d, omega, lambdasq(:), alpha_d, beta_d
    integer,  intent(in) :: k
    integer :: jj
    real(dp) :: sumloglam
    if (.not.(d>0.0_dp .and. omega>0.0_dp)) then
      lp = -1.0d300; return
    end if
    sumloglam = 0.0_dp
    do jj=1,size(lambdasq)
      if (lambdasq(jj) <= 0.0_dp) then
        lp = -1.0d300; return
      end if
      sumloglam = sumloglam + log(lambdasq(jj))
    end do
    ! Eq. (B.74) with p(omega) proportional to omega^(-1) and the Gamma(alpha_d,
    ! beta_d) hyperprior on delta contributing (alpha_d - 1)*log(delta) - beta_d*delta.
    lp = real(k,dp)*d*log(omega) - real(k,dp)*log_gamma(d) - d*sumloglam &
         + (alpha_d - 1.0_dp)*log(d) - beta_d*d
  end function logpost_delta

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

  ! The routine the v4 build used for the same GIG(1/2) draw, reproduced
  ! verbatim from bbqr at c651bca^ so that a v4-sampler arm is the historical
  ! code and not a reconstruction of it.  Draws inverse Gaussian(mu, lambda);
  ! v4's callers inverted the result to get the latent.
  !
  ! q here is the textbook expression.  The trailing sqrt(4*mu*lambda*nu +
  ! mu^2*nu^2) approaches the leading terms as nu*mu/lambda -> 0, so for a tiny
  ! chi -- equivalently a large mu -- q loses precision and can come back as 0
  ! or negative, whereupon 1/q is Inf or the latent goes negative.  rgig_half
  ! rationalizes the same algebra to avoid that.  This is not a bug being
  ! preserved for its own sake: v4's numbers came out of this formula, so
  ! attributing any part of the v4-vs-v6 gap to the sampler requires running it.
  !
  ! `uu` is declared and never used.  It is dead in the v4 original too, and is
  ! kept so that this routine diffs clean against it; -Wall warns about it.
  subroutine rinvgaus(mu, lambda, fn_val)
    implicit none
    real(dp), intent(in)  :: mu, lambda
    real(dp), intent(out) :: fn_val
    real(dp) :: nu, q, z, uu
    call rnorm(nu); nu = nu*nu
    q = mu + (nu*mu*mu)/(2.0_dp*lambda) - (mu/(2.0_dp*lambda))*sqrt(4.0_dp*mu*lambda*nu + mu*mu*nu*nu)
    call runif(z)
    if (z <= (mu/(mu+q))) then
      fn_val = q
    else
      fn_val = (mu*mu)/q
    end if
  end subroutine rinvgaus

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

end subroutine QRc_AL_mcmc
