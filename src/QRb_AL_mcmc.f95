! Binary quantile regression with adaptive lasso (ALD link)
subroutine QRb_AL_mcmc(n, k, r, keep, &
     y, tau, X, &
     a, b, delta_init, omega_init, mh_sd0, target_acc, burn_in, clip, zs_sampler, &
     fix_sigma, sigma_fixed, constrain_beta_norm, fix_beta1, norm_slopes_only, &   ! identification-restriction switches (fix_beta1: anchor beta_1 = 1)
     alpha_d, beta_d, kappa_lam, &                   ! Gamma(shape, rate) hyperprior on delta
     alpha_om, beta_om, b0_mean, b0_prec, &          ! omega hyperprior; intercept prior (precision form)
     betadraw, beta0draw, sigmadraw, lambdasqdraw, omegadraw, deltadraw, accdraw, &
     cliphits, &
     beta_init, beta0_init)
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! -------- per-site clamp switches --------
  ! One switch per numerical clamp, replacing the single use_boundaries flag.
  ! The index order is the order in which the sites are reached inside one
  ! sweep of the sampler, and is the contract shared with the R layer: see
  ! .BBQR_CLIP_SITES in internals.R, which must list the same eight names in
  ! the same order.  clip = .true. throughout reproduces the old
  ! use_boundaries = .true. path, and clip = .false. throughout reproduces
  ! use_boundaries = .false.
  ! Sites 1-8 are the clamps the old use_boundaries flag gated together.
  ! Sites 9-15 restore guards that v4 applied UNCONDITIONALLY and v5 deleted.
  ! v4 ran use_boundaries = FALSE, so in v4 sites 1-8 were all off and sites
  ! 9-15 were all on.  That combination, and only that one, is v4's clamp
  ! configuration.
  !
  ! The v4 sites are split by mechanism: flooring a sampler ARGUMENT
  ! (v4_zarg, v4_sarg) and bounding sampled STATE (v4_zbound, v4_sbound) are
  ! different interventions, and lumping them would hide which one carries
  ! the effect.  C_V4RSIG lives inside the sigma draw and is therefore inert
  ! under anchor = "sigma1"; it is kept apart from C_V4ROM so that no single
  ! switch mixes a live guard with a dead one.
  ! The v4 argument floors are split one switch per floor. They are distinct
  ! interventions: C_V4SCHI stops s_j collapsing and so caps how much
  ! shrinkage a coefficient can receive, bounding the effective penalty from
  ! above, while C_V4SLAM and C_V4SPSI act on the shrinkage hyperparameter.
  ! Lumped, a hit count cannot be attributed to any one of them.
  integer, parameter :: NCLIP = 18
  integer, parameter :: C_YSTAR = 1, C_Z = 2, C_S = 3, C_SIGMA = 4, &
                        C_LAMBDASQ = 5, C_OMEGA = 6, C_DELTA = 7, C_MHSD = 8
  integer, parameter :: C_V4SD = 9, C_V4ZPSI = 10, C_V4ZCHI = 11, &
                        C_V4ZBND = 12, &
                        C_V4SLAM = 13, C_V4SPSI = 14, C_V4SCHI = 15, &
                        C_V4SBND = 16, &
                        C_V4RSIG = 17, C_V4ROM = 18

  ! -------- which routine draws the latents z and s --------
  ! Both full conditionals are GIG(1/2).  Two routines implement that draw and
  ! they are NOT numerically equivalent, so which one ran is part of a result's
  ! provenance rather than an implementation detail.
  !
  !   ZS_GIG        rgig_half: q rationalized to avoid cancellation, plus a
  !                 chi == 0 limiting branch.  The current default.
  !   ZS_V4INVGAUS  the routine the v4 build shipped: draw the RECIPROCAL of the
  !                 latent as inverse Gaussian from the textbook q formula and
  !                 invert.  That formula subtracts two nearly equal quantities
  !                 when chi is tiny, which is what v4's guards were written to
  !                 contain -- see the v4 source comments quoted beside
  !                 .BBQR_CLIP_PRESETS in internals.R.
  !
  ! Both consume exactly one normal and one uniform per draw, so two arms that
  ! differ only in this switch stay aligned in the RNG stream and their paired
  ! difference is the routine and nothing else.  The single exception is
  ! rgig_half's chi == 0 branch, which calls rgamma instead; the v4 argument
  ! floors (C_V4ZCHI, C_V4SCHI) make chi strictly positive, so that branch is
  ! unreachable whenever those guards are on.
  integer, parameter :: ZS_GIG = 0, ZS_V4INVGAUS = 1

  ! -------- inputs --------
  integer, intent(in) :: n, k, r, keep
  integer, intent(in) :: y(n)
  real(dp), intent(in) :: tau
  real(dp), intent(in) :: X(n,k)
  real(dp), intent(in) :: a, b       ! Gamma(shape=a, rate=b) prior for inverse-scale sigma
  real(dp), intent(in) :: delta_init, omega_init
  real(dp), intent(in) :: mh_sd0, target_acc
  integer,  intent(in) :: burn_in
  logical,  intent(in) :: clip(NCLIP)
  integer,  intent(in) :: zs_sampler   ! ZS_GIG or ZS_V4INVGAUS; see above
  logical,  intent(in) :: fix_sigma, constrain_beta_norm, fix_beta1
  logical,  intent(in) :: norm_slopes_only   ! TRUE: ||beta||=1 over SLOPES only (Manski);
                                             ! FALSE: legacy joint (beta0,beta) norm   ! fix_beta1: anchor beta_1 = 1
  real(dp), intent(in) :: sigma_fixed
  ! Gamma(shape = alpha_d, rate = beta_d) hyperprior on delta.  Propriety of the
  ! joint posterior requires beta_d > 0: the marginal for delta tends to a
  ! positive constant, so a flat prior leaves it non-integrable.  See the
  ! propriety proposition in the appendix.
  real(dp), intent(in) :: alpha_d, beta_d
  ! Smooth barrier on lambda_j^2: the conditional is multiplied by
  ! exp(-kappa_lam/lambda_j^2), which is conjugate and so only shifts the
  ! inverse-gamma scale. kappa_lam = 0 recovers the unmodified model.
  real(dp), intent(in) :: kappa_lam
  ! -------- derivation-version parameters --------
  ! These four reals select which appendix hierarchy is being sampled.  There
  ! is deliberately no integer `model` switch and no branch on one: each
  ! version is a point in this parameter space, so the code has one code path
  ! and the version is a property of the arguments.  The R layer maps the
  ! names v3/v4/v5 onto these values and records which it used.
  !
  !   Gamma(shape = alpha_om, rate = beta_om) hyperprior on omega.  The omega
  !   full conditional is Gamma(alpha_om + k*delta, beta_om + sum_j
  !   1/lambda_j^2).  V3's p(omega) prop 1/omega is the point
  !   alpha_om = beta_om = 0, which reproduces Gamma(k*delta, sum_j
  !   1/lambda_j^2) exactly rather than approximately.  V4 and V5 require both
  !   strictly positive: with p(omega) prop 1/omega the joint posterior does
  !   not exist, because the wedge on which omega, lambda_j^2, s_j -> 0
  !   together and beta_j = O(sqrt(omega)) leaves the non-integrable measure
  !   d omega / omega.  See Proposition 1a of appendix V4.
  real(dp), intent(in) :: alpha_om, beta_om
  !   N(b0_mean, 1/b0_prec) prior on the intercept, carried as a PRECISION so
  !   that the flat prior of V3 and V4 is the exact interior point
  !   b0_prec = 0 and needs no separate branch.  V5 sets b0_prec = 1/V00 > 0,
  !   which makes the posterior proper unconditionally; under a flat intercept
  !   and a free sigma, propriety instead requires a > 1 in the Gamma(a, b)
  !   prior on sigma plus an overlap condition on the design.
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
  ! Number of times each clamp site BOUND.  Counted whether or not that site's
  ! switch is set, so a run with every clamp off still reports how often each
  ! clamp would have fired -- which is the screen for which sites matter at
  ! all.  Denominators: C_YSTAR, C_Z, C_V4ZPSI, C_V4ZCHI and C_V4ZBND are out
  ! of r*n; C_S, C_LAMBDASQ, C_V4SLAM, C_V4SPSI, C_V4SCHI and C_V4SBND out of
  ! r*k; C_MHSD out of
  ! burn_in; C_DELTA out of the accepted proposals; C_OMEGA and C_V4ROM out
  ! of r; C_SIGMA and C_V4RSIG out of r only when sigma is actually sampled;
  ! and C_V4SD out of n, since v4 floored that sd in the initialisation
  ! sweep alone.  Multiply by the number of guards in a group where it has
  ! more than one: each rate group has two, and the argument floors are now
  ! one switch per floor and so have one each.
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
  real(dp) :: mu_tr, sd_tr, chi, psi, muIG, lamIG, tmpIG
  real(dp) :: chi_i, psi_all
  real(dp) :: psi_g, chi_g          ! sampler arguments after the v4 floors
  real(dp) :: bnorm
  real(dp), allocatable :: beta(:), s(:), lambdasq(:)
  real(dp), allocatable :: z(:), ystar(:), Xbeta(:), resid(:), resid2(:)
  real(dp), parameter :: one=1.0_dp
  real(dp), parameter :: eps_small=1.0e-6_dp, big=1.0e6_dp
  real(dp), parameter :: s_min=1.0e-4_dp, s_max=1.0e4_dp
  real(dp), parameter :: mh_min=1.0e-4_dp, mh_max=2.0_dp

  if (n <= 0 .or. k <= 0) return          ! guarded R-side; never terminate R
  if (r <= 0 .or. keep <= 0) return       ! guarded R-side; never terminate R
  rk = r/keep; if (rk*keep /= r) return   ! guarded R-side; never terminate R

  call rndstart()

  ! Zeroed here rather than beside the other run state, because the
  ! initialisation sweep below already reaches the ystar guard and zeroing
  ! any later would discard its counts.
  cliphits = 0.0_dp

  theta = (one - 2.0_dp*tau) / (tau*(one - tau))
  phisq = 2.0_dp / (tau*(one - tau))

  allocate(beta(k), s(k), lambdasq(k))
  allocate(z(n), ystar(n), Xbeta(n), resid(n), resid2(n))

  beta  = beta_init
  beta0 = beta0_init
  if (fix_beta1) beta(1) = 1.0_dp   ! anchor: hold the first slope at 1, never sampled

  z = 1.0_dp; s = 1.0_dp; lambdasq = 1.0_dp
  sigma = merge(sigma_fixed, 1.0_dp, fix_sigma)
  delta = delta_init; omega = omega_init !

  Xbeta = matmul(X, beta)
  do i = 1, n
    mu_tr = beta0 + Xbeta(i) + theta*z(i)                 ! mean
    ! v4 floored this sd in the INITIALISATION sweep only; its main-loop
    ! counterpart below is deliberately left bare, matching v4, because by
    ! then v4's hard z bound already keeps phisq*z/sigma strictly positive.
    ! So C_V4SD fires at most n times in a whole run, not n per sweep.
    sd_tr = sqrt(floor_rec(phisq * z(i) / sigma, eps_small, &
                           clip(C_V4SD), cliphits(C_V4SD)))   ! sd
    if (y(i) == 1) then
      call rtnorm_geweke(0.0_dp, .true.,  mu_tr, sd_tr, ystar(i))   ! y* > 0
    else
      call rtnorm_geweke(0.0_dp, .false., mu_tr, sd_tr, ystar(i))   ! y* ≤ 0
    end if
  end do

  mh_sd = mh_sd0; accepted = 0; acceptance_rate = 0.0_dp

  outix = 0
  do it = 1, r
    Xbeta = matmul(X, beta)

    ! (1) y* | . -- Eq. (B.9)
    do i=1,n
      mu_tr = beta0 + Xbeta(i) + theta*z(i)
      sd_tr = sqrt(phisq * z(i) / sigma)
      if (y(i) == 1) then
        call rtnorm_geweke(0.0_dp, .true.,  mu_tr, sd_tr, ystar(i))
      else
        call rtnorm_geweke(0.0_dp, .false., mu_tr, sd_tr, ystar(i))
      end if
    end do

    ! (2) beta0 | . -- appendix V5 eq. (beta0 moments); V3/V4 are b0_prec = 0.
    !
    ! Written as prior precision + likelihood precision so that the flat prior
    ! is the exact point b0_prec = 0 and not a large-variance approximation to
    ! it.  With C0 = sum 1/z_i and S0 = sum (ystar_i - x_i'beta - theta z_i)/z_i,
    !
    !   1/V0* = b0_prec + (sigma/phisq)*C0
    !   m0*   = V0* * ( b0_prec*b0_mean + (sigma/phisq)*S0 )
    !
    ! At b0_prec = 0 these reduce ALGEBRAICALLY to V0* = phisq/(sigma*C0) and
    ! m0* = S0/C0, the flat-prior update.  They do not reduce to it in
    ! FLOATING POINT, because the sigma/phisq factor is multiplied in and then
    ! divided out again in a different order, which can move the last bit.  A
    ! last-bit difference is harmless in itself but would make a V3 run
    ! diverge from the earlier build over a few thousand sweeps and so would
    ! contaminate the historical baseline.  The flat case is therefore
    ! evaluated in exactly its historical form.  This is a branch on the value
    ! of one argument, not on a model identifier: the two limbs are the same
    ! distribution, and b0_prec = 0 is the flat prior exactly.
    denom = 0.0_dp; num = 0.0_dp
    do i=1,n
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

    ! (3) z_i | . -- Eq. (B.38), sampled as a reciprocal inverse Gaussian
    psi_all = sigma * (theta*theta + 2.0_dp*phisq) / phisq
    do i=1,n
      resid(i) = ystar(i) - (beta0 + Xbeta(i))
      chi_i    = sigma * (resid(i)*resid(i)) / phisq
      ! v4 floored both GIG arguments and then bounded z itself.  Flooring
      ! chi also suppresses rgig_half's chi == 0 branch, which v4's rinvgaus
      ! did not have -- so v4_zarg is not a pure no-op even where it never
      ! changes a value.
      psi_g = floor_rec(psi_all, eps_small, clip(C_V4ZPSI), cliphits(C_V4ZPSI))
      chi_g = floor_rec(chi_i,   eps_small, clip(C_V4ZCHI), cliphits(C_V4ZCHI))
      if (zs_sampler == ZS_V4INVGAUS) then
        ! v4 drew 1/z as inverse Gaussian and inverted.  The two floors above
        ! are exactly v4's own max(psi, eps) and max(chi, eps), so routing them
        ! through floor_rec leaves this path bit-faithful to v4 while keeping
        ! the guards addressable and their hits counted.
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

    ! (4) beta_j | . -- Eqs. (B.49)--(B.50)
    do j=1,k
      if (fix_beta1 .and. j == 1) cycle   ! anchor: skip beta_1 update, stays at 1
      beta_old_j = beta(j)
      muIG = 0.0_dp
      do i=1,n
        muIG = muIG + (X(i,j)*X(i,j)) / z(i)
      end do
      denom = muIG + (phisq / sigma) / s(j)
      num = 0.0_dp
      do i=1,n
        num = num + X(i,j) * ( ystar(i) - beta0 - theta*z(i) - (Xbeta(i) - X(i,j)*beta_old_j) ) / z(i)
      end do
      mu_j  = num / denom
      var_j = (phisq / sigma) / denom
      call rnorm(zstd); beta(j) = mu_j + sqrt(var_j) * zstd
      do i=1,n
        Xbeta(i) = Xbeta(i) + X(i,j) * (beta(j) - beta_old_j)
      end do
    end do

    ! (4b) optional identification restriction: project (beta0,beta) onto the unit
    ! L2 sphere and adjust the ALD inverse-scale. At the observed-data level,
    ! beta -> c*beta requires sigma -> sigma/c. This mode is an identification
    ! alternative, not part of the Appendix B Gibbs sampler.
    if (constrain_beta_norm) then
      if (norm_slopes_only) then
        bnorm = sqrt(sum(beta*beta))                     ! slope vector only
      else
        bnorm = sqrt(beta0*beta0 + sum(beta*beta))       ! legacy: intercept included
      end if
      if (bnorm > eps_small) then
        beta0 = beta0 / bnorm
        beta  = beta  / bnorm
        Xbeta = Xbeta / bnorm
        sigma = sigma * bnorm
        call clip_rec(sigma, eps_small, big, clip(C_SIGMA), cliphits(C_SIGMA))
      end if
    end if

    ! (5) s_j | . -- Eq. (B.54), sampled as a reciprocal inverse Gaussian
    do j=1,k
      chi   = beta(j)*beta(j)
      psi   = sigma / floor_rec(lambdasq(j), eps_small, clip(C_V4SLAM), &
                                cliphits(C_V4SLAM))
      psi_g = floor_rec(psi, eps_small, clip(C_V4SPSI), cliphits(C_V4SPSI))
      ! chi = beta_j^2. Flooring it stops s_j collapsing, which caps how much
      ! shrinkage a coefficient can take and so bounds the effective penalty.
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

    ! (6) sigma | . -- Eqs. (B.60)--(B.62)
    ! Skipped only in the optional fixed-sigma identification mode.
    if (.not. fix_sigma) then
      do i=1,n
        resid2(i) = ystar(i) - (beta0 + Xbeta(i)) - theta*z(i)
      end do
      scale_g = 0.0_dp
      do i=1,n
        scale_g = scale_g + ( (resid2(i)*resid2(i)) / (2.0_dp*phisq*z(i)) + z(i) )
      end do
      do j=1,k
        scale_g = scale_g + s(j) / (2.0_dp*floor_rec(lambdasq(j), eps_small, &
                                      clip(C_V4RSIG), cliphits(C_V4RSIG)))
      end do
      scale_g = scale_g + b
      muIG = a + k + 1.5_dp*real(n,dp)
      call rgamma(muIG, 1.0_dp/floor_rec(scale_g, eps_small, &
                                        clip(C_V4RSIG), cliphits(C_V4RSIG)), sigma)
      call clip_rec(sigma, eps_small, big, clip(C_SIGMA), cliphits(C_SIGMA))
    end if

    ! (7) lambda_j^2 | . -- Eq. (B.66)
    do j=1,k
      muIG = delta + 1.0_dp
      scale_g = 1.0_dp / ( 0.5_dp*sigma*s(j) + omega + kappa_lam )
      call rgamma(muIG, scale_g, tmpIG)
      lambdasq(j) = 1.0_dp / tmpIG
      call clip_rec(lambdasq(j), s_min, s_max, clip(C_LAMBDASQ), cliphits(C_LAMBDASQ))
    end do

    ! (8) omega | . -- Gamma(alpha_om + k*delta, beta_om + sum_j 1/lambda_j^2).
    !
    ! alpha_om = beta_om = 0 is V3's p(omega) prop 1/omega and reproduces
    ! Gamma(k*delta, sum_j 1/lambda_j^2) exactly.  V4 and V5 pass both
    ! strictly positive; see the declaration comment for why that is required
    ! rather than optional.
    muIG = alpha_om + real(k,dp) * delta
    scale_g = beta_om
    do j=1,k
      scale_g = scale_g + 1.0_dp / floor_rec(lambdasq(j), eps_small, &
                                             clip(C_V4ROM), cliphits(C_V4ROM))
    end do
    call rgamma(muIG, 1.0_dp/floor_rec(scale_g, eps_small, &
                                      clip(C_V4ROM), cliphits(C_V4ROM)), omega)
    call clip_rec(omega, 1.0e-3_dp, 1.0e3_dp, clip(C_OMEGA), cliphits(C_OMEGA))

    ! (9) delta | . -- Eq. (B.74), with delta ~ Gamma(alpha_d, beta_d)
    log_curr = log(delta)
    call rnorm(prop_std)
    log_prop  = log_curr + mh_sd * prop_std
    delta_prop = exp(log_prop)
    lp_curr = logpost_delta(delta,      k, omega, lambdasq, alpha_d, beta_d)
    lp_prop = logpost_delta(delta_prop, k, omega, lambdasq, alpha_d, beta_d)
    ! The proposal is symmetric in eta = log(delta), so add the Jacobian eta.
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

    ! progress ping
    if (mod(it, 100000) == 0) then
        call intpr('Current iteration :', -1, it, 1)
    end if

  end do

  call rndend()

  deallocate(beta, s, lambdasq, z, ystar, Xbeta, resid, resid2)

contains
  ! Clamp `v` into [lo, hi], recording whether the clamp actually binds.
  !
  ! The counter is incremented whether or not `on` is set.  That is deliberate:
  ! an arm run with a clamp OFF still reports how often that clamp WOULD have
  ! fired, which is what makes the activation table readable from the baseline
  ! arm rather than needing the all-clamps-on arm.  Note that once one clamp is
  ! applied the trajectory diverges, so counts for the remaining sites are
  ! comparable within an arm but not across arms.
  !
  ! The arithmetic is the same max(lo, min(v, hi)) the clamps used inline, so
  ! `on = .true.` at every site reproduces the old use_boundaries = .true.
  ! path.  A non-finite `v` compares unequal to everything and is therefore
  ! counted; with `on` false it is left alone and propagates, so the R layer
  ! still records a failed fit rather than a silently repaired one.
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

  subroutine rtnorm_geweke(a, lb, mu, sigma, fn_val)
    implicit none
    logical, intent(in) :: lb
    real(dp), intent(in) :: a, mu, sigma
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
      cliphits(C_YSTAR) = cliphits(C_YSTAR) + 1.0_dp
      if (clip(C_YSTAR)) then
        fn_val = mu
      else
        fn_val = mu + sigma*c
      end if
      return
    end if
    if (abs(c) > 1.0e12_dp) then
      cliphits(C_YSTAR) = cliphits(C_YSTAR) + 1.0_dp
      if (clip(C_YSTAR)) then
        fn_val = mu
        return
      end if
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

end subroutine QRb_AL_mcmc
