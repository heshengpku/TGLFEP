!------------------------------------------------------------
! EPtran_tglf_map.f90
!
! PURPOSE:
!  Mapping from EPtran internal variables to TGLF interface.
!------------------------------------------------------------

subroutine EPtran_tglf_map

  use tglf_interface
  use TGLFEP_interface
  use EPtran_to_tglf

  implicit none
  !integer :: ir = 26

  if(ir .lt. 1 .or. ir .gt. 51) then
    write(*,*) 'ir must be 1<=ir<=nr=51'
    ir = 26
  endif
  ! if(mode_flag_in .lt. 1 .or. mode_flag_in .gt. 4) then
  !   write(*,*) 'mode_flag_in must be 1. EP+ITG/TEM drive'
  !   write(*,*) '                  2. EP drive only'
  !   write(*,*) '                  3. ITG/TEM drive only'
  !   write(*,*) '                  4. ITG/TEM in EP+ITG/TEM drive'
  !   write(*,*) 'mode_flag_in has been set as default 1.'
  ! endif

  tglf_sign_bt_in = -1.0
  tglf_sign_it_in = -1.0

  tglf_zs_in(1) = -1.0
  tglf_zs_in(2) = 1.0
  tglf_zs_in(3) = z_alpha_tglf
  
  tglf_mass_in(1) = 1./1836./2.
  tglf_mass_in(2) = mi_tglf/2.
  tglf_mass_in(3) = m_alpha_tglf/2.

  tglf_taus_in(1) = 1.0
  tglf_taus_in(2) = Ti_hat_tglf_r(ir)
  tglf_taus_in(3) = T_EP_hat_r(ir)

  factor_max = 0.5*1.0/z_alpha_tglf/n_EP_hat_r(ir)
  if(factor_in .lt. 0.) factor_in = 0.0
  if(factor_in .gt. factor_max) factor_in = factor_max

  tglf_as_in(1) = 1.0
  !tglf_as_in(2) = ni_hat_tglf_r(ir)
  tglf_as_in(2) = 1.0 - z_alpha_tglf*n_EP_hat_r(ir)*factor_in
  tglf_as_in(3) = n_EP_hat_r(ir)*factor_in

  if(mode_flag_in .eq. 2) then !EP drive only
    tglf_rlns_in(1) = 0.
    tglf_rlns_in(2) = 0.

    tglf_rlts_in(1) = 0.
    tglf_rlts_in(2) = 0.
  else
    tglf_rlns_in(1) = aoLn_e_tglf_r(ir)
    tglf_rlns_in(2) = aoLn_i_tglf_r(ir)

    tglf_rlts_in(1) = aoLT_e_tglf_r(ir)
    tglf_rlts_in(2) = aoLT_i_tglf_r(ir)
  endif

  !EP species
  if(mode_flag_in .eq. 3) then !ITG/TEM drive only
    tglf_rlns_in(3) = 0.
    tglf_rlts_in(3) = 0.

    !tglf_filter_in = 2.0
  else
    tglf_rlns_in(3) = aoLn_EP_r(ir)*scan_factor
    tglf_rlts_in(3) = aoLT_EP_r(ir)

    tglf_filter_in = 0.0   !The frequency threshold off
  endif

  !----------------------------------------------------------------
  ! Geometry parameters:
  !
  tglf_geometry_flag_in = 1  ! Miller 
  ! s-alpha
  tglf_rmin_sa_in     = r_hat_tglf_r(ir)
  tglf_rmaj_sa_in     = rmaj_hat_tglf_r(ir)
  tglf_q_sa_in        = q_tglf_r(ir)*q_factor
  tglf_shat_sa_in     = s_hat_tglf_r(ir)*shat_factor
  tglf_alpha_sa_in    = 0.0
  tglf_xwell_sa_in    = 0.0
  tglf_theta0_sa_in   = 0.0
  tglf_b_model_sa_in  = 1
  tglf_ft_model_sa_in = 1
  ! Miller
  tglf_rmin_loc_in    = r_hat_tglf_r(ir)
  tglf_rmaj_loc_in    = rmaj_hat_tglf_r(ir)
  tglf_zmaj_loc_in    = 0.0
  tglf_drmajdx_loc_in = drmajdr_tglf_r(ir)
  tglf_dzmajdx_loc_in = 0.0
  tglf_kappa_loc_in   = kappa_tglf_r(ir)
  tglf_s_kappa_loc_in = s_kappa_tglf_r(ir)
  tglf_delta_loc_in   = delta_tglf_r(ir)
  tglf_s_delta_loc_in = s_delta_tglf_r(ir)
  tglf_zeta_loc_in    = 0.0
  tglf_s_zeta_loc_in  = 0.0
  tglf_q_loc_in       = abs(q_tglf_r(ir))*q_factor
  tglf_q_prime_loc_in = q_prime_tglf_r(ir)*q_factor**2.0*shat_factor
  tglf_p_prime_loc_in = p_prime_tglf_r(ir)*q_factor* &
                         (tglf_as_in(1)*tglf_taus_in(1)*(tglf_rlns_in(1)+tglf_rlts_in(1)) + &
                          tglf_as_in(2)*tglf_taus_in(2)*(tglf_rlns_in(2)+tglf_rlts_in(2)) + &
                          tglf_as_in(3)*tglf_taus_in(3)*(tglf_rlns_in(3)+tglf_rlts_in(3)))/ &
                         ((aoLn_e_tglf_r(ir)+aoLT_e_tglf_r(ir)) + &
                          (1.0 - z_alpha_tglf*n_EP_hat_r(ir))*Ti_hat_tglf_r(ir)* &
                          (aoLn_i_tglf_r(ir)+aoLT_i_tglf_r(ir)) + &
                          n_EP_hat_r(ir)*T_EP_hat_r(ir)*(aoLn_EP_r(ir)+aoLT_EP_r(ir)))
  !----------------------------------------------------------------

  tglf_use_bper_in = .true.
  tglf_betae_in = betae_unit_tglf_r(ir)

  !tglf_xnue_in = xnue_tglf_r(ir)
  tglf_xnue_in = 0.
  tglf_zeff_in = zeff_tglf

  if(mode_flag_in .eq. 4) then !ITG/TEM in EP+ITG/TEM drive
    tglf_filter_in = 2.0
  endif

  !ky_in = n_toroidal/3.*q_factor*kymark(ir)
  ky_in = n_toroidal*tglf_q_loc_in/tglf_rmin_loc_in*rho_star(ir) !ky = n*q/(r/a)*rho_star
  !ky_in = n_toroidal*0.1*tglf_zs_in(3)/sqrt(tglf_mass_in(3)*tglf_taus_in(3))
  freq_AE_upper = freq_cutoff/q_factor*omega_TAE(ir) !freq_cutoff*omega_TAE
  
end subroutine EPtran_tglf_map
