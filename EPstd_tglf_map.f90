!------------------------------------------------------------
! EPstd_tglf_map.f90
!
! PURPOSE:
!  Mapping from GA-std EP variables to TGLF interface.
!------------------------------------------------------------

subroutine EPstd_tglf_map

  use tglf_interface
  use TGLFEP_interface

  implicit none

  if(mode_flag_in .lt. 1 .or. mode_flag_in .gt. 4) then
    write(*,*) 'mode_flag must be 1. EP+ITG/TEM drive'
    write(*,*) '                  2. EP drive only'
    write(*,*) '                  3. ITG/TEM drive only'
    write(*,*) '                  4. ITG/TEM in EP+ITG/TEM drive'
    write(*,*) 'mode_flag has been set as default 1.'
  endif

  tglf_geometry_flag_in = 0 !s-alpha

  tglf_zs_in(1)   = -1.0  !electron
  tglf_mass_in(1) = 2.723D-4
  tglf_zs_in(2)   = 1.0   !deuterium
  tglf_mass_in(2) = 1.0
  tglf_zs_in(3)   = 2.0   !alpha particle
  tglf_mass_in(3) = 2.0

  tglf_taus_in(1) = 1.0
  tglf_taus_in(2) = 1.0
  tglf_taus_in(3) = 100.0

  factor_max = 0.5*1.0/2./0.025
  if(factor_in .lt. 0.) factor_in = 0.0
  if(factor_in .gt. factor_max) factor_in = factor_max

  tglf_as_in(1)   = 1.0
  tglf_as_in(2)   = 1. - 2.*0.025*factor_in
  tglf_as_in(3)   = 0.025*factor_in

  if(mode_flag_in .eq. 2) then !EP drive only
    tglf_rlns_in(1) = 0.
    tglf_rlns_in(2) = 0.

    tglf_rlts_in(1) = 0.
    tglf_rlts_in(2) = 0.
  else !default
    tglf_rlns_in(1) = 1.
    tglf_rlns_in(2) = 1.

    tglf_rlts_in(1) = 3.
    tglf_rlts_in(2) = 3.
    !-----------------------------------
  endif

  if(mode_flag_in .eq. 3) then !ITG/TEM drive only
    tglf_rlns_in(3) = 0.
    tglf_rlts_in(3) = 0.

    !tglf_filter_in = 2.0   !default
  else
    tglf_rlns_in(3) = 4.0! - scan_factor
    tglf_rlts_in(3) = 0.0! + scan_factor

    tglf_filter_in = 0.0   !The frequency threshold off
  endif

  tglf_q_sa_in    = 2.0*q_factor
  tglf_shat_sa_in = 1.0*shat_factor
  !tglf_alpha_sa_in = 2.0

  if(mode_flag_in .eq. 4) then !ITG/TEM in EP+ITG/TEM drive
    tglf_filter_in = 2.0
  endif

  tglf_use_bper_in = .true.
  tglf_betae_in    = 0.002
  !tglf_adiabatic_elec_in = .true.

  ky_in = 0.01*n_toroidal !0.04!n_toroidal*(2.0*q_factor)/0.5*0.00125 !ky = n*q/(r/a)*rho_star
  freq_AE_upper = freq_cutoff/q_factor*2.6352 !freq_cutoff*omega_TAE
  
end subroutine EPstd_tglf_map
