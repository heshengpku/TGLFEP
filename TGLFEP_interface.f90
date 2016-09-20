!---------------------------------------------------------
! TGLFEP_interface.f90
! 7.13.2016
! PURPOSE:
!  The parameters for TGLFEP
!---------------------------------------------------------

module TGLFEP_interface
  
  implicit none

  integer :: TGLFEP_COMM

  integer,parameter :: ns=3, nmodes=4, nky=30
  real,parameter :: ky=0.15

  integer :: mode_flag_in = 2
  real :: width_in = 1.87, ky_in, factor_in = 1.0, factor_max

  integer :: ir = 0, n_toroidal = 3
  real :: q_factor = 1.0, shat_factor = 1.0, scan_factor = 1.0

  real,parameter :: freq_cutoff = -0.2
  real :: freq_AE_upper

  character(16) :: suffix = ''
  
end module TGLFEP_interface
