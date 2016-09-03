!---------------------------------------------------------
! TGLFEP_interface.f90
! 7.13.2016
! PURPOSE:
!  The parameters to TGLFEP
!
!---------------------------------------------------------
module TGLFEP_interface
  
  implicit none

  integer :: TGLFEP_COMM

  integer,parameter :: ns=3,nmodes=4,nky=30
  real,parameter :: ky=0.15 !=1./5/sqrt(2.) !=0.15

  integer :: mode_flag_in
  real :: width_in,ky_in,factor_in,factor_max

  real :: gmark,fmark

  integer :: ir = 0,n_toroidal = 3
  real :: q_factor = 1.0, shat_factor = 1.0, scan_factor = 1.0

  real,parameter :: freq_cutoff = -0.2
  real :: freq_AE_upper

  character(9) :: suffix = ''
  
end module TGLFEP_interface
