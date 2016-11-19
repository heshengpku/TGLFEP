!---------------------------------------------------------
! TGLFEP_interface.f90
! 7.13.2016
! PURPOSE:
!  The parameters for TGLFEP
!---------------------------------------------------------

module TGLFEP_interface
  
  implicit none

  integer :: TGLFEP_COMM

  integer,parameter :: nmodes=4

  integer :: process_in
  integer :: threshold_flag

  integer :: mode_in
  logical :: width_in_flag
  real :: width_in, width_min,width_max
  real :: factor_in, factor_max, ky_in, kymark

  integer :: ir, n_toroidal,ky_model
  real :: q_factor = 1.0, shat_factor = 1.0, scan_factor = 1.0

  real,parameter :: freq_cutoff = -0.2
  real :: freq_AE_upper

  integer :: nn
  real,allocatable,dimension(:) :: factor_out

  character(16) :: suffix = ''
  
end module TGLFEP_interface
