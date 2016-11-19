!------------------------------------------------------------
! TGLFEP_ky.f90
!
! PURPOSE:
!  Calculate the growth rate and frequency at a single ky
!------------------------------------------------------------

subroutine TGLFEP_ky

  use tglf_interface
  use TGLFEP_interface

  implicit none
  logical :: iexist

  call TGLFEP_tglf_map

  tglf_use_transport_model_in = .false.

  if(process_in .eq. 0) then
    inquire(file='input.ky',exist=iexist)
    if(iexist) then
      open(unit=33,file='input.ky',status='old')
      read(33,*) ky_in
    endif
    
    tglf_write_wavefunction_flag_in = 1
  endif

  tglf_kygrid_model_in = 0
  tglf_ky_in           = ky_in

  tglf_nmodes_in = nmodes

  !tglf_nbasis_min_in = 16
  !tglf_nbasis_max_in = 16
  tglf_nbasis_min_in = 32
  tglf_nbasis_max_in = 32
  tglf_nxgrid_in     = 32
  
  tglf_width_in      = width_in
  tglf_find_width_in = .false.

  call tglf_run

end subroutine TGLFEP_ky
