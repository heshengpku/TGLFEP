subroutine TGLFEP_nEPscan

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface

  implicit none
  integer :: id,np,ierr,color,key,tglf_comm,tglf_id,STATUS(MPI_STATUS_SIZE)

  integer :: nfactor = 10

  character(3) :: str_f
  integer :: i,n
  real :: ky_out,growthrate_out(nmodes),frequency_out(nmodes)

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  key = id / (nfactor)
  color = id - key*(nfactor)

  call MPI_COMM_SPLIT(TGLFEP_COMM,color,key,tglf_comm,ierr)

  call MPI_COMM_RANK(tglf_comm,tglf_id,ierr)

  call tglf_init('',tglf_comm)

  factor_in = real(color+1.)/real(nfactor)

  call TGLFEP_tglf_map

  tglf_use_transport_model_in = .true.

  tglf_ns_in = ns

  tglf_kygrid_model_in = 0
  tglf_ky_in           = ky
  tglf_nky_in          = nky

  tglf_nmodes_in = nmodes

  tglf_nbasis_min_in = 32
  tglf_nbasis_max_in = 32
  tglf_nxgrid_in     = 32
  
  tglf_width_in      = width_in
  tglf_find_width_in = .false.
  
  call tglf_run_mpi

  if(tglf_id .eq. 0) then
    write(str_f,'(F3.1)')factor_in
    open(unit=33,file="out.eigenvalue_f"//str_f,status='replace')

    write(33,*)"gyro-bohm normalized eigenvalue spectra for mode_flag ",&
                mode_flag_in,"factor ",factor_in,"width ",width_in
    write(33,*)"ky,(gamma(n),freq(n),n=1,nmodes_in)"
    do i=1,nky
      ky_out = get_ky_spectrum_out(i)
      do n=1,nmodes
        growthrate_out(n) = get_eigenvalue_spectrum_out(1,i,n)
        frequency_out(n) = get_eigenvalue_spectrum_out(2,i,n)
      enddo
      write(33,10)ky_out,(growthrate_out(n),frequency_out(n),n=1,nmodes)
    enddo

    close(33)
  endif

10 format(F8.4,8F12.7)

end subroutine TGLFEP_nEPscan