!------------------------------------------------------------
! TGLFEP_ky_widthscan.f90
!
! PURPOSE:
!  Calculate the growth rate and frequency vs width
!  Find the suitable width for max gamma_AE
!------------------------------------------------------------

subroutine TGLFEP_ky_widthscan

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  logical :: iexist
  character(19) :: str_file
  integer :: i,n
  integer :: nwidth
  real,allocatable,dimension(:) :: width
  real,allocatable,dimension(:,:) :: growthrate,growthrate_out &
                                    ,frequency,frequency_out
  real :: w,g(nmodes),f(nmodes)
  logical :: find_max
  real :: gmark,fmark

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  nwidth = nint((width_max - width_min)/0.01)

  allocate(width(nwidth))
  allocate(growthrate(nwidth,nmodes))
  allocate(growthrate_out(nwidth,nmodes))
  allocate(frequency(nwidth,nmodes))
  allocate(frequency_out(nwidth,nmodes))
    
  growthrate     = 0.0
  growthrate_out = 0.0
  frequency      = 0.0
  frequency_out  = 0.0

  do i = 1,nwidth
    width(i) = width_min + 0.01*(i-1)
  enddo

  do i = 1+id,nwidth,np

    width_in = width(i)

    call TGLFEP_ky

    do n = 1,nmodes
      growthrate(i, n) = get_growthrate(n)
      frequency(i, n)  = get_frequency(n)
    enddo

  enddo

  call MPI_BARRIER(TGLFEP_COMM,ierr)
  call MPI_ALLREDUCE(growthrate                      &
                    ,growthrate_out                  &
                    ,nwidth*nmodes                   &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)
  
  call MPI_ALLREDUCE(frequency                       &
                    ,frequency_out                   &
                    ,nwidth*nmodes                   &
                    ,MPI_DOUBLE_PRECISION            &
                    ,MPI_SUM                         &
                    ,TGLFEP_COMM                     &
                    ,ierr)

  if(id .eq. 0) then
    write(str_file,'(A18,I1)')'out.ky_widthscan_m',mode_in
    open(unit=33,file=trim(str_file//suffix),status='replace')
    ! inquire(file=trim(str_file//suffix),exist=iexist)
    ! if(iexist) then
    !   open(unit=33,file=trim(str_file//suffix),status='old',position='append')
    ! else
    !   open(unit=33,file=trim(str_file//suffix),status='new')
    ! endif

    write(33,*)"widthscan at ky =",ky_in,'mode_flag ',mode_in,'factor ',factor_in
    write(33,*)"width,(gamma(n),freq(n),n=1,nmodes_in)"
    do i=1,nwidth
      w = width(i)
      do n=1,nmodes
        g(n) = growthrate_out(i,n)
        f(n) = frequency_out(i,n)
      enddo
      write(33,10)w,(g(n),f(n),n=1,nmodes)
    enddo

    close(33)
  endif

  find_max = .true. 
  !find_max = .false.
  if(find_max) then
    width_in = width_min !- 0.01
    gmark = 0.0
    fmark = 0.0
    do i = 1,nwidth
      do n = 1,nmodes
        if(frequency_out(i,n) .lt. freq_AE_upper .and. growthrate_out(i,n) .gt. gmark) then
          gmark = growthrate_out(i,n)
          fmark = frequency_out(i,n)
          width_in = width(i)
        endif
      enddo
    enddo
    !if(id .eq. 0) print *,ir,'Find the width = ',width_in,'for mode_flag ',mode_in,', factor ',factor_in,', ky ',ky_in
  endif

  deallocate(width)
  deallocate(growthrate)
  deallocate(growthrate_out)
  deallocate(frequency)
  deallocate(frequency_out)

10 format(F5.2,8F12.7)

end subroutine TGLFEP_ky_widthscan
