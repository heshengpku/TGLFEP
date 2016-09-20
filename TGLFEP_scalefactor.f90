!------------------------------------------------------------
! TGLFEP_scalefactor.f90
!
! PURPOSE:
!  Calculate the EP density threshold 
!  usually called after TGLFEP_ky_widthscan
!------------------------------------------------------------

subroutine TGLFEP_scalefactor

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface 

  implicit none
  integer :: id,np,ierr,STATUS(MPI_STATUS_SIZE)
  integer :: i,n,k,imark
  logical :: iexist
  real :: f0,f1,ft
  integer,parameter :: nfactor = 10
  real,dimension(nfactor) :: factor
  real,dimension(nfactor,nmodes) :: growthrate,growthrate_out &
                                   ,frequency,frequency_out
  real,dimension(nmodes) :: g,f
  real :: gmark,fmark

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  growthrate     = 0.0
  growthrate_out = 0.0
  frequency      = 0.0
  frequency_out  = 0.0

  f0 = 0.0
  f1 = factor_in

  do k = 1,5

    do i = 1,nfactor
      factor(i) = (f1-f0)/nfactor*i+f0
    enddo

    do i = 1+id,nfactor,np
      factor_in = factor(i)

      call TGLFEP_ky

      do n=1,nmodes
        growthrate(i,n) = get_growthrate(n)
        frequency(i,n)  = get_frequency(n)
      enddo
    enddo

    call MPI_BARRIER(TGLFEP_COMM,ierr)

    call MPI_ALLREDUCE(growthrate                      &
                      ,growthrate_out                  &
                      ,nfactor*nmodes                  &
                      ,MPI_DOUBLE_PRECISION            &
                      ,MPI_SUM                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)
    
    call MPI_ALLREDUCE(frequency                       &
                      ,frequency_out                   &
                      ,nfactor*nmodes                  &
                      ,MPI_DOUBLE_PRECISION            &
                      ,MPI_SUM                         &
                      ,TGLFEP_COMM                     &
                      ,ierr)

    if(id .eq. 0) then
      inquire(file=trim('out.scalefactor'//suffix),exist=iexist)
      if(iexist) then
        open(unit=33,file=trim('out.scalefactor'//suffix),status='old',position='append')
      else
        open(unit=33,file=trim('out.scalefactor'//suffix),status='new')
      endif     

      write(33,*) 'mode_flag ',mode_flag_in,'ky ',ky_in,'width ',width_in
      write(33,*)"factor,(gamma(n),freq(n),n=1,nmodes_in)"
      do i = 1,nfactor
        do n = 1,nmodes
          g(n) = growthrate_out(i,n)
          f(n) = frequency_out(i,n)
        enddo
        write(33,10)factor(i),(g(n),f(n),n=1,nmodes)
      enddo

      close(33)
    endif

    imark = 0
    do i = 1,nfactor
      do n = 1,nmodes
        gmark = growthrate_out(i,n)
        fmark = frequency_out(i,n)
        if(gmark .gt. 1.e-7 .and. fmark .lt. freq_AE_upper) then
          f1 = factor(i)
          imark = i
          exit
        endif
      enddo
      if(imark .gt. 0) exit
    enddo

    if(imark .gt. 1) f0 = factor(imark-1)

  enddo !k

  if(imark .eq. 0) then
    factor_in = 10000. !NaN, no threshold found
  else
    factor_in = factor(imark)
  endif
  if(id .eq. 0) then
    print *, 'At',trim(suffix),' Scale Factor is ',factor_in,'with width = ',width_in,'ky = ',ky_in,'and freq = ',fmark
  endif

10 format(F8.4,8F12.7)

end subroutine TGLFEP_scalefactor
