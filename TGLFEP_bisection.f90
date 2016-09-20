!------------------------------------------------------------
! TGLFEP_bisection.f90
!
! PURPOSE:
!  Calculate the density threshold using the bisection method
!  The width_in should be given
!------------------------------------------------------------

subroutine TGLFEP_bisection

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface 

  implicit none
  integer :: i,n
  logical :: imark,iexist
  real :: f0,f1
  real,dimension(nmodes) :: g,f

  mode_flag_in = 2 !gamma_AE only

  f0 = 0.0
  f1 = factor_in

  iexist = .false.
  do while(f1-f0 .gt. 1.e-2)
  !do i = 1,10    !when |f1-f0|<err uses too much time to find

    call TGLFEP_ky

    do n=1,nmodes
      g(n) = get_growthrate(n)
      f(n) = get_frequency(n)
    enddo

    !write(*,10)factor_in,(g(n),f(n),n=1,nmodes)

    imark = .false.
    do n = 1,nmodes
      if(g(n) .gt. 1.e-7 .and. f(n) .lt. freq_AE_upper) then
        imark = .true.
        iexist = .true.
        exit
      endif
    enddo
    
    if(imark) then
      f1 = factor_in
      factor_in = (f0 + f1)/2.
    else
      if(iexist) then
        f0 = factor_in
        factor_in = (f0 + f1)/2.
      else
        factor_in = 2.*f1 - f0
        f0 = f1
        f1 = factor_in
      endif
    endif

    if(f0 .gt. factor_max) then
      f1 = 10000. !NaN, no threshold found
      exit
    endif

  !enddo !i
  end do

  print *, 'At',trim(suffix),' Scale Factor is ',f1,'with width = ',width_in,'ky = ',ky_in

10 format(F8.4,8F12.7)

end subroutine TGLFEP_bisection
