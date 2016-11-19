!------------------------------------------------------------
! TGLFEP_mainsub.f90
! 7.13.2016
! PURPOSE:
!  Call different subroutines
!------------------------------------------------------------

subroutine TGLFEP_mainsub(COMM_IN)

  use mpi
  use TGLFEP_interface

  implicit none

  integer,intent(in) :: COMM_IN
  integer :: id,np,color,key,ierr,STATUS(MPI_STATUS_SIZE)
  integer :: i
  logical :: spectrum_flag

  character(1) :: str_n
  character(3) :: str_c

  TGLFEP_COMM = COMM_IN

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  select case(process_in)
  case(0)
    call TGLFEP_ky !using width_in, ky_in, &mode_in, factor_in

  case(1)
    call TGLFEP_TM !using width_in, &mode_in, factor_in

  case(2)
    call TGLFEP_ky_nEPscan !using width_in,ky_in

  case(3)
    if(.not. width_in_flag) then
      mode_in = 2 !The width is usually found without background plasma
      call TGLFEP_ky_widthscan !get the width for the maximun growthrate
      kymark = ky_in
    endif

    spectrum_flag = .true.
    if(spectrum_flag) then  !To print the spectrum if needed
      key = id / 3
      color = id - key*3
      ! key = id / 4    !case(3) TGLFEP_ky_nEPscan included
      ! color = id - key*4
      call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

      select case(color)
      case(0) !The spectrum with both background plasma and EPs
        mode_in = 1
        call TGLFEP_TM
      case(1) !The spectrum with only EPs
        mode_in = 2
        call TGLFEP_TM
      case(2) !The spectrum with only EPs, but usually only ITG/TEM, excluded TAE/EPM
        mode_in = 4
        call TGLFEP_TM
      ! case(3) !The gamma and frequency vs nEP, like FIG. 4 in BassWaltz2010
      !   call TGLFEP_ky_nEPscan 
      end select
    endif

  case(4)
    if(.not. width_in_flag) then
      !using factor_in, ky_in (default)
      mode_in = 2 !The width is usually found without background plasma
      call TGLFEP_ky_widthscan !get the width for all toroidal n
      kymark = ky_in
    endif

    spectrum_flag = .false.
    if(spectrum_flag) then  !To print the spectrum if needed
      key = id / 2
      color = id - key*2
      call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)
      
      select case(color)
      case(0)
       mode_in = 1
       call TGLFEP_TM
      case(1)
       mode_in = 2
       call TGLFEP_TM
      end select

      call MPI_BARRIER(COMM_IN,ierr)
    endif

    key = id / nn
    color = id - key*nn
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    if(threshold_flag .eq. 0) then
      n_toroidal = color + 1
      write(str_n,'(I1)')n_toroidal
      suffix = trim(suffix)//'_n'//str_n
    else
      scan_factor = real(color+1)/10.0   !for a/Ln_EP scan, usually to calculate C_R
      write(str_c,'(F3.1)')scan_factor
      suffix = trim(suffix)//'_c'//str_c
    endif

    mode_in = 2 !The density threshold only for gamma_AE > 0 recipe in TGLF
    call TGLFEP_scalefactor

    if(id .eq. 0) then
      factor_out(1) = factor_in
      do i = 1,nn-1
        call MPI_RECV(factor_out(i+1),1,MPI_DOUBLE_PRECISION,i,i,COMM_IN,STATUS,ierr)
      enddo

      factor_in = minval(factor_out)

      ! print *, 'At',trim(suffix),'Scale Factor Min is ',factor_in,&
      !          'with width = ',width_in,'at each n:',(factor_out(i),i=1,nn)

    else if(id .lt. nn) then
      call MPI_SEND(factor_in,1,MPI_DOUBLE_PRECISION,0,id,COMM_IN,STATUS,ierr)
    endif

  case default
    print *,'No process_in'
  end select

end subroutine TGLFEP_mainsub
