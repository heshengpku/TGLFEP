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
  integer :: nn = 5
  character(1) :: str_n
  character(3) :: str_c
  integer :: code_in = 6

  TGLFEP_COMM = COMM_IN

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  select case(code_in)
  case(0)
    call TGLFEP_ky !using width_in, ky_in, &mode_flag_in, factor_in

  case(1)
    call TGLFEP_TM !using width_in, &mode_flag_in, factor_in

  case(2)
    call TGLFEP_ky_widthscan !using ky_in, &mode_flag_in, factor_in

  case(3)
    call TGLFEP_ky_nEPscan !using width_in,ky_in

  case(4)
    mode_flag_in = 2 !The width is usually found without background plasma
    call TGLFEP_ky_widthscan !get the width for the maximun growthrate

    key = id / 3
    color = id - key*3
    ! key = id / 4    !case(3) TGLFEP_ky_nEPscan included
    ! color = id - key*4
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    select case(color)
    case(0) !The spectrum with both background plasma and EPs
      mode_flag_in = 1
      call TGLFEP_TM
    case(1) !The spectrum with only EPs
      mode_flag_in = 2
      call TGLFEP_TM
    case(2) !The spectrum with only EPs, but usually only ITG/TEM, excluded TAE/EPM
      mode_flag_in = 4
      call TGLFEP_TM
    case(3) !The gamma and frequency vs nEP, like FIG. 4 in BassWaltz2010
      call TGLFEP_ky_nEPscan 
    end select

  case(5)
    !using width_in, ky_in (default)

    key = id / nn
    color = id - key*nn
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    ! n_toroidal = color + 1
    ! write(str_n,'(I1)')n_toroidal
    ! suffix = trim(suffix)//'_n'//str_n
    
    scan_factor = real(color+1)/10.0   !for a/Ln_EP scan, usually to calculate C_R
    write(str_c,'(F3.1)')scan_factor
    suffix = trim(suffix)//'_c'//str_c

    call TGLFEP_bisection
  case(6)
    !using factor_in, ky_in (default)
    mode_flag_in = 2 !The width is usually found without background plasma
    call TGLFEP_ky_widthscan !get the width for all toroidal n

    ! key = id / 2     !To print the spectrum if needed
    ! color = id - key*2
    ! call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)
    
    ! select case(color)
    ! case(0)
    !  mode_flag_in = 1
    !  call TGLFEP_TM
    ! case(1)
    !  mode_flag_in = 2
    !  call TGLFEP_TM
    ! end select
    
    ! call MPI_BARRIER(COMM_IN,ierr)

    key = id / nn
    color = id - key*nn
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    n_toroidal = color + 1
    write(str_n,'(I1)')n_toroidal
    suffix = trim(suffix)//'_n'//str_n
    
    ! scan_factor = real(color+1)/10.0   !for a/Ln_EP scan, usually to calculate C_R
    ! write(str_c,'(F3.1)')scan_factor
    ! suffix = trim(suffix)//'_c'//str_c

    call TGLFEP_scalefactor
  case default
    print *,'No code_in'
  end select

end subroutine TGLFEP_mainsub
