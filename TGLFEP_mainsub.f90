!------------------------------------------------------------
! TGLFEP_mainsub.f90
! 7.13.2016
! PURPOSE:
!  Mapping parameters to TGLF interface.
!  Finding the width for TAE
!------------------------------------------------------------

subroutine TGLFEP_mainsub(COMM_IN)

  use mpi
  !use tglf_interface
  !use tglf_pkg
  use TGLFEP_interface

  implicit none
  integer,intent(in) :: COMM_IN
  integer :: id,np,color,key,ierr,STATUS(MPI_STATUS_SIZE)
  integer :: nn = 4
  character(1) :: str_n
  integer :: code_in = 5

  TGLFEP_COMM = COMM_IN

  call MPI_COMM_RANK(TGLFEP_COMM,id,ierr)
  call MPI_COMM_SIZE(TGLFEP_COMM,np,ierr)

  select case(code_in)
  case(0)
    call TGLFEP_ky
  case(1)
    call TGLFEP_TM !using width_in, &mode_flag_in, factor_in
  case(2)
    call TGLFEP_nEPscan !using width_in, &mode_flag_in, factor_in scan
  case(3)
    call TGLFEP_ky_widthscan !using ky_in, &mode_flag_in, factor_in
  case(4)
    call TGLFEP_ky_nEPscan !using width_in,ky_in, &mode_flag_in scan, factor_in scan
  case(5)
    call TGLFEP_ky_widthscan !get the width for the maximun growthrate

    key = id / 4
    color = id - key*4
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    select case(color)
    case(0)
      mode_flag_in = 1
      call TGLFEP_TM
    case(1)
      mode_flag_in = 2
      call TGLFEP_TM
    case(2)
      mode_flag_in = 4
      call TGLFEP_TM
    case(3)
      call TGLFEP_ky_nEPscan 
    end select
  case(6)
    call TGLFEP_bisection
  case(7)
    call TGLFEP_scalefactor
  case(8)
    call TGLFEP_ky_widthscan !get the width for all toroidal n

    key = id / 2
    color = id - key*2
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    select case(color)
    case(0)
      mode_flag_in = 1
      call TGLFEP_TM
    case(1)
      mode_flag_in = 2
      call TGLFEP_TM
    end select

    call MPI_BARRIER(COMM_IN,ierr)

    key = id / nn
    color = id - key*nn
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    n_toroidal = color + 1

    write(str_n,'(I1)')n_toroidal
    suffix = trim(suffix)//'_n'//str_n

    call TGLFEP_bisection
  case(9)
    key = id / nn
    color = id - key*nn
    call MPI_COMM_SPLIT(COMM_IN,color,key,TGLFEP_COMM,ierr)

    n_toroidal = color + 1

    write(str_n,'(I1)')n_toroidal
    suffix = trim(suffix)//'_n'//str_n

    call TGLFEP_scalefactor
  case default
    print *,'No code_in'
  end select

end subroutine TGLFEP_mainsub
