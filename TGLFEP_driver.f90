!------------------------------------------------------------------
! TGLFEP_driver.f90
! 7.13.2016
! PURPOSE:
!  TGLFEP_driver calls TGLFEP_mainsub
!------------------------------------------------------------------

program TGLFEP_driver

  use mpi
  use TGLFEP_interface

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------
   
  integer :: id, np, color, key, ierr, STATUS(MPI_STATUS_SIZE)
  integer :: TGLFEP_COMM_IN
  logical :: multiple_flag = .true.
  integer :: scan_mode_flag = 1 !1 = radial, 2 = safety factor q, 3 = shear s, 4 = other scan
  integer :: scan_n = 25 !just scan one at one time. !q_scan_n = 8, shat_scan_n = 9
  character(2) :: str_r
  real :: q_data(8)
  data q_data /0.5, 1.0, 1.5, 2.0, 2.5/
  !data q_data /0.55, 0.625, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0/
  !data q_data /0.875, 1.125, 1.375, 1.625, 1.875/
  real :: shat_data(9)
  data shat_data /-1.0, -0.5, -0.2, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5/

  !For DIII-D NBI case with GYRO inputs
  real :: factor(17)
  data factor /2., 2., 1., 1., 1., 1., 1., 1., 1., 1., &
               1., 1., 2., 4., 6., 8., 12./
  real :: width(17)
  data width /0.99,0.86,0.87,0.85,0.82,0.81,0.83,0.87,0.98,1.09, &
              1.12,0.90,0.91,0.84,0.90,0.83,0.85/

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

  call EPtran_mainsub

  mode_flag_in = 2
  factor_in = 1.0
  !width_in = 1.87

  !ir = 9 !r/a = 0.5

  if(multiple_flag) then

    key = id / (scan_n)
    color = id - key*(scan_n)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,TGLFEP_COMM_IN,ierr)

    select case(scan_mode_flag)
    case(1)
      ir = color + 11
      !factor_in = factor(ir)
      !width_in = width(ir)
      str_r = achar(ir/10+iachar("0"))//achar(mod(ir,10)+iachar("0"))
      write(suffix,'(A2,A2)')'_r',str_r
    case(2)
      q_factor = q_data(color+1)
      write(suffix,'(A2,F4.2)')'_q',q_factor
    case(3)
      shat_factor = shat_data(color+1)
      write(suffix,'(A2,SP,F4.1)')'_s',shat_factor
    case(4)
      !scan_factor = real(color)*0.5
      scan_factor = real(color+1)/10.0
      write(suffix,'(A2,F3.1)')'_c',scan_factor
    end select

  else
    TGLFEP_COMM_IN = MPI_COMM_WORLD
  endif

  call TGLFEP_mainsub(TGLFEP_COMM_IN)

  call MPI_FINALIZE(ierr)

end program TGLFEP_driver
