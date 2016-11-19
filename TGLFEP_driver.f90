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

  logical :: iexist, scan_flag, factor_in_profile
  integer :: scan_parameter_flag = 0 !0 = radial, 1 = safety factor q, 2 = shear s, 
                                     !3 = EP density gradient scan (or other scan)
  integer :: scan_n = 1 !q scan = 13, shat scan = 9 for GA-std case
                        !radial scan for DIII-D NBI (GYRO inputs) = 17
                        !radial scan for EPtran = 50 (start at irs = 2, r/a = 0 excluded)
  integer :: i,ii,irs
  real,allocatable,dimension(:) :: scan_p, factor, width, kymark_out
  character(2) :: str_r

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

  !read input
  inquire(file='input.TGLFEP',exist=iexist)
  if(iexist) then
    open(unit=33,file='input.TGLFEP',status='old')
    read(33,'(I1)') process_in

    if(process_in .le. 1) read(33,'(I1)') mode_in
    if(process_in .eq. 4) read(33,'(I1)') threshold_flag

    read(33,'(I1)') ky_model

    read(33,*) scan_flag
    if(scan_flag) then
      read(33,'(I1)') scan_parameter_flag
      if(scan_parameter_flag .eq. 0) then
        read(33,*) scan_n
        read(33,*) irs
      else
        read(33,*) ir
        read(33,*) scan_n

        allocate(scan_p(scan_n))
        do i = 1,scan_n
          read(33,*) scan_p(i)
        enddo
      endif

      read(33,*) factor_in_profile
      allocate(factor(scan_n))
      if(factor_in_profile) then
        do i = 1,scan_n
          read(33,*) factor(i)
        enddo
      else
        read(33,*) factor_in
        factor(:) = factor_in
      endif

      read(33,*) width_in_flag
      allocate(width(scan_n))
      allocate(kymark_out(scan_n))
      if(width_in_flag) then
        do i = 1,scan_n
          read(33,*) width(i)
        enddo
      else
        width = 0.00
        read(33,*) width_min
        read(33,*) width_max
      endif
    else
      read(33,*) ir

      allocate(factor(1))

      read(33,*) factor_in
      factor(1) = factor_in

      read(33,*) width_in_flag
      if(width_in_flag) then
        read(33,*) width_in
      else
        read(33,*) width_min
        read(33,*) width_max
      endif
    endif

    close(33)
  else
    print *, 'input.TGLFEP file not found'
    stop
  endif

  call read_input_profile

  if(ky_model .eq. 0) then
    n_toroidal = 4
  else
    n_toroidal = 3
  endif

  if(process_in .eq. 4) then
    if(threshold_flag .eq. 0) then
      nn = 5
    else
      nn = 15
    endif
    allocate(factor_out(nn))
  endif

  !run mainsub
  if(scan_flag) then

    key = id / (scan_n)
    color = id - key*(scan_n)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,TGLFEP_COMM_IN,ierr)

    select case(scan_parameter_flag)
    case(0)
      ir = color + irs
      str_r = achar(ir/10+iachar("0"))//achar(mod(ir,10)+iachar("0"))
      write(suffix,'(A2,A2)')'_r',str_r
    case(1)
      q_factor = scan_p(color+1)
    case(2)
      shat_factor = scan_p(color+1)
    case(3)
      scan_factor = scan_p(color+1)
    end select

    if(scan_parameter_flag .gt. 0) then
      str_r = achar((color+1)/10+iachar("0"))//achar(mod(color+1,10)+iachar("0"))
      write(suffix,'(A2,A2)')'_s',str_r
    endif

    factor_in = factor(color+1)
    if(width_in_flag) width_in = width(color+1)

  else
    TGLFEP_COMM_IN = MPI_COMM_WORLD
  endif

  call TGLFEP_mainsub(TGLFEP_COMM_IN)

  !write output
  if(scan_flag .and. .not. width_in_flag) then !get width_out and kymark_out
    if(id .eq. 0) then
      width(1) = width_in
      kymark_out(1) = kymark
      do i = 1,scan_n-1
        call MPI_RECV(width(i+1),1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,STATUS,ierr)
        call MPI_RECV(kymark_out(i+1),1,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
      enddo
    else if(id .lt. scan_n) then
      call MPI_SEND(width_in,1,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,STATUS,ierr)
      call MPI_SEND(kymark,1,MPI_DOUBLE_PRECISION,0,id+scan_n,MPI_COMM_WORLD,STATUS,ierr)
    endif
  endif

  if(id .eq. 0) then

    open(unit=33,file='out.TGLFEP',status='replace')
    write(33,*) 'process_in = ',process_in

    if(process_in .le. 1) write(33,*) 'mode_in = ',mode_in
    if(process_in .eq. 4) write(33,*) 'threshold_flag = ',threshold_flag

    write(33,*) 'ky_model = ',ky_model

    write(33,*) '--------------------------------------------------------------'
    write(33,*) 'scan_flag = ',scan_flag
    if(scan_flag) then
      write(33,*) 'scan_parameter_flag = ',scan_parameter_flag
      if(scan_parameter_flag .eq. 0) then
        write(33,*) 'scan_n = ',scan_n
        write(33,*) 'irs = ',irs

        if(width_in_flag) then
          write(33,*) 'ir,  width'
          do i = 1,scan_n
            write(33,'(I3,F8.2)') irs+i-1,width(i)
          enddo
        else
          write(33,*) 'ir,  width,  kymark'
          do i = 1,scan_n
            write(33,'(I3,F8.2,F9.3)') irs+i-1,width(i),kymark_out(i)
          enddo
        endif

      else
        write(33,*) 'ir = ',ir
        write(33,*) 'scan_n = ',scan_n

        if(width_in_flag) then
          write(33,*) 'scan parameter values,  width'
          do i = 1,scan_n
            write(33,'(F6.3,F8.2)') scan_p(i),width(i)
          enddo
        else
          write(33,*) 'scan parameter values,  width,  kymark'
          do i = 1,scan_n
            write(33,'(F6.3,F8.2,F9.3,F8.2)') scan_p(i),width(i),kymark_out(i)
          enddo
        endif

      endif

      write(33,*) '--------------------------------------------------------------'
      write(33,*) 'factor_in_profile = ',factor_in_profile
      if(factor_in_profile) then
        do i = 1,scan_n
          write(33,'(F7.2)') factor(i)
        enddo
      else
        write(33,*) 'factor_in = ',factor(1)
      endif

    else
      write(33,*) 'ir = ',ir
      write(33,*) 'width = ',width_in
      write(33,*) 'kymark = ',kymark

      write(33,*) '--------------------------------------------------------------'
      write(33,*) 'factor_in = ',factor(1)
    endif

    write(33,*) 'width_in_flag = ',width_in_flag
    if(.not. width_in_flag) write(33,*) 'width_min = ',width_min,'width_max = ',width_max

  endif

  if(process_in .eq. 4) then !print out 'density threshold'
    if(id .eq. 0) then
      open(unit=33,file='out.density_threshold',status='replace')
      write(33,*)'scan_parameter_flag =',scan_parameter_flag,' scan_n =',scan_n

      write(33,*) 'SFmin'
      write(33,10) factor_in
      do i = 1,scan_n-1
        call MPI_RECV(factor_in,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,STATUS,ierr)
        write(33,10) factor_in
      enddo

      write(33,*) '--------------------------------------------------------------'
      write(33,*) 'At each n:'
      write(33,*) (factor_out(ii),ii=1,nn)
      do i = 1,scan_n-1
        call MPI_RECV(factor_out,nn,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
        write(33,*) (factor_out(ii),ii=1,nn)
      enddo

      close(33)
    else if(id .lt. scan_n) then
      call MPI_SEND(factor_in,1,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,STATUS,ierr)
      call MPI_SEND(factor_out,nn,MPI_DOUBLE_PRECISION,0,id+scan_n,MPI_COMM_WORLD,STATUS,ierr)
    endif
  endif

  if(process_in .eq. 4) deallocate(factor_out)
  if(scan_flag) then
    if(scan_parameter_flag .gt. 0) deallocate(scan_p)
    deallocate(width)
    deallocate(kymark_out)
  endif
  deallocate(factor)

  ! if(id .eq. 0) call dump_profile

  call deallocate_profile

  call MPI_FINALIZE(ierr)

10 format(F8.4)

end program TGLFEP_driver
