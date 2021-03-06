!------------------------------------------------------------------
! TGLFEP_driver.f90
! 7.13.2016
! PURPOSE:
!  TGLFEP_driver calls TGLFEP_mainsub
!------------------------------------------------------------------

program TGLFEP_driver

  use mpi
  use TGLFEP_interface
  use TGLFEP_profile

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------
   
  integer :: id, np, color, key, ierr, STATUS(MPI_STATUS_SIZE)
  integer :: TGLFEP_COMM_IN

  logical :: iexist, factor_in_profile
  integer :: scan_n = 1
  integer :: i,ii,irs
  real,allocatable,dimension(:) :: factor, width, kymark_out, SFmin
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

    read(33,*) scan_n
    read(33,*) irs

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
  key = id / (scan_n)
  color = id - key*(scan_n)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,TGLFEP_COMM_IN,ierr)

  ir = color + irs
  str_r = achar(ir/10+iachar("0"))//achar(mod(ir,10)+iachar("0"))
  write(suffix,'(A2,A2)')'_r',str_r

  factor_in = factor(color+1)
  if(width_in_flag) width_in = width(color+1)

  call TGLFEP_mainsub(TGLFEP_COMM_IN)

  !write output
  if(.not. width_in_flag) then !get width_out and kymark_out
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

    write(33,*) '--------------------------------------------------------------'
    write(33,*) 'factor_in_profile = ',factor_in_profile
    if(factor_in_profile) then
      do i = 1,scan_n
        write(33,'(F7.2)') factor(i)
      enddo
    else
      write(33,*) 'factor_in = ',factor(1)
    endif

    write(33,*) 'width_in_flag = ',width_in_flag
    if(.not. width_in_flag) write(33,*) 'width_min = ',width_min,'width_max = ',width_max

    close(33)

  endif

  if(process_in .eq. 4) then !print out 'density threshold'
    if(id .eq. 0) then
      open(unit=33,file='out.TGLFEP',status='old',position='append')
      write(33,*) '**************************************************************'
      write(33,*) '************** The critical EP density gradient **************'
      write(33,*) '**************************************************************'

      allocate(SFmin(scan_n))

      if(threshold_flag .eq. 0) then

        SFmin(1) = factor_in
        do i = 1,scan_n-1
          call MPI_RECV(factor_in,1,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,STATUS,ierr)
          SFmin(i+1) = factor_in
        enddo

        write(33,*) '--------------------------------------------------------------'
        write(33,*) 'SFmin'
        write(33,10) SFmin

        write(33,*) '--------------------------------------------------------------'
        write(33,*) 'The EP density threshold n_EP/n_e (%) for gamma_AE = 0'
        do i = 1,scan_n
          write(33,10) SFmin(i)*as(irs+i-1,is)*100. !percent
        enddo

        write(33,*) '--------------------------------------------------------------'
        write(33,*) 'The EP beta crit (%) = beta_e*(n_EP_th/n_e)*(T_EP/T_e)'
        do i = 1,scan_n
          if(geometry_flag .eq. 0) then
            write(33,10) SFmin(i)*betae(irs+i-1)*100.*as(irs+i-1,is)*taus(irs+i-1,is) !percent
          else
            write(33,10) SFmin(i)*betae(irs+i-1)*100.*as(irs+i-1,is)*taus(irs+i-1,is)*kappa(irs+i-1)**2 !percent
          endif
        enddo

        deallocate(SFmin)

        write(33,*) '--------------------------------------------------------------'
        write(33,*) 'The scale factor for EP density threshold at each n:'
        write(33,*) (factor_out(ii),ii=1,nn)
        do i = 1,scan_n-1
          call MPI_RECV(factor_out,nn,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
          write(33,*) (factor_out(ii),ii=1,nn)
        enddo

      else

        write(33,*) '--------------------------------------------------------------'
        write(33,*) 'a/Ln_EP multiplied by: 0.1, 0.2, ..., 1.5'
        write(33,*) 'The scale factor for EP density threshold:'
        write(33,*) (factor_out(ii),ii=1,nn)
        do i = 1,scan_n-1
          call MPI_RECV(factor_out,nn,MPI_DOUBLE_PRECISION,i,i+scan_n,MPI_COMM_WORLD,STATUS,ierr)
          write(33,*) (factor_out(ii),ii=1,nn)
        enddo

      endif

      close(33)
    else if(id .lt. scan_n) then
      if(threshold_flag .eq. 0) then
        call MPI_SEND(factor_in,1,MPI_DOUBLE_PRECISION,0,id,MPI_COMM_WORLD,STATUS,ierr)
      endif

      call MPI_SEND(factor_out,nn,MPI_DOUBLE_PRECISION,0,id+scan_n,MPI_COMM_WORLD,STATUS,ierr)
    endif
  endif

  if(process_in .eq. 4) deallocate(factor_out)
  deallocate(width)
  deallocate(kymark_out)
  deallocate(factor)

  ! if(id .eq. 0) call dump_profile

  call deallocate_profile

  call MPI_FINALIZE(ierr)

10 format(F8.4)

end program TGLFEP_driver
