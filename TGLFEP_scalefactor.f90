subroutine TGLFEP_scalefactor

  use mpi
  use tglf_interface
  use tglf_pkg
  use TGLFEP_interface 

  implicit none
  real :: f0,f1

  call TGLFEP_ky_widthscan  !get the width_in

  f0 = 0.0
  f1 = factor_in

  do while(fmark .gt. freq_AE_upper)

    !f0 = factor_in
    f1 = 2.0*f1

    if(f1 .gt. factor_max) then
      f1 = factor_max
      factor_in = factor_max
      call TGLFEP_ky_widthscan  !get the width_in
      exit
    else
      factor_in = f1
      call TGLFEP_ky_widthscan  !get the width_in
    endif

  end do

  call TGLFEP_bisection

end subroutine TGLFEP_scalefactor
