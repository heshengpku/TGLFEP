subroutine TGLFEP_tglf_map 

  implicit none
  integer :: model_in = 0

  select case(model_in)
  case(0)
    call EPstd_tglf_map
  case(1)
    call EPGYRO_tglf_map
  case(2)
    call EPtran_tglf_map
  case default
    print *,'No model_in'
  end select

end subroutine TGLFEP_tglf_map
