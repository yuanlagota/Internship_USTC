module moose_units
  use iso_fortran_env
  implicit none


  character(len=*), parameter :: &
     KILOMETER  = "km", &
     METER      = "m", &
     CENTIMETER = "cm", &
     MILLIMETER = "mm", &
     MICROMETER = "μm", &
     NANOMETER  = "nm", &
     PICOMETER  = "pm", &
     FEMTOMETER = "fm"


  character(len=*), parameter :: &
     DEGREE     = "deg", &
     RADIAN     = "rad"



  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function length_scale(units)
  character(len=*), intent(in) :: units
  real(real64)                 :: length_scale


  select case(units)
  case(KILOMETER)
     length_scale = 1.d3
  case(METER)
     length_scale = 1.d0
  case(CENTIMETER)
     length_scale = 1.d-2
  case(MILLIMETER)
     length_scale = 1.d-3
  case(MICROMETER)
     length_scale = 1.d-6
  case(NANOMETER)
     length_scale = 1.d-9
  case(PICOMETER)
     length_scale = 1.d-12
  case(FEMTOMETER)
     length_scale = 1.d-15
  case default
     write (6, 9000) trim(units);   stop
  end select
 9000 format("ERROR: invalid units ",a,"!")

  end function length_scale
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine assert_angular_units(angular_units, procedure_name)
  character(len=*), intent(in) :: angular_units, procedure_name


  select case(angular_units)
  case(DEGREE, RADIAN)
  case default
     call ANGULAR_UNITS_ERROR(angular_units, procedure_name)
  end select

  end subroutine assert_angular_units
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine ANGULAR_UNITS_ERROR(angular_units, procedure_name)
  use moose_error
  character(len=*), intent(in) :: angular_units, procedure_name


  call ERROR("invalid angular units '"//trim(angular_units)//"'")

  end subroutine ANGULAR_UNITS_ERROR
  !-----------------------------------------------------------------------------

end module moose_units
