module flare_m3dc1
  use iso_fortran_env
  use moose_error
  use flare_bfield
  use flare_equi2d
  implicit none
  private


  ! M3D-C1 magnetic field ......................................................
  type, extends(magnetic_field), public :: m3dc1_bfield
     contains
     procedure :: eval, jac
  end type m3dc1_bfield


  interface m3dc1_bfield
     procedure :: load_m3dc1_bfield
  end interface
  ! m3dc1_bfield ...............................................................



  ! M3D-C1 equilibrium field ...................................................
  type, extends(equi2d), public :: m3dc1_equi2d
     contains
     procedure :: FpsiN, FdF
  end type m3dc1_equi2d
  ! m3dc1_equi2d ...............................................................



  public :: &
     load_m3dc1_equi2d

  contains
  !-----------------------------------------------------------------------------


! type m3dc1_bfield ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function load_m3dc1_bfield(filename, timeslice, amplitude, phase) result(this)
  use moose_utils, only: basename
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: timeslice
  real(real64),     intent(in), optional :: amplitude, phase
  type(m3dc1_bfield)           :: this


  call ERROR("no fusion-io support")

  end function load_m3dc1_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(m3dc1_bfield), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: B(this%mdim)


  call ERROR("no fusion-io support")

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(m3dc1_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: jac(this%mdim, this%ndim)


  call ERROR("no fusion-io support")

  end function jac
  !-----------------------------------------------------------------------------
! type m3dc1_bfield ============================================================



! type m3dc1_equi2d ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function load_m3dc1_equi2d(filename, factor) result(this)
  use moose_geometry, only: hypersurf2d
  use moose_utils,    only: dirname
  character(len=*),  intent(in) :: filename
  real(real64),      intent(in), optional :: factor
  type(m3dc1_equi2d)            :: this


  call ERROR("no fusion-io support")

  end function load_m3dc1_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  function FpsiN(this, psiN) result(F)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: psiN
  real(real64)                    :: F

  real(real64) :: B(3)


  call ERROR("no fusion-io support")

  end function FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function FdF(this, r)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: FdF(0:1)


  call ERROR("no fusion-io support")

  end function FdF
  !-----------------------------------------------------------------------------
! type m3dc1_equi2d ============================================================

end module flare_m3dc1
