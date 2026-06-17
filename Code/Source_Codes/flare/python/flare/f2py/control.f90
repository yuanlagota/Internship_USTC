module control
  implicit none

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init(greeting)
  use flare_control, only: flare_init
  logical, intent(in) :: greeting


  call flare_init(standalone=.false., greeting=greeting)

  end subroutine init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function set_parameter(key, expr) result(istat)
  use moose_utils, only: startswith
  character(len=*), intent(in) :: key, expr
  integer                      :: istat

  character(len=len(key)) :: param


  istat = 1
  if     (startswith(key, "bspline3d.", rest=param)) then
     istat = set_bspline3d_parameter(param, expr)
  elseif  (startswith(key, "separatrix2d.", rest=param)) then
     istat = set_separatrix2d_parameter(param, expr)
  elseif (startswith(key, "fieldline.", rest=param)) then
     istat = set_fieldline_parameter(param, expr)
  elseif (startswith(key, "fluxsurf2d.", rest=param)) then
     istat = set_fluxsurf2d_parameter(param, expr)
  elseif (startswith(key, "fluxsurf3d.", rest=param)) then
     istat = set_fluxsurf3d_parameter(param, expr)
  elseif (startswith(key, "melnikov_function.", rest=param)) then
     istat = set_melnikov_function_parameter(param, expr)
  elseif (startswith(key, "screen_output.", rest=param)) then
     istat = set_screen_output_parameter(param, expr)
  elseif (startswith(key, "rpath2d.", rest=param)) then
     istat = set_rpath2d_parameter(param, expr)
  elseif (startswith(key, "base_mesh.", rest=param)) then
     istat = set_base_mesh_parameter(param, expr)
  elseif (startswith(key, "mmesh_generator.", rest=param)) then
     istat = set_mmesh_generator_parameter(param, expr)
  elseif (startswith(key, "task.", rest=param)) then
     istat = set_task_parameter(param, expr)
  endif

  end function set_parameter
  !-----------------------------------------------------------------------------
  function set_bspline3d_parameter(key, expr) result(istat)
  use flare_bspline3d
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("bmax")
     read (expr, *) bmax
  case default
     istat = 2
  end select

  end function set_bspline3d_parameter
  !-----------------------------------------------------------------------------
  function set_separatrix2d_parameter(key, expr) result(istat)
  use flare_control
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("step_size")
     read (expr, *) separatrix2d_step_size
  case("offset")
     read (expr, *) separatrix2d_offset
  case("fx")
     read (expr, *) separatrix2d_fX
  case("epsabs")
     read (expr, *) separatrix2d_epsabs
  case("nmax")
     read (expr, *) separatrix2d_nmax
  case("alpha")
     read (expr, *) separatrix2d_alpha
  case default
     istat = 2
  end select

  end function set_separatrix2d_parameter
  !-----------------------------------------------------------------------------
  function set_fieldline_parameter(key, expr) result(istat)
  use flare_fieldline
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("hstart")
     read (expr, *) hstart
  case("hmin")
     read (expr, *) hmin
  case("hmax")
     read (expr, *) hmax
  case("epsabs")
     read (expr, *) epsabs
  case("epsabs_xsect")
     read (expr, *) epsabs_xsect
  case("step_type")
     step_type = expr
  case("edom")
     read (expr, *) edom
  case("diffusion")
     read (expr, *) diffusion
  case default
     istat = 2
  end select

  end function set_fieldline_parameter
  !-----------------------------------------------------------------------------
  function set_fluxsurf2d_parameter(key, expr) result(istat)
  use flare_fluxsurf2d
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("step_size")
     read (expr, *) step_size
  case("epsabs")
     read (expr, *) epsabs
  case default
     istat = 2
  end select

  end function set_fluxsurf2d_parameter
  !-----------------------------------------------------------------------------
  function set_fluxsurf3d_parameter(key, expr) result(istat)
  use flare_fluxsurf3d
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("npoints")
     read (expr, *) npoints
  case("nctrl")
     read (expr, *) nctrl
  case("k")
     read (expr, *) k
  case("fit_method")
     read (expr, *) fit_method
  case("eps")
     read (expr, *) eps
  case("lambda1")
     read (expr, *) lambda1
  case("lambda2")
     read (expr, *) lambda2
  case("knot_balancing")
     read (expr, *) knot_balancing
  case default
     istat = 2
  end select

  end function set_fluxsurf3d_parameter
  !-----------------------------------------------------------------------------
  function set_melnikov_function_parameter(key, expr) result(istat)
  use flare_melnikov_function
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("epsabs")
     read (expr, *) epsabs
  case("mstart")
     read (expr, *) mstart
  case("mmax")
     read (expr, *) mmax
  case default
     istat = 2
  end select

  end function set_melnikov_function_parameter
  !-----------------------------------------------------------------------------
  function set_screen_output_parameter(key, expr) result(istat)
  use flare_control
  character(len=*), intent(in) :: key, expr
  integer                      :: istat

  integer :: verbosity


  istat = 0
  select case(key)
  case("verbosity")
     read (expr, *) verbosity
     verbose = .false.;   very_verbose = .false.
     if (verbosity >= 1) verbose = .true.
     if (verbosity >= 2) very_verbose = .true.
     if (verbosity <  0) report = .false.
  case default
     istat = 2
  end select

  end function set_screen_output_parameter
  !-----------------------------------------------------------------------------
  function set_rpath2d_parameter(key, expr) result(istat)
  use flare_rpath2d
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("epsabs")
     read (expr, *) rpath2d_epsabs
  case("hstart")
     read (expr, *) rpath2d_hstart
  case("hmin")
     read (expr, *) rpath2d_hmin
  case("hmax")
     read (expr, *) rpath2d_hmax
  case("xoffset")
     read (expr, *) rpath2d_Xoffset
  case default
     istat = 2
  end select

  end function set_rpath2d_parameter
  !-----------------------------------------------------------------------------
  function set_base_mesh_parameter(key, expr) result(istat)
  use flare_mmesh_base_generator
  use moose_qmesh_generator
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("inner_boundary")
     read (expr, *) inner_boundary_method
  case("updown_symmetry_tolerance")
     read (expr, *) updown_symmetry_tolerance
  case("points_per_segment")
     read (expr, *) points_per_segment
  case("min_width")
     read (expr, *) min_width
  case("max_squeeze")
     read (expr, *) max_squeeze
  case default
     istat = 2
  end select

  end function set_base_mesh_parameter
  !-----------------------------------------------------------------------------
  function set_mmesh_generator_parameter(key, expr) result(istat)
  use flare_mmesh_generator
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  istat = 0
  select case(key)
  case("npoints")
     read (expr, *) npoints
  case("nctrl")
     read (expr, *) nctrl
  case("epsabs")
     read (expr, *) epsabs
  case("fit_method")
     read (expr, *) fit_method
  case("lambda1")
     read (expr, *) lambda1
  case("lambda2")
     read (expr, *) lambda2
  case default
     istat = 2
  end select

  end function set_mmesh_generator_parameter
  !-----------------------------------------------------------------------------
  function set_task_parameter(key, expr) result(istat)
  use flare_control
  character(len=*), intent(in) :: key, expr
  integer                      :: istat


  select case(key)
  case("diagnostic_mode")
     read (expr, *, iostat=istat) diagnostic_mode
  case default
     istat = 2
  end select

  end function set_task_parameter
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine reset_counter(key)
  use flare_control
  integer, intent(in) :: key


  if (key == 0) then
     counter(:) = 0
  else
     counter(key) = 0
  endif

  end subroutine reset_counter
  !-----------------------------------------------------------------------------
  function get_counter(key)
  use flare_control
  integer, intent(in) :: key
  integer(kind=8)     :: get_counter


  get_counter = counter(key)

  end function get_counter
  !-----------------------------------------------------------------------------

end module control
