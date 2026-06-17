module moose_geometry
  use moose_polygon
  use moose_polygon2d
  use moose_quad
  use moose_curve
  use moose_bspline_curve
  use moose_crspline
  use moose_ellipse
  use moose_interp_curve
  use moose_fourier_curve
  use moose_surface
  use moose_axisurf
  use moose_fourier_surface
  use moose_interp_surface
  use moose_rzplane
  use moose_torosurf
  use moose_trisurf
  use moose_hypersurface
  use moose_hypermesh3d
  use moose_contours
  implicit none


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_curve(instructions) result(C)
  use moose_error
  use moose_dict
  use moose_analysis, only: INTERP_LINEAR
  use moose_utils,    only: split
  use moose_r3grid,   only: length_scale
  character(len=*), intent(in) :: instructions
  class(curve), allocatable :: C

  type(dict)         :: metadata
  type(polygon)      :: P
  character(len=256) :: dtype, filename, units
  integer            :: iu


  call split(instructions, filename, units, set=':', default='m')
  open  (newunit=iu, file=filename, action="read")
  metadata = readtxt_dict(iu)


  ! determine curve type from header
  if (metadata%has_key("TYPE")) then
     dtype = metadata%get("TYPE")
     select case(dtype)
     case('bspline_curve')
        allocate (C, source=readtxt_bspline_curve(iu, metadata, length_scale(units)))

     case('fourier_curve')
        allocate (C, source=readtxt_fourier_curve(iu, metadata, length_scale(units)))

     case('interp_curve')
        allocate (C, source=readtxt_interp_curve(iu, metadata, length_scale(units)))

     case('polygon', 'polygon2d')

     case default
        call ERROR("unkown curve type '"//trim(dtype)//"'", "loadtxt_curve")
     end select
  endif
  close (iu)


  ! fallback: linear interpolation
  if (.not.allocated(C)) then
     P = polygon(filename, length_scale(units))
     allocate (C, source=interp_polygon(P))
  endif

  end function loadtxt_curve
  !-----------------------------------------------------------------------------

end module moose_geometry
