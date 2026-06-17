!===============================================================================
! Implementation of axisymmetric surfaces
!===============================================================================
module moose_axisurf
  use iso_fortran_env
  use moose_surface
  use moose_polygon2d
  implicit none
  private


  type, extends(hypersurf3d_patch), public :: axisurf
     class(polygon2d), allocatable :: P
     real(real64), allocatable :: v(:)
     logical :: vdef
     integer :: nv ! number of surface segments

     real(real64) :: phi0, dphi ! surface coordinate u = (phi - phi0) / dphi  mod  1
                                ! dphi = 2 * pi / symmetry

     contains
     procedure :: broadcast
     procedure :: free
     procedure :: savetxt, writenc
     procedure :: checksum => axisurf_checksum

     procedure :: eval
     procedure :: jac
     procedure :: interp, vphi, normal_vector

     procedure :: vcurve, vfallback, area

     procedure :: intersect, rzslice_intersect => intersect
     procedure :: ray_intersect => axisurf_ray_intersect
     procedure :: winding_number_phi
     procedure :: rmesh3d => make_rmesh3d
  end type axisurf


  interface axisurf
     procedure :: init
     procedure :: loadtxt
  end interface axisurf


  public :: &
     decode_header, readnc_axisurf

  contains
  !---------------------------------------------------------------------


! supplemental procedures:
  !-----------------------------------------------------------------------------
  subroutine decode_header(header, description, vlabel, phi0, symmetry)
  use moose_utils
  character(len=*),           intent(in)    :: header
  character(len=len(header)), intent(  out) :: description, vlabel
  real(real64),               intent(  out), optional :: phi0
  integer,                    intent(  out), optional :: symmetry

  character(len=len(header)) :: s, tmp
  integer :: i, n


  vlabel = ""
  if (present(phi0)) phi0 = 0.d0
  if (present(symmetry)) symmetry = 1
  n = nsubstrings(header)
  if (n <= 1) then
     description = lstrip(header, "# ")
     vlabel = ""
     return
  endif


  description = lstrip(substring(header, 1), '# ')
  do i=2,n
     s = substring(header, i)
     if (nsubstrings(s, '=') < 2) cycle

     if (substring(s, 1, '=') == "vlabel") then
        vlabel = lstrip(rstrip(substring(s, 2, '='), '" '), '" ')
     endif

     if (present(phi0)  .and.  substring(s, 1, '=') == "phi0") then
        tmp = substring(s, 2, '=')
        read (tmp, *) phi0
     endif

     if (present(symmetry)  .and.  substring(s, 1, '=') == "symmetry") then
        tmp = substring(s, 2, '=')
        read (tmp, *) symmetry
     endif
  enddo

  end subroutine decode_header
  !-----------------------------------------------------------------------------


! constructors:
  !---------------------------------------------------------------------
  function init(P, phi0, symmetry, metadata, v) result(this)
  use moose_error
  use moose_dict
  use moose_math, only: pi2
  class(polygon2d), intent(in) :: P
  real(real64),     intent(in), optional :: phi0, v(:)
  integer,          intent(in), optional :: symmetry
  type(dict),       intent(in), optional :: metadata
  type(axisurf)                :: this

  real(real64) :: vmin, vmax
  logical      :: bound


  this%vdef = present(v)
  if (this%vdef) then
     vmin = minval(v)
     vmax = maxval(v)
  else
     vmin = 0.d0
     vmax = P%length()
  endif


  bound = P%is_closed()
  call init_surface(this, [0.d0, pi2], [vmin, vmax], [.true., bound], metadata)
  allocate (this%P, source=P)
  this%nv = P%segments()
  this%phi0 = 0.d0;   if (present(phi0)) this%phi0 = phi0
  this%symmetry = 1;   if (present(symmetry)) this%symmetry = symmetry
  this%dphi = pi2 / this%symmetry
  if (this%vdef) then
     if (size(v) /= this%nv+1) call ERROR("size(v) /= size(P)", "init_axisurf")
     allocate (this%v(0:this%nv), source = v)
  endif

  end function init
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function loadtxt(filename, units, convert_units) result(this)
  !
  ! load surface contour in given units from text file (option: convert to new units)
  !
  use moose_dict
  use moose_table, only: table
  use moose_units, only: length_scale
  use moose_utils, only: user_option
  character(len=*), intent(in) :: filename, units
  character(len=*), intent(in), optional :: convert_units
  type(axisurf)                :: this


  character(len=256) :: header, description, vlabel
  type(polygon2d) :: P
  type(table) :: T
  type(dict) :: metadata
  real(real64) :: scale_factor, phi0
  integer :: iu, symmetry


  ! set scale factor if units need to be converted
  if (present(convert_units)) then
     scale_factor = length_scale(units) / length_scale(convert_units)
     call metadata%set("units", convert_units)
  else
     scale_factor = 1.d0
     call metadata%set("units", units)
  endif


  ! read text file
  open  (newunit=iu, file=filename)
  read  (iu, '(a)') header
  close (iu)
  T = table(filename, transposed=.true.)
  P = polygon2d(T%values(1:2,:) * scale_factor)

  if (header(1:1) == "#") then
     call decode_header(header, description, vlabel, phi0, symmetry)
     if (description /= "") call metadata%set("description", description)
     if (vlabel /= "") call metadata%set("vlabel", vlabel)
  else
     call metadata%set("description", trim(filename))
     vlabel = ""
     phi0 = 0.d0
     symmetry = 1
  endif


  ! initialize axisurf
  if (T%rows() == 3) then
     this = init(P, phi0, symmetry, metadata, T%values(3,:))
  else
     this = init(P, phi0, symmetry, metadata)
  endif

  end function loadtxt
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function readnc_axisurf(grp, convert_units) result(this)
  !
  ! read from netcdf group (optional: convert to given units)
  !
  use moose_netcdf
  use moose_dict
  use moose_units, only: length_scale
  use moose_utils, only: user_option
  class(netcdf_dataset), intent(in) :: grp
  character(len=*),      intent(in), optional :: convert_units
  type(axisurf)                     :: this

  type(dict) :: metadata
  character(len=128) :: units
  real(real64), allocatable :: x(:,:), v(:)
  real(real64) :: phi0
  integer :: nv, symmetry


  ! read metadata
  metadata = readnc_dict(grp)
  call metadata%pop("phi0", phi0, 0.d0)
  call metadata%pop("symmetry", symmetry, 1)
  units = metadata%get("units", "m")


  ! read data
  nv = grp%dim("nv")
  allocate (x(2, nv), v(nv))
  call grp%get_var("rz", x)
  call grp%get_var("v", v)
  if (present(convert_units)) then
     x = x * length_scale(units) / length_scale(convert_units)
     call this%metadata%set("units", convert_units)
  endif


  ! initialize axisurf
  this = axisurf(polygon2d(x), phi0, symmetry, metadata, v)

  end function readnc_axisurf
  !---------------------------------------------------------------------


! type-bound procedures:
  !---------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(axisurf), intent(inout) :: this


  call this%surface_broadcast()
  if (rank > 0) allocate (this%P)
  call this%P%broadcast()
  call proc(0)%broadcast(this%nv)
  call proc(0)%broadcast(this%vdef)
  if (this%vdef) call proc(0)%broadcast_allocatable(this%v)
  call proc(0)%broadcast(this%phi0)
  call proc(0)%broadcast(this%dphi)
  call proc(0)%broadcast(this%symmetry)

  end subroutine broadcast
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine free(this)
  class(axisurf), intent(inout) :: this


  call this%P%free()
  deallocate (this%P)
  if (this%vdef) deallocate (this%v)
  call this%surface_free()

  end subroutine free
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine savetxt(this, filename)
  use moose_utils, only: str
  class(axisurf),   intent(in) :: this
  character(len=*), intent(in) :: filename

  character(len=256) :: header
  integer :: i, iu


  open  (newunit=iu, file=filename, action="write")
  header = "# " // this%description()
  if (this%metadata%has_key("vlabel")) header = trim(header) // '; vlabel = "' // this%vlabel() // '"'
  if (this%phi0 /= 0.d0) header = trim(header)//'; phi0 = '//str(this%phi0)
  if (this%symmetry /= 1) header = trim(header)//'; symmetry = '//str(this%symmetry)
  write (iu, '(a)') trim(header)
  if (this%vdef) then
     do i=0,this%nv
        write (iu, *) this%P%node(i), this%v(i)
     enddo
  else
     write (iu, '(dt)') this%P
  endif
  close (iu)
 1000 format("# ",a)
 1001 format("# ",a,'; vlabel = "',a,'"')

  end subroutine savetxt
  !---------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, ncfile)
  use moose_netcdf
  use moose_math, only: pi
  class(axisurf),        intent(in) :: this
  class(netcdf_dataset), intent(in) :: ncfile

  integer :: ndim, nv


  call ncfile%def_dim("dim_0002", 2, ndim)
  call ncfile%def_dim("nv", this%nv + 1, nv)

  if (this%phi0 /= 0.d0) call ncfile%put_att("phi0", this%phi0)
  if (this%symmetry /= 1) call ncfile%put_att("symmetry", this%symmetry)
  call this%metadata%writenc(ncfile)

  call ncfile%def_var("v",    NF90_DOUBLE, [nv])
  call ncfile%def_var("rz",   NF90_DOUBLE, [ndim, nv])
  call ncfile%enddef()

  call ncfile%put_var("rz",   transpose(this%P%nodes()))
  if (this%vdef) then
     call ncfile%put_var("v", this%v)
  else
     call ncfile%put_var("v", this%vfallback())
  endif

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function axisurf_checksum(this)
  use moose_utils, only: checksum
  class(axisurf), intent(in) :: this
  integer                    :: axisurf_checksum


  axisurf_checksum = checksum([this%phi0, 1.d0*this%symmetry, this%P%nodes()])

  end function axisurf_checksum
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function eval(this, x) result(x3d)
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: x(this%ndim)
  real(real64)               :: x3d(this%mdim)


  !x3d(1:2) = ...
  x3d(3)   = x(1)
  ! @todo: implement sampling along polygon
  write (6, *) "ERROR: axisurf%eval not implemented yet!"
  stop

  end function eval
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function jac(this, x)
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: x(this%ndim)
  real(real64)               :: jac(this%mdim, this%ndim)


  ! @todo: implement sampling along polygon
  write (6, *) "ERROR: axisurf%jac not implemented yet!"
  stop

  end function jac
  !---------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp(this, u) result(x)
  !
  ! interpolate (r, z, phi) at u
  !
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: u(2)
  real(real64)               :: x(3)

  integer :: j


  j = int(u(2))
  x(1:2) = this%P%node(j) + (this%P%node(j+1) - this%P%node(j)) * (u(2) - j)
  x(3) = this%phi0 + u(1) * this%dphi

  end function interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vphi(this, u)
  !
  ! convert u to (v, phi [deg])
  !
  use moose_math, only: pi
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: u(2)
  real(real64)               :: vphi(2)

  integer :: j


  if (this%vdef) then
     j = int(u(2))
     vphi(1) = this%v(j) + (this%v(j+1) - this%v(j)) * (u(2) - j)
  else
     vphi(1) = u(2)
  endif
  vphi(2) = (this%phi0 + u(1) * this%dphi) * 180 / pi

  end function vphi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function normal_vector(this, u) result(v)
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: u(2)
  real(real64)               :: v(3)


  v(1:2) = this%p%normal(int(u(2)))
  v(3) = 0.d0

  end function normal_vector
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function vcurve(this, phi)
  !
  ! generate curve along v-coordinate
  !
  use moose_curve, only: curve
  use moose_bspline_curve, only: bspline_polygon
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: phi
  class(curve), allocatable   :: vcurve


  vcurve = bspline_polygon(this%P)

  end function vcurve
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function vfallback(this)
  class(axisurf), intent(in) :: this
  real(real64)               :: vfallback(0:this%nv)


  vfallback = this%P%accumulated_lengths()

  end function vfallback
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function area(this, j)
  !
  ! surface area of segment *j*
  !
  use moose_math, only: pi2
  class(axisurf), intent(in) :: this
  integer,        intent(in) :: j
  real(real64)               :: area

  real(real64) :: x1(2), x2(2), r0


  x1 = this%P%node(j)
  x2 = this%P%node(j+1)
  r0 = (x1(1) + x2(1)) / 2

  area = sqrt(sum((x2-x1)**2)) * pi2 * r0

  end function area
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function intersect(this, p1, p2, px, t, u)
  use moose_math, only: pi2
  class(axisurf), intent(in)  :: this
  real(real64),   intent(in)  :: p1(this%mdim), p2(this%mdim)
  real(real64),   intent(out) :: px(this%mdim), t, u(this%ndim)
  logical                     :: intersect

  real(real64) :: r, s
  integer      :: n


  intersect = .false.
  call this%P%intersect(p1(1:2), p2(1:2), XSECT_SEGMENT, px(1:2), t, s, n)
  if (n == -1) return


  ! calculate toroidal position of intersection
  px(3) = p1(3) + t * (p2(3) - p1(3))
  r     = (px(3) - this%phi0) / this%dphi
  u(1)  = r - floor(r)
  u(2)  = n + s
  intersect = .true.

  end function intersect
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function axisurf_ray_intersect(this, p0, v, px, t, n)
  !
  ! check for intersection of ray (p0 + t * v) with axisurf
  !
  use moose_math, only: pi2
  class(axisurf), intent(in)  :: this
  real(real64),   intent(in)  :: p0(this%mdim), v(this%mdim)
  real(real64),   intent(out) :: px(this%mdim), t
  integer,        intent(out) :: n
  logical                     :: axisurf_ray_intersect

  real(real64) :: pxi(3), ti, xa(2), xb(2), dx(2)
  integer :: i


  axisurf_ray_intersect = .false.
  t = huge(1.d0)
  xa = this%P%node(0)
  do i=1,this%P%nnodes()-1
     xb = this%P%node(i)
     dx = xb - xa

     if (ray_intersect(p0, v, xa, dx, pxi, ti)) then
        if (ti < t) then
           axisurf_ray_intersect = .true.
           n = i
           t = ti
           px = pxi
        endif
     endif

     xa = xb
  enddo

  end function axisurf_ray_intersect
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  pure function winding_number_phi(this, p) result(wn)
  class(axisurf), intent(in) :: this
  real(real64),   intent(in) :: p(3)
  integer                    :: wn


  wn = this%P%winding_number(p(1:2))

  end function winding_number_phi
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function make_rmesh3d(this, nphi) result(grid)
  use moose_math,  only: linspace, pi
  use moose_grids, only: rmesh3d
  class(axisurf), intent(in) :: this
  integer,        intent(in) :: nphi
  type(rmesh3d)              :: grid

  real(real64) :: phi(nphi)


  phi = linspace(this%phi0, this%phi0 + this%dphi, nphi) * 180 / pi
  if (this%vdef) then
     grid = rmesh3d(phi, this%v, transpose(this%P%nodes()), &
        "Toroidal Angle [deg]", this%vlabel(), "R", "Z")
  else
     grid = rmesh3d(phi, this%vfallback(), transpose(this%P%nodes()), &
        "Toroidal Angle [deg]", "Arclength along surface contour", "R", "Z")
  endif

  end function make_rmesh3d
  !---------------------------------------------------------------------


! module procedures:
  !---------------------------------------------------------------------
  function ray_intersect(p0, v, xa, dx, px, t)
  !
  ! check for intersection of ray (p0 + t * v) with axisymmetric
  ! surface element xa -> xa + dx
  !
  use moose_algorithms, only: solve_quadratic_equation
  real(real64),   intent(in   ) :: p0(3), v(3), xa(2), dx(2)
  real(real64),   intent(  out) :: px(3), t
  logical                       :: ray_intersect

  real(real64), parameter :: eps = 1.d-7

  real(real64) :: a, b, c, det, drdz, r2, r_zdrdz, vdrdz, s, t12(2)
  integer :: n


  ray_intersect = .false.

  ! planar surface element
  if (abs(dx(2)) < eps) then
     ! ray is parallel to plane
     if (abs(v(3)) < eps) return

     ! intersection is behind ray
     t = (xa(2) - p0(3)) / v(3)
     if (t < eps) return

     ! ray intersects plane
     px(1:2) = p0(1:2) + v(1:2) * t
     px(3) = xa(2)
     r2 = px(1)**2 + px(2)**2
     s = (sqrt(r2) - xa(1)) / dx(1)
     if (s < 0.d0  .or.  s > 1.d0) return

     ! ray intersects surface element
     ray_intersect = .true.


  ! non-planar surface element
  else
     ! prepare parameters for quadratic equation
     drdz = dx(1) / dx(2)
     vdrdz = v(3) * drdz
     r_zdrdz = xa(1) + (p0(3) - xa(2)) * drdz

     ! set up quadratic equation
     a = v(1)**2 + v(2)**2 - vdrdz**2
     b = 2 * (p0(1) * v(1) + p0(2) * v(2) - r_zdrdz * vdrdz)
     c = p0(1)**2 + p0(2)**2 - r_zdrdz**2

     ! solve quadratic equation
     det = b**2 - 4 * a * c
     if (det <= 0.d0) return
     call solve_quadratic_equation(a, b, c, 0.d0, huge(1.d0), t12, n)
     if (n == 0) return
     t = t12(1)
     px = p0 + v * t
     s = (px(3) - xa(2)) / dx(2)
     if (s < 0.d0  .or.  s > 1.d0) return

     ! ray intersects surface element
     ray_intersect = .true.
  endif

  end function ray_intersect
  !---------------------------------------------------------------------

end module moose_axisurf
