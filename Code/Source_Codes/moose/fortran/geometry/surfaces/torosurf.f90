!===============================================================================
! Implementation of toroidal, non-axisymmetric surfaces
!===============================================================================
module moose_torosurf
  use iso_fortran_env
  use moose_surface
  implicit none
  private


  type, extends(hypersurf3d_patch), public :: torosurf
     ! number of surface segments in toroidal and poloidal direction
     integer :: nu, nv

     ! geometry
     real(real64), allocatable :: phi(:), v(:,:), rz(:,:,:)
     integer :: phi_order
     logical :: vdef

     ! precalculated geometry parameters for speedup of intersection
     real(real64), private :: dphiSym ! = 2 * pi / symmetry
     integer, private :: m0

     ! bounding box and shape coefficients for surface elements
     real(real64), private, dimension(:,:,:,:), allocatable :: bbox
     real(real64), private, dimension(:,:,:),   allocatable :: B, C, D
     real(real64), private, dimension(:),       allocatable :: dphi

     contains
     procedure :: reverse_phi_order
     procedure :: broadcast
     procedure :: free
     procedure :: savetxt, savenc, writenc
     procedure :: checksum => torosurf_checksum

     procedure :: toroidal_coordinates
     procedure :: varray
     procedure :: poloidal_coordinates
     procedure :: aux_eval
     procedure :: eval
     procedure :: jac
     procedure :: interp, vphi, normal_vector, is_closed

     procedure :: area, lhshift

     procedure :: intersect, rzslice_intersect
     procedure :: includes
     generic :: winding_number => winding_number_iphi
     procedure :: winding_number_phi, winding_number_iphi
     procedure :: vcurve
     procedure :: rzslice, polygon2d => polygon2d_slice
     procedure :: tpzmesh3d => make_tpzmesh3d
     procedure :: r3grid => make_r3grid
  end type torosurf


  interface torosurf
     procedure :: new
     procedure :: init
     procedure :: loadtxt
     procedure :: axisymmetric_torosurf
  end interface torosurf


  public :: &
     loadnc_torosurf, readnc_torosurf, setup_torosurf

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(nu, nv, symmetry, phi_order, vdef) result(this)
  integer, intent(in) :: nu, nv, symmetry
  integer, intent(in), optional :: phi_order
  logical, intent(in), optional :: vdef
  type(torosurf)      :: this


  this%nu   = nu
  this%nv   = nv
  this%symmetry = symmetry
  allocate (this%phi(0:nu), this%v(0:nv, 0:nu), source=0.d0)
  allocate (this%rz(2, 0:nv, 0:nu), source=0.d0)
  this%phi_order = 1;   if (present(phi_order)) this%phi_order = phi_order
  this%vdef = .false.;   if (present(vdef)) this%vdef = vdef

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine reverse_phi_order(this)
  use moose_algorithms
  class(torosurf), intent(inout) :: this

  real(real64), allocatable :: tmp(:,:)
  integer :: i


  ! reverse poloidal direction along with toroidal direction for same direction of normal vector
  this%phi = reverse_array(this%phi)
  allocate (tmp(2,0:this%nv))
  do i=0,this%nu/2
     tmp = this%rz(:,:,i)
     this%rz(:,:,i) = reverse_matrix(this%rz(:,:,this%nu-i), 2)
     this%rz(:,:,this%nu-i) = reverse_matrix(tmp, 2)
  enddo
  deallocate (tmp)
  this%phi_order = -this%phi_order

  end subroutine reverse_phi_order
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init(mesh, angle) result(this)
  use moose_error
  use moose_math,  only: pi
  use moose_grids, only: tpzmesh3d
  type(tpzmesh3d),  intent(in) :: mesh
  type(torosurf)               :: this
  character(len=*), intent(in), optional :: angle

  real(real64) :: angle_scale
  integer :: nu, nv


  angle_scale = 1.d0
  if (present(angle)) then
     select case (angle)
     case ('rad')
     case ('deg')
        angle_scale = pi / 180.d0
     case default
        call ERROR("invalid angle type")
     end select
  endif

  nu   = mesh%n(1) - 1
  nv   = mesh%n(2) - 1
  this = new(nu, nv, 1, vdef=.true.)
  this%phi       = mesh%domain%u * angle_scale
  this%v         = transpose(mesh%domain%v)
  this%rz(1,:,:) = transpose(mesh%x(1,:,:))
  this%rz(2,:,:) = transpose(mesh%x(2,:,:))
  this%phi_order = 1;   if (this%phi(nu) - this%phi(0) < 0.d0) call reverse_phi_order(this)
  call setup_torosurf(this)

  end function init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt(filename, units, convert_units, vfallback) result(this)
  !
  ! load surface in given units from text file (optional: convert to new units)
  !
  use moose_error
  use moose_math, only: linspace, pi
  use moose_units, only: length_scale
  use moose_utils, only: user_option
  use moose_axisurf, only: decode_header
  character(len=*), intent(in) :: filename, units
  character(len=*), intent(in), optional :: convert_units, vfallback
  type(torosurf)               :: this

  integer, parameter :: iu = 99

  character(len=128) :: header, line, name, vlabel
  real(real64) :: phi, x(2), tmp, R0, Z0
  integer :: i, ierr, j, nu, nv


  ! read header
  open  (iu, file=filename)
  read  (iu, '(a)') header
  do
     read  (iu, '(a)') line
     if (line /= "") exit
  enddo
  if (line(1:1) == "#") line = line(2:)
  read  (line, *) nu, nv, this%symmetry
  ! optional reference point
  R0 = 0.d0
  Z0 = 0.d0
  read  (line, *, end=10) nu, nv, this%symmetry, R0, Z0
 10 continue
  this%nu = nu - 1
  this%nv = nv - 1


  ! set name and v-label
  call decode_header(header, name, vlabel)
  if (name == "") name = trim(filename)
  call this%metadata%set("name", name)


  ! check layout
  if (nu == 1) then
     write (6, 9001);   stop
  endif
 9001 format("error: non-axisymmetric surface must have toroidal resolution nu > 1")


  ! @todo: check if surface is closed in toroidal & poloidal direction
  ! initialize variables
  allocate (this%phi(0:this%nu), this%v(0:this%nv, 0:this%nu))
  allocate (this%rz(2, 0:this%nv, 0:this%nu))
  this%phi_order = 1


  ! read data
  this%vdef = .true.
  do i=0,this%nu
     ! toroidal position of slice i
     read  (iu, *) phi;   this%phi(i) = phi / 180.d0 * pi

     ! RZ-outline of slice i
     do j=0,this%nv
        read  (iu, '(a)') line
        read  (line, *) this%rz(:,j,i)
        read  (line, *, iostat=ierr) this%rz(:,j,i), this%v(j,i);   if (ierr /= 0) this%vdef = .false.
        this%rz(1,j,i) = this%rz(1,j,i) + R0
        this%rz(2,j,i) = this%rz(2,j,i) + Z0
     enddo
  enddo
  close (iu)
  if (present(convert_units)) then
     this%rz = this%rz * length_scale(units) / length_scale(convert_units)
     call this%metadata%set("units", convert_units)
  else
     call this%metadata%set("units", units)
  endif


  ! fallback definition for v
  if (.not. this%vdef) then
     select case(user_option("index", vfallback))
     case ("index")
        this%v = spread(linspace(0.d0, 1.d0 * this%nv, this%nv+1), 2, this%nu+1)
        vlabel = "Node index"

     case ("arclength")
        call set_arclength_parametrization(this, .false.)
        vlabel = "Arc length [" // this%units() // "]"

     case default
        call ERROR("invalid choice vfallback = '"//vfallback//"'")
     end select
  endif
  if (vlabel /= "") call this%metadata%set("vlabel", vlabel)


  ! verify phi is in increasing order
  do i=1,this%nu-1
     if ((this%phi(i+1)-this%phi(i)) * (this%phi(i)-this%phi(i-1)) < 0.d0) then
        write (6, 9002);   stop
     endif
  enddo
 9002 format("error: phi must be monotonic")
  ! reverse direction (if necessary)
  if (this%phi(0) > this%phi(this%nu)) call reverse_phi_order(this)


  ! set up dependent coefficients
  call setup_torosurf(this)

  end function loadtxt
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  function loadnc_torosurf(filename, convert_units) result(this)
  !
  ! load torosurf from netcdf file (optional: convert units)
  !
  use moose_netcdf
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: convert_units
  type(torosurf)               :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  this = readnc_torosurf(nc, convert_units)
  call nc%close()

  end function loadnc_torosurf
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function readnc_torosurf(grp, convert_units) result(this)
  !
  ! read torosurf from netcdf group (optional: convert units)
  !
  use moose_netcdf
  use moose_dict
  use moose_math,  only: pi
  use moose_units, only: length_scale
  use moose_utils, only: user_option
  class(netcdf_dataset), intent(in) :: grp
  character(len=*),      intent(in), optional :: convert_units
  type(torosurf)                    :: this

  type(dict) :: metadata
  character(len=128) :: units
  integer :: ierr, nu, nv, symmetry


  ! read metadata
  metadata = readnc_dict(grp)
  if (metadata%has_key("symmetry")) then
     call metadata%pop("symmetry", symmetry)
  else
     ! legacy ...
     call metadata%pop("nsym", symmetry)
  endif
  units = metadata%get("units", "m")


  ! read data
  nu = grp%dim("nu") - 1
  nv = grp%dim("nv") - 1
  this = new(nu, nv, symmetry)
  call grp%get_var("phi", this%phi)
  call grp%get_var("v", this%v)
  call grp%get_var("rz", this%rz)
  this%metadata = metadata
  this%vdef = .true.

  this%phi = this%phi / 180.d0 * pi
  if (present(convert_units)) then
     this%rz = this%rz * length_scale(units) / length_scale(convert_units)
     call this%metadata%set("units", convert_units)
  endif

  if (this%phi(nu) - this%phi(0) < 0.d0) call reverse_phi_order(this)
  call setup_torosurf(this)

  end function readnc_torosurf
  !---------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function axisymmetric_torosurf(A, symmetry, nphi) result(this)
  !
  ! construct torosurf object from axisurf *A*
  !
  use moose_utils, only: user_option
  use moose_math,  only: linspace, pi2
  use moose_axisurf
  class(axisurf), intent(in) :: A
  integer,        intent(in), optional :: symmetry, nphi
  type(torosurf)             :: this

  integer :: i


  this = new(user_option(1, nphi), A%nv, user_option(1, symmetry))
  this%phi = linspace(0.d0, pi2 / this%symmetry, this%nu + 1)
  do i=0,A%nv
     this%rz(:,i,0) = A%P%node(i)
     this%rz(:,i,1:) = transpose(spread(this%rz(:,i,0), 1, this%nu))
  enddo
  ! TODO: vdef
  call set_arclength_parametrization(this, .true.)
  call setup_torosurf(this)

  end function axisymmetric_torosurf
  !-----------------------------------------------------------------------------


! supplemental constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine set_arclength_parametrization(this, normalized)
  !
  ! set poloidal coordinate from arc length
  !
  class(torosurf), intent(inout) :: this
  logical,         intent(in   ) :: normalized

  real(real64) :: dv
  integer :: i, j


  do i=0,this%nu
     this%v(0,i) = 0.d0
     do j=1,this%nv
        dv = sqrt(sum((this%rz(:,j,i)-this%rz(:,j-1,i))**2))
        this%v(j,i) = this%v(j-1,i) + dv
     enddo
     if (normalized) this%v(:,i) = this%v(:,i) / this%v(this%nv,i)
  enddo

  end subroutine set_arclength_parametrization
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(torosurf), intent(inout) :: this


  call proc(0)%broadcast(this%nu)
  call proc(0)%broadcast(this%nv)
  call proc(0)%broadcast(this%symmetry)
  call proc(0)%broadcast(this%vdef)
  call proc(0)%broadcast_allocatable(this%phi)
  call proc(0)%broadcast_allocatable(this%v)
  call proc(0)%broadcast_allocatable(this%rz)
  call this%metadata%broadcast()
  if (rank == 0) return

  call setup_torosurf(this)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(torosurf), intent(inout) :: this


  deallocate (this%phi, this%v, this%rz)
  deallocate (this%bbox, this%B, this%C, this%D, this%dphi)
  call this%surface_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savetxt(this, filename, iotype, scale)
  use moose_error, only: VALUE_ERROR
  use moose_math,  only: pi
  use moose_tpzmesh3d
  class(torosurf),  intent(in) :: this
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: iotype
  real(real64),     intent(in), optional :: scale

  integer, parameter :: &
     IOTYPE_DEFAULT = 0, &
     IOTYPE_LEGACY  = 1

  character(len=256) :: err, title
  type(tpzmesh3d) :: T
  real(real64)    :: x(2, 0:this%nu, 0:this%nv), xscale
  integer :: i


  ! select output type
  i = IOTYPE_DEFAULT
  if (present(iotype)) then
     select case (iotype)
     case ('', 'default')
        i = IOTYPE_DEFAULT

     case ('legacy')
        i = IOTYPE_LEGACY

     case default
        write (err, 9001) trim(iotype)
        call VALUE_ERROR(err)
     end select
  endif
 9001 format("invalid iotype ",a)


  ! scale factor
  xscale = 1.d0;   if (present(scale)) xscale = scale


  ! prepare data for output
  x(1,:,:) = transpose(this%rz(1,:,:)) * xscale
  x(2,:,:) = transpose(this%rz(2,:,:)) * xscale
  T = tpzmesh3d(transpose(this%v), this%phi / pi * 180.d0, x)


  ! write data
  select case (i)
  case (IOTYPE_DEFAULT)
     call T%savetxt(filename)
  case (IOTYPE_LEGACY)
     title = this%description()
     if (this%vdef .and. this%vlabel() /= "") write (title, 1001) trim(title), trim(this%vlabel())
     call T%write_legacy_format(filename, this%symmetry, title=trim(title), phi_order=this%phi_order, vdef=this%vdef)
  end select
 1001 format(a,'; vlabel = "',a,'"')

  end subroutine savetxt
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  class(torosurf),       intent(in) :: this
  character(len=*),      intent(in) :: filename

  type(netcdf_dataset) :: nc


  nc = netcdf_create(filename, "torosurf")
  call this%writenc(nc)
  call nc%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, ncfile)
  use moose_netcdf
  use moose_math, only: pi
  class(torosurf),       intent(in) :: this
  class(netcdf_dataset), intent(in) :: ncfile

  integer :: ndim, nu, nv


  call ncfile%def_dim("dim_0002", 2, ndim)
  call ncfile%def_dim("nu", this%nu + 1, nu)
  call ncfile%def_dim("nv", this%nv + 1, nv)
  call ncfile%put_att("symmetry",  this%symmetry)
  call this%metadata%writenc(ncfile)

  call ncfile%def_var("phi",  NF90_DOUBLE, [nu])
  call ncfile%def_var("v",    NF90_DOUBLE, [nv, nu])
  call ncfile%def_var("rz",   NF90_DOUBLE, [ndim, nv, nu])
  call ncfile%enddef()

  call ncfile%put_var("phi", this%phi * 180.d0 / pi)
  call ncfile%put_var("v",   this%v)
  call ncfile%put_var("rz",  this%rz)

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function torosurf_checksum(this)
  use moose_utils, only: checksum
  class(torosurf), intent(in) :: this
  integer                     :: torosurf_checksum


  torosurf_checksum = checksum([1.d0*this%symmetry, this%phi, this%rz])

  end function torosurf_checksum
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  subroutine setup_torosurf(this)
  !
  ! set up bounding boxes and shape coefficients for given geometry
  !
  use moose_math, only: pi2
  class(torosurf), intent(inout) :: this

  real(real64) :: vrange(2), phic
  integer :: i, j, k
  associate (rz => this%rz)


  vrange = [minval(this%v), maxval(this%v)]
  call init_surface(this, [this%phi(0), this%phi(this%nu)], vrange, [.false., .false.])
  ! set toroidal domain of symmetry
  phic = (this%phi(0) + this%phi(this%nu)) / 2
  this%dphiSym = pi2 / abs(this%symmetry)
  this%m0 = int(floor(phic / this%dphiSym))


  ! set up bounding boxes and shape coefficients
  allocate (this%bbox(2, 2, 0:this%nv-1, 0:this%nu-1))
  allocate (this%B(2, 0:this%nv-1, 0:this%nu-1))
  allocate (this%C(2, 0:this%nv-1, 0:this%nu-1))
  allocate (this%D(2, 0:this%nv-1, 0:this%nu-1))
  allocate (this%dphi(0:this%nu-1))
  do i=0,this%nu-1
     this%dphi(i) = this%phi(i+1) - this%phi(i)
     do j=0,this%nv-1
        this%bbox(1,1,j,i) = minval(rz(1,j:j+1,i:i+1))
        this%bbox(2,1,j,i) = maxval(rz(1,j:j+1,i:i+1))
        this%bbox(1,2,j,i) = minval(rz(2,j:j+1,i:i+1))
        this%bbox(2,2,j,i) = maxval(rz(2,j:j+1,i:i+1))

        do k=1,2
           this%B(3-k,j,i) = rz(k,j+1,i  ) - rz(k,j,i)
           this%C(k  ,j,i) = rz(k,j  ,i+1) - rz(k,j,i)
           this%D(3-k,j,i) = rz(k,j+1,i+1) - rz(k,j,i+1) - rz(k,j+1,i) + rz(k,j,i)
        enddo
        this%B(1,j,i) = - this%B(1,j,i)
        this%D(1,j,i) = - this%D(1,j,i)
     enddo
  enddo

  end associate
  end subroutine setup_torosurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine toroidal_coordinates(this, phi, it, t)
  !
  ! toroidal index *it* with this%phi(it) <= phi [rad] <= this%(it+1), or it = -1
  ! and local coordinate *t*
  !
  use moose_algorithms, only: binary_search_L
  class(torosurf), intent(in   ) :: this
  real(real64),    intent(in   ) :: phi
  integer,         intent(  out) :: it
  real(real64),    intent(  out) :: t


  if (phi < this%phi(0)  .or.  phi > this%phi(this%nu)) then
     it = -1;   t = -1.d0
     return
  endif


  if (phi == this%phi(0)) then
     it = 0
  else
     it = binary_search_L(this%phi, phi) - 1
  endif
  t = (phi - this%phi(it)) / (this%phi(it+1) - this%phi(it))

  end subroutine toroidal_coordinates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function varray(this, it, t)
  !
  ! interpolated v at (it, t)
  !
  class(torosurf), intent(in) :: this
  integer,         intent(in) :: it
  real(real64),    intent(in) :: t
  real(real64)                :: varray(0:this%nv)


  varray = (1.d0-t) * this%v(:,it)  +  t * this%v(:,it+1)

  end function varray
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine poloidal_coordinates(this, it, t, v, ip, p)
  !
  ! poloidal index *ip* with this%varray(it,t)(ip) <= v <= this%varray(it,t)(ip+1), or ip = -1
  ! and local coordinate *p*
  !
  use moose_algorithms
  class(torosurf), intent(in   ) :: this
  integer,         intent(in   ) :: it
  real(real64),    intent(in   ) :: t, v
  integer,         intent(  out) :: ip
  real(real64),    intent(  out) :: p

  real(real64) :: varray(0:this%nv)


  varray = this%varray(it, t)
  if (v < varray(0)  .or.  v > varray(this%nv)) then
     ip = -1
     return
  endif


  if (v == varray(0)) then
     ip = 0
  else
     ip = binary_search_L(varray, v) - 1
  endif
  p = (v - varray(ip)) / (varray(ip+1) - varray(ip))

  end subroutine poloidal_coordinates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function aux_eval(this, it, t, v) result(rz)
  class(torosurf), intent(in) :: this
  integer,         intent(in) :: it
  real(real64),    intent(in) :: t, v
  real(real64)                :: rz(2)

  real(real64) :: rztmp(2,2), p
  integer :: ip


  rz = 0.d0
  call this%poloidal_coordinates(it, t, v, ip, p)
  if (ip == -1) return

  rztmp = (1.d0-p) * this%rz(:,ip,it:it+1)  +  p * this%rz(:,ip+1,it:it+1)
  rz = (1.d0-t) * rztmp(:,1)  +  t * rztmp(:,2)

  end function aux_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(rzphi)
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: x(this%ndim)
  real(real64)                :: rzphi(this%mdim)

  real(real64) :: t
  integer :: it


  rzphi(3) = x(1)
  call this%toroidal_coordinates(x(1), it, t)
  if (it == -1) return

  rzphi(1:2) = this%aux_eval(it, t, x(2))

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: x(this%ndim)
  real(real64)                :: jac(this%mdim, this%ndim)


  ! @todo: implement torosurf%jac
  write (6, *) "ERROR: torosurf%jac not implemented yet!"
  stop

  end function jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp(this, u) result(x)
  !
  ! interpolate (r, z, phi) at u
  !
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: u(2)
  real(real64)                :: x(3)

  real(real64) :: rz(2,2)
  integer :: i, j


  i = int(u(1))
  j = int(u(2))
  rz = this%rz(:,j:j+1,i) + (this%rz(:,j:j+1,i+1) - this%rz(:,j:j+1,i)) * (u(1) - i)
  x(1:2) = rz(:,1) + (rz(:,2) - rz(:,1)) * (u(2) - j)
  x(3) = this%phi(i) + (this%phi(i+1) - this%phi(i)) * (u(1) - i)

  end function interp
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vphi(this, u)
  !
  ! convert u to (v, phi [deg])
  !
  use moose_math, only: pi
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: u(2)
  real(real64)                :: vphi(2)

  real(real64) :: v(2)
  integer :: i, j


  i = int(u(1))
  j = int(u(2))
  v = this%v(j:j+1,i) + (this%v(j:j+1,i+1) - this%v(j:j+1,i)) * (u(1) - i)
  vphi(1) = v(1) + (v(2) - v(1)) * (u(2) - j)
  vphi(2) = (this%phi(i) + (this%phi(i+1) - this%phi(i)) * (u(1) - i)) * 180 / pi

  end function vphi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function normal_vector(this, u) result(v)
  use moose_math, only: cross_product
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: u(2)
  real(real64)                :: v(3)

  real(real64) :: b(2), c(2), d(2), e1(3), e2(3), r, t, p
  integer :: it, ip


  it = int(u(1));   t = u(1) - it
  ip = int(u(2));   p = u(2) - ip
  b(1) = - this%B(2,ip,it)
  b(2) =   this%B(1,ip,it)
  c    =   this%C(:,ip,it)
  d(1) = - this%D(2,ip,it)
  d(2) =   this%D(1,ip,it)
  r = this%rz(1,ip,it) + p * b(1) + t * c(1) + p * t * d(1)

  ! e1 = dx / du
  e1(1:2) = c + p * d
  e1(3) = r * (this%phi(it+1) - this%phi(it))

  ! e2 = dx / dv
  e2(1:2) = b + t * d
  e2(3) = 0.d0

  ! n = e1 x e2
  v = cross_product(e2, e1)
  v = v / norm2(v)

  end function normal_vector
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function is_closed(this, phi)
  !
  ! check if surface slice is closed at *phi*
  !
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: phi
  logical                     :: is_closed

  real(real64) :: phi_mod, t, x0(2), x1(2)
  integer      :: i


  is_closed = .false.
  phi_mod = modulo(phi, this%dphiSym) + this%m0 * this%dphiSym
  call this%toroidal_coordinates(phi_mod, i, t)
  if (i < 0  .or.  i >= this%nu) return


  x0 = this%rz(:,0,i) + t * (this%rz(:,0,i+1) - this%rz(:,0,i))
  x1 = this%rz(:,this%nv,i) + t * (this%rz(:,this%nv,i+1) - this%rz(:,this%nv,i))
  is_closed = all(abs(x1 - x0) < 1.d2 * maxval(max(abs(x0), abs(x1))) * epsilon(1.d0))

  end function is_closed
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function area(this, i, j)
  !
  ! surface area of patch *(i, j)*
  !
  class(torosurf), intent(in) :: this
  integer,         intent(in) :: i, j
  real(real64)                :: area

  real(real64) :: r(2,2), z(2,2), r0, dphi, ds(2)


  r  = this%rz(1,j:j+1,i:i+1)
  z  = this%rz(2,j:j+1,i:i+1)
  r0 = sum(r) / 4
  ds = sqrt((r(2,:) - r(1,:))**2  +  (z(2,:) - z(1,:))**2)

  dphi = this%phi(i+1) - this%phi(i)
  area = dphi*r0 * sum(ds)/2

  end function area
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine lhshift(this, ds, ierr)
  !
  ! left-hand shift R-Z contours of torosurf
  !
  use moose_polygon2d
  class(torosurf), intent(inout) :: this
  real(real64),    intent(in   ) :: ds
  integer,         intent(  out) :: ierr

  type(polygon2d) :: P
  integer :: i


  P = polygon2d(this%nv)
  do i=0,this%nu
     call P%set_nodes(this%rz(:,:,i))
     call P%shift(ds, ierr)
     if (ierr /= 0) return
     if (P%segments() /= this%nv) then
        ierr = 2
        return
     endif

     this%rz(:,:,i) = transpose(P%nodes())
  enddo

  end subroutine lhshift
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function intersect(this, p1, p2, px, t, u)
  !
  ! calculate intersection between trajectory p1 -> p2 with non-zero toroidal
  ! increment and toroidal surface
  !
  use moose_algorithms
  use moose_math, only: pi2
  class(torosurf), intent(in)  :: this
  real(real64),    intent(in)  :: p1(this%mdim), p2(this%mdim)
  real(real64),    intent(out) :: px(this%mdim), t, u(this%ndim)
  logical                      :: intersect

  real(real64), parameter :: eps = 1.d-8

  real(real64) :: dp(2), dphi, mdphi, p1m(3), p2m(3), bbox(2,2), ulocal, vlocal, xi(2)
  integer :: m, m1, m2, dm, i, iA, iB, j


  ! set RZ-bounding box for p1 -> p2
  bbox(1,1) = min(p1(1), p2(1))
  bbox(2,1) = max(p1(1), p2(1))
  bbox(1,2) = min(p1(2), p2(2))
  bbox(2,2) = max(p1(2), p2(2))

  ! increments in RZ- and toroidal-direction
  dp    = p2(1:2) - p1(1:2)
  dphi  = p2(3)   - p1(3)

  ! p1 -> p2 may go across boundary period
  p1m   = p1;   m1  = int(floor(p1(3) / this%dphiSym))
  p2m   = p2;   m2  = int(floor(p2(3) / this%dphiSym))
  dm    = 1;   if (dphi < 0.d0) dm = -1
  ! scan boundary periods (should only be 1 or 2)
  intersect = .false.
  boundary_periods: do m=m1,m2,dm
     ! 1. calculate local toroidal angles (with respect to boundary period)
     mdphi  = (m - this%m0) * this%dphiSym
     p1m(3) = p1(3) - mdphi
     p2m(3) = p2(3) - mdphi


     ! 2. find first and last index in toroidal direction for segments which may intersect
     ! -1 since lower index = 0 in phi
     if (dm > 0) then
        iA = binary_search_L(this%phi, p1m(3)) - 1
        iB = binary_search_R(this%phi, p2m(3)) - 1
     else
        iA = binary_search_R(this%phi, p1m(3)) - 1
        iB = binary_search_L(this%phi, p2m(3)) - 1
     endif

     ! p1 and p2 are both below lower boundary of surface -> no intersection possible
     if (iA == -1  .and.  iB == -1) cycle

     ! p1 and p2 are both above upper boundary of surface -> no intersection possible
     if (iA == this%nu  .and.  iB == this%nu) cycle

     ! trim indices (either of p1 and p2 may still be outside toroidal boundaries)
     iA = max(min(iA,this%nu-1),0)
     iB = max(min(iB,this%nu-1),0)


     ! 3. loop over all toroidal segments in range
     toroidal_scan: do i=iA,iB,dm
        ! check if toroidal width of surface segment is small
        ! this includes the case when both this%dphi(i) = dphi = 0
        if (this%dphi(i) <= eps * abs(dphi)) then
           call poloidal_scan_planar_elements(i, j, t, xi, px, intersect)
           ulocal = (xi(1) + 1.d0) / 2
           vlocal = (xi(2) + 1.d0) / 2
        else
           call default_poloidal_scan(i, j, t, ulocal, vlocal, px, intersect)
        endif

        if (intersect) then
           u(1) = i + ulocal
           u(2) = j + vlocal
           return
        endif
     enddo toroidal_scan
  enddo boundary_periods

  contains
  !.............................................................................
  subroutine poloidal_scan_planar_elements(i, j, t, xi, px, intersect)
  !
  ! check for intersection with planar-like elements at phi(i)
  !
  use moose_polygon2d, only: winding_number
  use moose_quad
  integer,      intent(in)  :: i
  integer,      intent(out) :: j
  real(real64), intent(out) :: t, xi(2), px(3)
  logical,      intent(out) :: intersect

  type(quad) :: Q
  integer :: jj


  intersect = .false.
  if (dphi == 0.d0) return

  t = (this%phi(i) - p1m(3)) / dphi
  px = p1 + (p2 - p1) * t
  scan: do jj=0,this%nv-1
     Q%x(:,1) = this%rz(:,jj,i)
     Q%x(:,2) = this%rz(:,jj+1,i)
     Q%x(:,3) = this%rz(:,jj+1,i+1)
     Q%x(:,4) = this%rz(:,jj,i+1)
     if (Q%winding_number(px(1:2)) /= 0) then
        j = jj
        xi = Q%inverse_transform(px(1:2))
        intersect = .true.
        exit
     endif
  enddo scan

  end subroutine poloidal_scan_planar_elements
  !.............................................................................
  !pure subroutine default_poloidal_scan(i, j, u, v, px, istat)
  subroutine default_poloidal_scan(i, j, t, u, v, px, intersect)
  !
  ! calculate intersection with toroidal segment i
  !
  use moose_algorithms, only: solve_quadratic_equation
  integer,      intent(in)  :: i
  integer,      intent(out) :: j
  real(real64), intent(out) :: t, u, v, px(3)
  logical,      intent(out) :: intersect

  real(real64) :: alpha0, alpha1, r0(2), r1(2), w0(2), w1(2)
  real(real64) :: c1, c2, c3, tmin, tmax, tA, tB, tij(2), uij, vij, pxij(3)
  integer :: jj, k, n
  associate (A => this%rz, B => this%B, C => this%C, D => this%D)


  ! initialize t-domain of interest
  tmin = 0.d0
  tmax = 1.d0
  t    = 1.d0 + epsilon(1.d0)

  ! mapping t <-> u (direction of trajectory (T) and toroidal direction of surface element)
  alpha0 = (p1m(3) - this%phi(i)) / this%dphi(i)
  alpha1 =                   dphi / this%dphi(i)

  ! map toroidal boundaries of surface elements to t-domain
  if (dphi /= 0) then
     tA = -alpha0 / alpha1
     tB = (1.d0 - alpha0) / alpha1
     tmin = max(min(tA,tB), 0.d0)
     tmax = min(max(tA,tB), 1.d0)
  endif


  intersect = .false.
  ! poloidal scan within toroidal segment i
  scan_segments: do jj=0,this%nv-1
     ! check bounding box for segment (i,j)
     if (bbox(1,1) > this%bbox(2,1,jj,i)) cycle
     if (bbox(2,1) < this%bbox(1,1,jj,i)) cycle
     if (bbox(1,2) > this%bbox(2,2,jj,i)) cycle
     if (bbox(2,2) < this%bbox(1,2,jj,i)) cycle


     ! calculate intersection with segment (i,j)
     r0 = B(:,jj,i)  +  alpha0 * D(:,jj,i)
     r1 =               alpha1 * D(:,jj,i)
     w0 = p1m(1:2)  -  A(:,jj,i)  -  alpha0 * C(:,jj,i)
     w1 =                     dp  -  alpha1 * C(:,jj,i)

     c1 = sum(r1 * w1)
     c2 = sum(r1 * w0  +  r0 * w1)
     c3 = sum(r0 * w0)
     call solve_quadratic_equation(c1, c2, c3, tmin, tmax, tij, n)
     ! no intersection found -> continue with next segment
     check_roots: do k=1,n
        ! continue with next root/segment if intersection is further away than previous one
        if (tij(k) >= t) cycle

        ! evaluate this intersection point ................................
        pxij = p1     + tij(k) * (p2 - p1)
        uij  = alpha0 + tij(k) * alpha1

        ! calculate poloidal coordinate v
        r0 = B(:,jj,i)             + uij * D(:,jj,i)
        w0 = pxij(1:2) - A(:,jj,i) - uij * C(:,jj,i)
        if (abs(r0(1)) > abs(r0(2))) then
           vij = - w0(2) / r0(1)
        else
           vij =   w0(1) / r0(2)
        endif

        ! continue with next segment if intersection is outside of poloidal range
        if (vij < 0.d0  .or.  vij > 1.d0) cycle
        !..................................................................


        ! this is a good intersection point
        j     = jj
        u     = uij
        v     = vij
        t     = tij(k)
        px    = pxij
        intersect = .true.

        ! continue with remaining segments in case of multiple intersection points
        exit check_roots
     enddo check_roots
  enddo scan_segments

  end associate
  end subroutine default_poloidal_scan
  !.............................................................................
  end function intersect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rzslice_intersect(this, p1, p2, px, t, u) result(intersect)
  !
  ! intersection check for p1(3) = p2(3)
  !
  use moose_math, only: sign_test
  class(torosurf), intent(in   ) :: this
  real(real64),    intent(in   ) :: p1(this%mdim), p2(this%mdim)
  real(real64),    intent(  out) :: px(this%mdim), t, u(this%ndim)
  logical                        :: intersect

  real(real64) :: phi_mod, d1, d2, tj, ui, uj, vp(2), vn(2), x1(2), x2(2), xj(2), rmin, rmax, zmin, zmax, sumvp2
  integer      :: i, ierr, j


  intersect = .false.
  phi_mod = modulo(p1(3), this%dphiSym) + this%m0 * this%dphiSym
  call this%toroidal_coordinates(phi_mod, i, ui)
  if (i < 0  .or.  i >= this%nu) return
  px(3) = p1(3)
  u(1) = i + ui
  t = 1.d0


  ! bounding box, tangent and normal vectors for line segment L: p1 -> p2
  rmin = min(p1(1), p2(1))
  rmax = max(p1(1), p2(1))
  zmin = min(p1(2), p2(2))
  zmax = max(p1(2), p2(2))
  vp = p2(1:2) - p1(1:2)
  vn = [-vp(2), vp(1)]
  sumvp2 = sum(vp**2)
  if (sumvp2 == 0.d0) return
  vp = vp / sumvp2


  do j=1,this%nv
     x1 = this%rz(:,j-1,i) + ui * (this%rz(:,j-1,i+1) - this%rz(:,j-1,i))
     x2 = this%rz(:,j,i) + ui * (this%rz(:,j,i+1) - this%rz(:,j,i))
     if (max(x1(1), x2(1)) < rmin  .or.  min(x1(1), x2(1)) > rmax  .or. &
         max(x1(2), x2(2)) < zmin  .or.  min(x1(2), x2(2)) > zmax) cycle

     d1 = sum((x1 - p1(1:2)) * vn)
     d2 = sum((x2 - p1(1:2)) * vn)

     ! x1 and x2 are on the same side of L
     if (sign_test(d1, d2) == 1) cycle

     ! x1 is located on L
     if (d1 == 0.d0) then
        uj = 0.d0

     ! x2 is located on L
     elseif (d2 == 0.d0) then
        uj = 1.d0

     ! intersection between x1 and x2
     else
        uj = d1 / (d1 - d2)
     endif
     xj = x1 + uj * (x2 - x1)

     ! intersection between p1 and p2
     tj = sum((xj - p1(1:2)) * vp)
     if (tj >= 0.d0  .and.  tj <= t) then
        intersect = .true.
        t = tj
        u(2) = j - 1 + uj
        px(1:2) = xj
     endif
  enddo

  end function rzslice_intersect
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function includes(this, phi, allow_symmetry)
  !
  ! return true if toroidal domain of torosurf includes phi [rad]
  !
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: phi
  logical,         intent(in) :: allow_symmetry
  logical                     :: includes

  real(real64) :: phi_
  integer :: m


  if (allow_symmetry) then
     m = this%m0 - int(floor(phi / this%dphiSym))
     phi_ = phi + m * this%dphiSym
  else
     phi_ = phi
  endif


  includes = .true.
  if (phi_ < this%phi(0)  .or.  phi_ > this%phi(this%nu)) includes = .false.

  end function includes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function winding_number_phi(this, p) result(wn)
  !
  ! compute winding number for p = (r, z, phi)
  !
  use moose_polygon2d, only: winding_number
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: p(3)
  integer                     :: wn


  wn = winding_number(this%rzslice(p(3)), p(1:2), .false.)

  end function winding_number_phi
  !-----------------------------------------------------------------------------
  pure function winding_number_iphi(this, p, iphi) result(wn)
  !
  ! compute winding number for p = (r, z) at iphi-th toroidal index
  !
  use moose_polygon2d, only: winding_number
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: p(2)
  integer,         intent(in) :: iphi
  integer                     :: wn


  wn = winding_number(this%rz(:, :, iphi), p, .false.)

  end function winding_number_iphi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function vcurve(this, phi)
  !
  ! generate curve along v-coordinate at fixed phi
  !
  use moose_curve, only: curve
  use moose_bspline_curve, only: bspline_polygon
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: phi
  class(curve), allocatable   :: vcurve


  vcurve = bspline_polygon(this%polygon2d(phi))

  end function vcurve
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function rzslice(this, phi) result(x)
  !
  ! generate slice at toroidal position phi [rad]
  !
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: phi
  real(real64)                :: x(2, 0:this%nv)

  real(real64) :: phi_mod, t
  integer      :: i, j


  phi_mod = modulo(phi, this%dphiSym) + this%m0 * this%dphiSym
  call this%toroidal_coordinates(phi_mod, i, t)
  if (i < 0  .or.  i >= this%nu) then
     x = 0.d0
     return
  endif


  do j=0,this%nv
     x(:,j) = this%rz(:,j,i) + t * (this%rz(:,j,i+1) - this%rz(:,j,i))
  enddo

  end function rzslice
  !-----------------------------------------------------------------------------
  function polygon2d_slice(this, phi) result(P)
  use moose_error
  use moose_utils, only: str
  use moose_polygon2d
  class(torosurf), intent(in) :: this
  real(real64),    intent(in) :: phi
  type(polygon2d)             :: P


  if (.not.this%includes(phi, .true.)) call ERROR("profile is not defined at phi = "//str(phi))
  P = polygon2d(this%rzslice(phi))

  end function polygon2d_slice
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_tpzmesh3d(this) result(grid)
  use moose_math,  only: pi
  use moose_grids, only: tpzmesh3d
  class(torosurf), intent(in) :: this
  type(tpzmesh3d)             :: grid

  real(real64), allocatable :: rz(:,:,:)


  allocate (rz(2, 0:this%nu, 0:this%nv))
  rz(1,:,:) = transpose(this%rz(1,:,:))
  rz(2,:,:) = transpose(this%rz(2,:,:))

  grid = tpzmesh3d(transpose(this%v), this%phi * 180.d0/pi, rz, &
     "Toroidal Angle [deg]", this%vlabel(), "R", "Z")
  deallocate (rz)

  end function make_tpzmesh3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_r3grid(this) result(G)
  use moose_units
  use moose_grids, only: r3grid, cylindrical_r3grid
  class(torosurf), intent(in) :: this
  type(r3grid)                :: G


  G = cylindrical_r3grid(this%tpzmesh3d(), METER, DEGREE, 1, 2, 3, .true.)

  end function make_r3grid
  !-----------------------------------------------------------------------------

end module moose_torosurf
