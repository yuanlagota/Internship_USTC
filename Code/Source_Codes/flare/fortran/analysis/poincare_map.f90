#include <txtio.h>

module flare_poincare_map
  use iso_fortran_env
  use moose_txtio
  use moose_rlist
  use flare_fieldline
  implicit none
  private



  integer, parameter :: &
     BSPLINE_FIT_LINEAR    = 0, &
     BSPLINE_FIT_NONLINEAR = 1



  ! Poincare map for reference point p0 ........................................
  type, extends(txtio), public :: poincare_map
     ! initial point (r[m], z[m], phi[rad])
     real(real64) :: p0(3)
     integer      :: idir

     ! toroidal location of Poincare section [rad]
     real(real64) :: phiX
     ! toroidal symmetry
     integer :: nsym

     ! return points in Poincare section
     type(rlist) :: points
     integer :: istat

     contains
     procedure :: broadcast, free, send, write_formatted, savenc, writenc

     procedure :: bspline_multifit, bspline_autofit, harmonic_autofit
  end type poincare_map


  interface poincare_map
     procedure :: new
     procedure :: load
  end interface poincare_map
  ! poincare_map ...............................................................



  public :: &
     loadnc_poincare_map, readnc_poincare_map, &
     poincare_maps, &
     recv_poincare_map


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(p0, idir, phiX, nsym, npoints) result(this)
  real(real64), intent(in) :: p0(3), phiX
  integer,      intent(in) :: idir, nsym, npoints
  type(poincare_map)       :: this


  call init_txtio(this, "poincare_map")
  this%p0     = p0
  this%idir   = idir
  this%phiX   = phiX
  this%nsym   = nsym
  this%points = rlist(4, chunk_size=npoints)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename) result(this)
  use moose_dict
  use moose_math, only: pi
  character(len=*), intent(in) :: filename
  type(poincare_map)           :: this

  type(dict) :: metadata
  integer    :: iu


  open  (newunit=iu, file=filename)
  metadata = read_metadata(iu, "poincare_map")
  this%p0   = metadata%getreal_rank1("P0", 3);   this%p0(3) = this%p0(3) / 180.d0 * pi
  this%idir = metadata%getint("DIRECTION")
  this%phiX = metadata%getreal("PHIX") / 180.d0 * pi
  this%nsym = metadata%getint("NSYMMETRY")
  close (iu)

  this%points = rlist(filename)

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_poincare_map(filename) result(this)
  !
  ! load poincare_map from netcdf file
  !
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(poincare_map)           :: this

  type(netcdf_dataset) :: nc


  nc = netcdf_open(filename)
  this = readnc_poincare_map(nc)
  call nc%close()

  end function loadnc_poincare_map
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function readnc_poincare_map(nc) result(this)
  !
  ! read poincare_map from netcdf dataset
  !
  use moose_math, only: pi
  use moose_netcdf
  type(netcdf_dataset), intent(in) :: nc
  type(poincare_map)               :: this

  real(real64), pointer :: values(:,:)
  real(real64) :: p0(3), phiX
  integer :: idir, npoints, nsym


  npoints = nc%dim("npoints")
  call nc%get_att("direction", idir)
  call nc%get_att("symmetry", nsym)
  call nc%get_att("r0", p0(1))
  call nc%get_att("z0", p0(2))
  call nc%get_att("phi0", p0(3));   p0(3) = p0(3) / 180 * pi
  call nc%get_att("phiX", phiX)

  call init_txtio(this, "poincare_map")
  this%p0     = p0
  this%idir   = idir
  this%phiX   = phiX
  this%nsym   = nsym
  this%points = readnc_rlist(nc, "points")

  end function readnc_poincare_map
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function poincare_maps(p0, idir, phiX, nsym, npoints, nsections, bounded, fdriver) result(P)
  !
  ! generate Poincare map for reference point p0 by tracing field line in
  ! direction idir to R-Z plane at phiX with toroidal symmetry nsym with npoints
  ! return points
  !
  use moose_error,  only: DOMAIN_ERROR
  use moose_math,   only: pi2, deg3
  use flare_model,  only: bfield
  use flare_fieldline, only: fieldline_driver => fdriver
  real(real64),                    intent(in) :: p0(3), phiX
  integer,                         intent(in) :: idir, nsym, npoints, nsections
  logical,                         intent(in), optional :: bounded
  type(fieldline_driver), pointer, intent(in), optional :: fdriver
  type(poincare_map)                          :: P(0:nsections-1)

  type(fieldline_driver), pointer :: F
  real(real64) :: dphi, dphiX, phi, phiXX, psiN, theta, y(3)
  integer :: i, j


  ! initialize field line driver (or link to user provided one)
  if (present(fdriver)) then
     F => fdriver
     call F%reset()
  else
     allocate (F, source=fieldline_driver(bounded))
  endif


  ! initialize parameters for Poincare map
  dphi = idir * pi2 / nsym
  do i=0,nsections-1
     phi  = phiX + dphi / nsections * i
     P(i) = poincare_map(p0, idir, phi, nsym, npoints)
  enddo


  ! trace field line from p0 to phiX (mod 2 pi / nsym)
  y     = p0
  dphiX = phiX - y(3)
  if (abs(dphiX) > 0.d0) then
     phiXX = y(3) + modulo(dphiX, dphi)
     P(:)%istat = F%evolve3(y, phiXX)
     if (P(0)%istat == INTERSECT_BOUNDARY  .or.  P(0)%istat > 0) then
        call free_fdriver();   return
     endif
  endif
  phiXX = y(3)


  ! generate Poincare maps
  trace: do i=0,npoints-1
     do j=1,nsections
        phi = phiXX  +  i * dphi  +  j * dphi / nsections
        P(:)%istat = F%evolve3(y, phi)
        if (P(0)%istat == INTERSECT_BOUNDARY  .or.  P(0)%istat > 0) then
           exit trace
        endif

        theta = bfield%equi%poloidal_angle(y);   if (theta < 0.d0) theta = theta + pi2
        psiN = bfield%equi%psiN(y)
        call P(mod(j, nsections))%points%append([y(1), y(2), theta/pi2*360.d0, psiN])
     enddo
  enddo trace
  call free_fdriver()

  contains
  !.............................................................................
  subroutine free_fdriver()


  if (.not.present(fdriver)) then
     call F%free()
     deallocate (F)
  endif

  end subroutine free_fdriver
  !.............................................................................
  end function poincare_maps
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function recv_poincare_map(from) result(this)
  use moose_mpi
  integer, intent(in) :: from
  type(poincare_map)  :: this


  call proc(from)%recv(this%p0)
  call proc(from)%recv(this%idir)
  call proc(from)%recv(this%phiX)
  call proc(from)%recv(this%nsym)
  this%points = recv_rlist(from)

  end function recv_poincare_map
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(poincare_map), intent(inout) :: this


  call this%txtio_broadcast()
  call proc(0)%broadcast(this%p0)
  call proc(0)%broadcast(this%idir)
  call proc(0)%broadcast(this%phiX)
  call proc(0)%broadcast(this%nsym)
  call proc(0)%broadcast(this%istat)
  call this%points%broadcast()

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(poincare_map), intent(inout) :: this


  call this%points%free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine send(this, to)
  use moose_mpi
  class(poincare_map), intent(inout) :: this
  integer,             intent(in   ) :: to


  call proc(to)%send(this%p0)
  call proc(to)%send(this%idir)
  call proc(to)%send(this%phiX)
  call proc(to)%send(this%nsym)
  call this%points%send(to)

  end subroutine send
  !-----------------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  use moose_math, only: pi
  class(poincare_map), intent(in   ) :: this
  integer,             intent(in   ) :: unit, vlist(:)
  character(len=*),    intent(in   ) :: iotype
  integer,             intent(  out) :: iostat
  character(len=*),    intent(inout) :: iomsg

  real(real64) :: p0(3)
  integer :: n


  n = this%points%nelements()
  p0 = this%p0;   p0(3) = p0(3) / pi * 180.d0
  WRITETXT(metadata_fmt("P0",        "3(x,g0.14)"), p0)
  WRITETXT(metadata_fmt("DIRECTION", "i0"        ), this%idir)
  WRITETXT(metadata_fmt("PHIX",      "g0.14"     ), this%phiX / pi * 180.d0)
  WRITETXT(metadata_fmt("NSYMMETRY", "i0"        ), this%nsym)
  WRITETXT(metadata_fmt("POINTS",    "i0"        ), n)
  call this%points%write_formatted(unit, iotype, vlist, iostat, iomsg)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  class(poincare_map), intent(in) :: this
  character(len=*),    intent(in) :: filename

  type(netcdf_dataset) :: nc


  nc = netcdf_create(filename)
  call this%writenc(nc)
  call nc%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine writenc(this, nc)
  use moose_math, only: pi
  use moose_netcdf
  class(poincare_map),  intent(in) :: this
  type(netcdf_dataset), intent(in) :: nc

  integer :: ndim, npoints


  call nc%def_dim("dim_0004", 4, ndim)
  call nc%def_dim("npoints", this%points%nelements(), npoints)
  call nc%put_att("r0", this%p0(1))
  call nc%put_att("z0", this%p0(2))
  call nc%put_att("phi0", this%p0(3) / pi * 180)
  call nc%put_att("direction", this%idir)
  call nc%put_att("phiX", this%phiX / pi * 180)
  call nc%put_att("symmetry", this%nsym)
  call nc%def_var("points", NF90_DOUBLE, [ndim, npoints])
  call nc%enddef()

  call nc%put_var("points", this%points%values())

  end subroutine writenc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline_multifit(this, nctrl, k, theta0, knot_balancing, method, epsabs, lambda1, lambda2) result(B)
  !
  ! Fit (order *k*) B-Spline curve with *nctrl* control points to Poincare map
  !
  use moose_error
  use moose_utils,    only: str, user_option
  use moose_math,     only: pi
  use moose_geometry, only: multifit => bspline_multifit
  use moose_bfit
  class(poincare_map), intent(in) :: this
  integer,             intent(in) :: nctrl, k
  real(real64),        intent(in), optional :: theta0, epsabs, lambda1, lambda2
  logical,             intent(in), optional :: knot_balancing
  integer,             intent(in), optional :: method
  type(bfit)                         :: B

  real(real64), pointer :: rz(:,:)
  real(real64), allocatable :: theta(:)
  real(real64) :: chisq, theta0_


  allocate (theta(this%points%nelements()))
  theta = this%points%column(3);   theta0_ = 0.d0
  rz => this%points%columns(1,2)
  if (present(theta0)) then
     theta0_ = theta0
     where (theta < theta0) theta = theta + 360.d0
     where (theta > theta0 + 360.d0) theta = theta - 360.d0
  endif

  select case(user_option(BSPLINE_FIT_LINEAR, method))
  case (BSPLINE_FIT_LINEAR)
     B%bspline_curve = multifit(theta, transpose(rz), nctrl, theta0_, theta0_+360.d0, .true., k, knot_balancing, chisq=chisq)
     allocate (B%e(0:0), source=chisq)

  case (BSPLINE_FIT_NONLINEAR)
     B = bfit(rz, nctrl, k, user_option(1.d-5,epsabs), lambda1, lambda2)

  case default
     call ERROR("invalid method = "//str(method), "bspline_multifit")
  end select
  deallocate (theta)

  end function bspline_multifit
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline_autofit(this, k, eps, theta0, knot_balancing, method) result(B)
  !
  ! Fit (order *k*) B-Spline curve to Poincare map. The number of control points
  ! is automaticaly increased until the accuracy *eps* is reached.
  !
  use moose_error
  use moose_bfit
  class(poincare_map), intent(in   ) :: this
  integer,             intent(in   ) :: k
  real(real64),        intent(in   ) :: eps
  real(real64),        intent(in   ), optional :: theta0
  logical,             intent(in   ), optional :: knot_balancing
  integer,             intent(in   ), optional :: method
  type(bfit)                         :: B

  real(real64) :: chisq_dof
  integer :: n, nmax


  nmax = log(1.d0 * this%points%nelements() / k) / log(2.d0)
  if (nmax < 3) call ERROR("too few points in Poincare map")
  do n=3,nmax
     B = this%bspline_multifit(2**n, k, theta0, knot_balancing, method, eps)
     chisq_dof = B%e(0) / 2**n
     if (chisq_dof < eps) exit
  enddo

  end function bspline_autofit
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function harmonic_autofit(this, updown_symmetry, eps, npoints, report, debug, max_order, max_iterations, &
     damping, footpoint_accuracy, chisq_dof) result(C)
  use moose_geometry, only: fourier_curve, fourier_autofit
  class(poincare_map), intent(in   ) :: this
  logical,             intent(in   ) :: updown_symmetry
  real(real64),        intent(in   ) :: eps
  logical,             intent(in   ), optional :: report, debug
  integer,             intent(in   ), optional :: npoints, max_order, max_iterations
  real(real64),        intent(in   ), optional :: damping, footpoint_accuracy
  real(real64),        intent(  out), optional :: chisq_dof
  type(fourier_curve)                :: C

  real(real64), pointer :: x(:,:)
  integer :: n


  x => this%points%columns(1,2)
  n = size(x,2);   if (present(npoints)) n = npoints
  C = fourier_autofit(x(:,0:n-1), updown_symmetry, eps, report, debug, max_order, &
     max_iterations, damping, footpoint_accuracy, chisq_dof)

  end function harmonic_autofit
  !-----------------------------------------------------------------------------

end module flare_poincare_map
