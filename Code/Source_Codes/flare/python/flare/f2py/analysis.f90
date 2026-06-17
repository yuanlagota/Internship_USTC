module analysis
  use kinds
  implicit none


  ! workspace for tracing field lines
  real(real64), allocatable :: fieldline_x(:,:), fieldline_s(:)
  logical :: fieldline_bounded


  ! workspace for tracing radial paths
  real(real64), allocatable :: rpath2d_x(:,:), rpath2d_s(:)
  integer :: rpath2d_idir


  ! workspace for fluxsurf2d
  real(real64), allocatable :: fluxsurf2d_x(:,:), fluxsurf2d_arcl(:), fluxsurf2d_theta(:), &
     fluxsurf2d_alpha(:)
  real(real64) :: fluxsurf2d_params(4)


  ! workspace for polygons
  real(real64), allocatable :: polygon2d_x(:,:)


  contains
  !-----------------------------------------------------------------------------


include "_analysis.f90"


! bfield
  !-----------------------------------------------------------------------------
  function bfield_eval(n, x) result(B)
  use flare_model
  integer,      intent(in) :: n
  real(real64), intent(in) :: x(3,*)
  real(real64)             :: B(3,n)

  integer :: i


  ! TODO: parallelization
  do i=1,n
     B(:,i) = bfield%eval(x(:,i))
  enddo

  end function bfield_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bfield_perturbation_eval(x) result(B)
  use flare_model
  real(real64), intent(in) :: x(3)
  real(real64)             :: B(3)


  B = bfield%perturbation_eval(x)

  end function bfield_perturbation_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bfield_jac(x) result(J)
  use flare_model
  real(real64), intent(in) :: x(3)
  real(real64)             :: J(3,3)


  J = bfield%jac(x)

  end function bfield_jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function iquery_equilibrium(varname) result(i)
  use flare_model
  character(len=*), intent(in) :: varname
  integer                      :: i


  i = 0.d0
  call assert_model()
  select case(varname)
  case("Bt_sign")
     i = bfield%equi%Bt_sign

  case("Bp_sign")
     i = bfield%equi%Bp_sign
  end select

  end function iquery_equilibrium
  !-----------------------------------------------------------------------------


! boundary
  !-----------------------------------------------------------------------------
  subroutine boundary_firstwall_rzslice(phi)
  use moose_geometry, only: polygon2d
  use flare_boundary
  real(real64), intent(in) :: phi

  type(polygon2d) :: P


  if (allocated(polygon2d_x)) deallocate (polygon2d_x)
  P = firstwall_rzslice(phi)
  allocate (polygon2d_x, source=P%nodes())

  end subroutine boundary_firstwall_rzslice
  !-----------------------------------------------------------------------------


! fieldline
  !-----------------------------------------------------------------------------
  subroutine fieldline_trace(x0, idir, ds, nsteps, stop_at_boundary, coordinates, angular_units)
  use flare_model
  use flare_fieldline
  real(real64),     intent(in) :: x0(3), ds
  integer,          intent(in) :: idir, nsteps
  logical,          intent(in) :: stop_at_boundary
  character(len=*), intent(in) :: coordinates, angular_units

  type(fieldline) :: F
  integer :: n


  call assert_model()
  if (allocated(fieldline_x)) deallocate(fieldline_x, fieldline_s)
  F = fieldline(x0, idir, ds, nsteps, stop_at_boundary, coordinates, angular_units)


  n = F%trace%nelements()
  allocate (fieldline_x(3, 0:n-1), source=0.d0)
  allocate (fieldline_s(n), source=0.d0)
  fieldline_x = F%trace%columns(1, 3)
  fieldline_s = F%trace%column(4)
  fieldline_bounded = F%at_boundary

  end subroutine fieldline_trace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine trace_dphi(x0, dphi, istat, x1)
  use flare_model
  use flare_fieldline
  real(real64), intent(in   ) :: x0(3), dphi
  integer,      intent(  out) :: istat
  real(real64), intent(  out) :: x1(3)

  type(fdriver) :: F
  real(real64) :: phi1


  call assert_model()
  F = fdriver()

  x1 = x0
  phi1 = x1(3) + dphi
  istat = F%evolve3(x1, phi1)

  call F%free()

  end subroutine trace_dphi
  !-----------------------------------------------------------------------------


! rpath2d
  !-----------------------------------------------------------------------------
  subroutine export_rpath2d(values, idir)
  real(real64), intent(in) :: values(:,:)
  integer,      intent(in) :: idir

  integer :: n


  if (allocated(rpath2d_x)) deallocate(rpath2d_x, rpath2d_s)
  n = size(values, 2)
  allocate (rpath2d_x(2, 0:n-1), source=values(1:2,:))
  allocate (rpath2d_s(n), source=values(3,:))
  rpath2d_idir = idir

  end subroutine export_rpath2d
  !-----------------------------------------------------------------------------
  subroutine make_rpath2d_trace(x0, param, t1, bounded)
  use flare_model, only: assert_equi2d, boundary2d
  use flare_rpath2d
  real(real64),     intent(in) :: x0(2), t1
  character(len=*), intent(in) :: param
  logical,          intent(in) :: bounded

  type(rpath2d_trace) :: trace


  call assert_equi2d("rpath2d_trace")
  if (bounded) then
     trace = rpath2d_trace(x0, RPATH2D_PARAM(param), t1, boundary=boundary2d)
  else
     trace = rpath2d_trace(x0, RPATH2D_PARAM(param), t1)
  endif
  call export_rpath2d(trace%values(), trace%idir)

  end subroutine make_rpath2d_trace
  !-----------------------------------------------------------------------------
  subroutine make_rpath2d_tracex(ix, xdir, param, t1, bounded)
  use flare_model, only: assert_equi2d, boundary2d
  use flare_rpath2d
  integer,          intent(in) :: ix, xdir
  character(len=*), intent(in) :: param
  real(real64),     intent(in) :: t1
  logical,          intent(in) :: bounded

  type(rpath2d_trace) :: trace


  call assert_equi2d("rpath2d_traceX")
  if (bounded) then
     trace = rpath2d_traceX(ix, xdir, RPATH2D_PARAM(param), t1, boundary=boundary2d)
  else
     trace = rpath2d_traceX(ix, xdir, RPATH2D_PARAM(param), t1)
  endif
  call export_rpath2d(trace%values(), trace%idir)

  end subroutine make_rpath2d_tracex
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine make_poincare_maps(p0, idir, phiX, nsymmetry, npoints, nsections, bounded, x_out, phi_out, npoints_out)
  use moose_math, only: pi
  use flare_model
  use flare_poincare_map, only: poincare_map, poincare_maps
  integer,      intent(in   ) :: idir, nsymmetry, npoints, nsections
  real(real64), intent(in   ) :: p0(3), phiX
  logical,      intent(in   ) :: bounded
  real(real64), intent(  out) :: x_out(0:nsections-1,4,npoints), phi_out(0:nsections-1)
  integer,      intent(  out) :: npoints_out(0:nsections-1)

  type(poincare_map), allocatable :: P(:)
  integer :: i


  call assert_model()
  allocate (P(0:nsections-1))
  P = poincare_maps(p0, idir, phiX/180.d0*pi, nsymmetry, npoints, nsections, bounded)
  do i=0,nsections-1
     npoints_out(i) = P(i)%points%nelements()
     phi_out(i)     = P(i)%phiX / pi * 180.d0
     x_out(i,:,:)   = P(i)%points%values()
     call P(i)%free()
  enddo
  deallocate (P)

  end subroutine make_poincare_maps
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function iquery_equi2d(varname) result(i)
  use flare_model
  character(len=*), intent(in) :: varname
  integer                      :: i


  call assert_equi2d("query_equi2d")
  select case(varname)
  case("nx")
     i = equi2d%X%nelements()
  end select

  end function iquery_equi2d
  !-----------------------------------------------------------------------------
  function rquery_equi2d(varname) result(r)
  use flare_model
  character(len=*), intent(in) :: varname
  real(real64)                 :: r


  r = 0.d0
  call assert_equi2d("query_equi2d")
  select case(varname)
  case("Bt_axis")
     r = equi2d%Bt_axis

  case("R_axis")
     r = equi2d%r0(1)

  case("Z_axis")
     r = equi2d%r0(2)

  case("poloidal_flux")
     r = equi2d%delta_psi

  case("Psi_axis")
     r = equi2d%Psi_axis
  end select

  end function rquery_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function equi2d_psin(rz) result(psiN)
  use flare_model
  real(real64), intent(in) :: rz(2)
  real(real64)             :: psin


  call assert_equi2d("equi2d_psiN")
  psin = equi2d%psiN(rz)

  end function equi2d_psin
  !-----------------------------------------------------------------------------
  function equi2d_grad_psi(rz) result(grad_psi)
  use flare_model
  real(real64), intent(in) :: rz(2)
  real(real64)             :: grad_psi(2)


  call assert_equi2d("equi2d_grad_psi")
  grad_psi = equi2d%Psi%deriv(rz)

  end function equi2d_grad_psi
  !-----------------------------------------------------------------------------
  function equi2d_poloidal_angle(rz) result(theta)
  use flare_model
  real(real64), intent(in) :: rz(2)
  real(real64)             :: theta


  call assert_equi2d("equi2d_poloidal_angle")
  theta = equi2d%poloidal_angle(rz)

  end function equi2d_poloidal_angle
  !-----------------------------------------------------------------------------
  function equi2d_xpoint(ix) result(x)
  use flare_model
  integer, intent(in) :: ix
  real(real64)        :: x(2)


  call assert_equi2d("equi2d_xpoint")
  x = equi2d%xpoint(ix)

  end function equi2d_xpoint
  !-----------------------------------------------------------------------------
  subroutine equi2d_xpoint_hessian(ix, lambda1, lambda2, v1, v2)
  use flare_model
  integer,      intent(in   ) :: ix
  real(real64), intent(  out) :: lambda1, lambda2, v1(2), v2(2)

  real(real64) :: v(2,-1:1)


  call assert_equi2d("equi2d_xpoint_hessian")
  v = equi2d%xpoint_hessian(ix)

  lambda1 = v(1,0)
  lambda2 = v(2,0)
  v1 = v(:,-1)
  v2 = v(:, 1)

  end subroutine equi2d_xpoint_hessian
  !-----------------------------------------------------------------------------
  subroutine equi2d_xpoint_stability(ix, lambda1, lambda2, v1, v2)
  use flare_model
  integer,      intent(in   ) :: ix
  real(real64), intent(  out) :: lambda1, lambda2, v1(2), v2(2)

  real(real64) :: v(2,-1:1)


  call assert_equi2d("equi2d_xpoint_stability")
  v = equi2d%xpoint_stability(ix)

  lambda1 = v(1,0)
  lambda2 = v(2,0)
  v1 = v(:,-1)
  v2 = v(:, 1)

  end subroutine equi2d_xpoint_stability
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function equi2d_rzcoords(psin, theta, x0, use_x0) result(rz)
  use flare_model
  real(real64), intent(in) :: psin, theta, x0(2)
  logical,      intent(in) :: use_x0
  real(real64)             :: rz(2)

  integer :: ierr


  call assert_equi2d("equi2d_rzcoords")
  if (use_x0) then
     rz = equi2d%rzcoords(psin, theta, x0=x0, ierr=ierr)
  else
     rz = equi2d%rzcoords(psin, theta, ierr=ierr)
  endif

  end function equi2d_rzcoords
  !-----------------------------------------------------------------------------
  function equi2d_rzcoordsx(psin) result(rz)
  use flare_model
  real(real64), intent(in) :: psin
  real(real64)             :: rz(2)


  call assert_equi2d("equi2d_rzcoordsX")
  rz = equi2d%rzcoordsX(psin)

  end function equi2d_rzcoordsx
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function equi3d_psin(r) result(psiN)
  use flare_model, only: bfield
  real(real64), intent(in) :: r(3)
  real(real64)             :: psin


  psin = bfield%equi%psiN(r)

  end function equi3d_psin
  !-----------------------------------------------------------------------------


! fluxsurf2d
  !-----------------------------------------------------------------------------
  subroutine init_fluxsurf2d_workspace(F)
  use flare_fluxsurf2d, only: fluxsurf2d
  class(fluxsurf2d), intent(in) :: F


  ! cleanup workspace, if necessary
  if (allocated(fluxsurf2d_x)) then
     deallocate (fluxsurf2d_x, fluxsurf2d_arcl, fluxsurf2d_theta, fluxsurf2d_alpha)
  endif

  ! copy data to workspace
  allocate (fluxsurf2d_x,     source=transpose(F%u))
  allocate (fluxsurf2d_arcl,  source=F%arcl)
  allocate (fluxsurf2d_theta, source=F%theta)
  allocate (fluxsurf2d_alpha, source=F%alpha)
  fluxsurf2d_params = [F%q, F%area, F%Vprime, F%current]

  end subroutine init_fluxsurf2d_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_fluxsurf2d(x, idir, param)
  use flare_fluxsurf2d, only: fluxsurf2d
  use flare_model,      only: assert_equi2d, boundary2d
  real(real64),     intent(in) :: x(2)
  integer,          intent(in) :: idir
  character(len=*), intent(in) :: param

  type(fluxsurf2d) :: F


  call assert_equi2d("construct_fluxsurf2d")
  F = fluxsurf2d(x, idir, param, boundary2d)
  call init_fluxsurf2d_workspace(F)

  end subroutine construct_fluxsurf2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_last_closed_fluxsurf2d(idir)
  use flare_fluxsurf2d, only: fluxsurf2d, last_closed_fluxsurf2d
  integer, intent(in) :: idir

  type(fluxsurf2d) :: F


  F = last_closed_fluxsurf2d(idir)
  call init_fluxsurf2d_workspace(F)

  end subroutine construct_last_closed_fluxsurf2d
  !-----------------------------------------------------------------------------

end module analysis
