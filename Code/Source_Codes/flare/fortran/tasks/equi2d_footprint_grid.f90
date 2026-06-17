subroutine flare_task_equi2d_footprint_grid(strike_point, phic, nsym, nu, nv, vstart, vend, xpoint, fmt, output, label)
  !
  ! Construct grid for magnetic footprint analysis from toroidally symmetric equilibrium.
  !
  ! **Parameters:**
  !
  ! :strike_point:  Strike point tag:
  !
  !                 :fwd, forward:            in direction of poloidal field.
  !
  !                 :bwd, backward:           against direction of poloidal field.
  !
  !                 :cw, clockwise:           in clockwise direction.
  !
  !                 :ccw, counter-clockwise:  in counter-clockwise direction.
  !
  ! :phic:          Center of toroidal domain [deg].
  !
  ! :nsym:          Symmetry for toroidal domain [phic - 180/nsym, phic + 180/nsym].
  !
  ! :nu:            Number of nodes in toroidal direction.
  !
  ! :nv:            Number of nodes along target.
  !
  ! :vstart:        Begin of target domain [cm] relative to strike point.
  !
  ! :vend:          End of target domain [cm] relative to strike point.
  !
  ! :xpoint:        Select X-point.
  !
  ! :fmt:           Either ``r3grid`` (for :func:`fieldline_connection <fieldline_connection>` task) or ``torosurf`` (for EMC3 post-processing).
  !
  ! :output:        Filename for grid output.
  !
  ! :label:         User defined label for ``torosurf`` file.
  !
  use iso_fortran_env
  use moose_error
  use moose_contours, only: xcontour, XCONTOUR_BRANCH, XCONTOUR_POSITIVE_ORIENTATION
  use moose_geometry, only: shifted_polygon2d, bspline_curve, bspline_polygon, torosurf
  use moose_grids,    only: cgrid, rmesh3d, r3grid, cylindrical_r3grid
  use moose_math,     only: linspace, pi, CYLINDRICAL_COORDINATES
  use moose_units
  use moose_utils,    only: str
  use flare_control
  use flare_model,    only: boundary2d, equi2d, assert_equi2d
  use flare_fluxsurf2d
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: strike_point, fmt, output, label
  real(real64),     intent(in) :: phic, vstart, vend
  integer,          intent(in) :: nsym, nu, nv, xpoint

  real(real64), parameter :: offset = 1.d-6

  type(xcontour)      :: separatrix
  type(bspline_curve) :: C
  type(cgrid)         :: V
  type(rmesh3d)       :: M
  type(r3grid)        :: G
  type(torosurf)      :: S
  real(real64) :: sp(2), u(2), ub, t(2), n(2), ds, bounds(2)
  integer      :: i, ib, idir, k, vdir


  if (rank > 0) return
  ! greeting
  call begin_task()
  if (report) then
     print *, "Constructing grid for magnetic footprint analysis from toroidally symmetric equilibrium."
     print *
  endif
  call assert_equi2d("equi2d_footprint_grid")
  idir = -make_idir(strike_point, -equi2d%Bp_sign)
  if (report) then
     print 1001, xpoint
     print 1002, vstart, vend
  endif
 1001 format(8x,"X-point: ",i0/)
 1002 format(3x,"- Target domain: ",f0.3," cm -> ",f0.3," cm from strike point"/)


  ! 1. construct separatrix strike point on boundary, and set up geometry parameters
  separatrix = separatrix2d(xpoint)
  k  = XCONTOUR_BRANCH(idir, XCONTOUR_POSITIVE_ORIENTATION)
  if (idir == 1) call separatrix%branch(k)%L%reverse()
  sp = separatrix%branch(k)%point(-1)
  u  = separatrix%branch(k)%point(-2) - sp ! vector pointing upstream (along separatrix leg)
  ib = -separatrix%iconnect(k)
  if (ib <= 0) then
     print *, "ib = ", ib
     call separatrix%savetxt("ERROR_separatrix")
     print 9001, k
     stop
  endif
 9001 format("ERROR in equi2d_footprint_grid: branch ",i0," of separatrix does not connect to boundary")


  ! 2. boundary geometry at strike point (tangent, normal & orientation)
  ub   = separatrix%uconnect(k)
  t    = boundary2d%P(ib)%tangent(int(ub))
  n    = boundary2d%P(ib)%normal(int(ub))
  vdir = int(sign(1.d0, sum(t * equi2d%Psi%deriv(sp) / equi2d%delta_Psi)))


  ! 3. shift boundary contour in upstream direction & sample points along target
  ds   = int(sign(1.d0, sum(n*u))) * offset
  C    = bspline_polygon(shifted_polygon2d(boundary2d%P(ib), ds), origin=ub)
  bounds(1) = vdir * vstart / 100.d0
  bounds(2) = vdir * vend   / 100.d0
  if (any(bounds < C%a)  .or.  any(bounds > C%b)) then
     print 9002
     stop
  endif
 9002 format("ERROR in equi2d_footprint_grid: vstart or vend out of bounds!")

  V   = cgrid(nv, 2, "Distance along target [cm]", ["R [m]", "Z [m]"])
  V%t = linspace(vstart, vend, nv)
  do i=0,nv-1
     V%x(:,i) = C%eval(vdir * V%t(i) / 100.d0)
  enddo


  ! 4. construct 3D grid
  ! ... from cgrid2d in R-Z plane at phic
  if (nu == 0) then
     if (report) print 4001, phic
     G = cylindrical_r3grid(V, METER, DEGREE, 1, 2, phic, .true.)


  ! ... from product of linspace and cgrid2d
  elseif (nu > 0) then
     if (report) print 4002, nu, phic - 180.d0/nsym, phic + 180.d0/nsym
     M = rmesh3d(linspace(phic-180.d0/nsym, phic+180.d0/nsym, nu), V, "Toroidal Angle [deg]")
     G = cylindrical_r3grid(M, METER, DEGREE, 1, 2, 3, .true.)

  else
     print 9003
  endif
 4001 format(3x,"- Constructing grid at toroidal position ",f0.3," deg"/)
 4002 format(3x,"- Constructing mesh with ",i0," nodes in toroidal domain ",f0.3," -> ",f0.3," deg"/)
 9003 format("ERROR in equi2d_footprint_grid: nu >= 0 required")


  ! 5. save output
  select case(fmt)
  case("r3grid")
     call G%savetxt(output)

  case("torosurf")
     S = torosurf(max(nu-1,0), nv-1, nsym, vdef=.true.)

     ! set label
     if (label == "") then
        call S%metadata%set("name", "x"//str(xpoint)//" "//trim(strike_point)//" footprint grid")
     else
        call S%metadata%set("name", label)
     endif

     ! set phi array
     if (nu == 0) then
        S%phi = phic / 180.d0 * pi
     else
        S%phi = M%domain%u / 180.d0 * pi
     endif

     ! set R,Z-contour
     do i=0,nu-1
        S%rz(:,:,i) = V%x * 100.d0
        S%v(:,i) = V%t
     enddo
     call S%savetxt(output, iotype='legacy')

  case default
     call ERROR("undefined output format: "//fmt)
  end select
  call finalize_task()

end subroutine flare_task_equi2d_footprint_grid
