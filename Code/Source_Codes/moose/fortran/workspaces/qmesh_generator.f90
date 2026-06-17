module moose_qmesh_generator
  use moose_error
  use moose_grids,    only: rmesh, qmesh
  use moose_geometry, only: polygon2d, curve
  use moose_quantiles
  use moose_mfunc
  use moose_bspline2d
  use iso_fortran_env
  implicit none


  character(len=256) :: debug_dummyU = ""

  real(real64) :: bbox_margin = 0.25d0

  integer :: points_per_segment = 4


  ! auxiliary potential for quasi-orthogonal mesh construction -----------------
  type, extends(scalar_mfunc2d), public :: aux_potential
     real(real64), allocatable :: q(:), w(:), x(:,:)
     integer :: n

     contains
     procedure :: eval => aux_potential_eval
     procedure :: deriv => aux_potential_deriv
     procedure :: hessian => aux_potential_hessian
  end type aux_potential


  interface aux_potential
     procedure :: compute_aux_potential
  end interface aux_potential
  ! aux_potential --------------------------------------------------------------



  ! workspace for qmesh construction -------------------------------------------
  type, extends(qmesh) :: qwork
     contains
     procedure :: contour, contour_length, delta

     procedure :: insert_nodes

     procedure :: submesh
  end type qwork
  ! qwork ----------------------------------------------------------------------



  interface set_qmesh_contour
     procedure :: set_qmesh_contour_array
     procedure :: set_qmesh_contour_polygon2d
     procedure :: set_qmesh_contour_curve
  end interface set_qmesh_contour


  interface aux_qmesh_distance_contour_generator
     procedure :: aux_qmesh_distance_contour_generator_curve
     procedure :: aux_qmesh_distance_contour_generator_poly2d
  end interface aux_qmesh_distance_contour_generator


  contains
  !-----------------------------------------------------------------------------


! constructors (aux_potential):
  !-----------------------------------------------------------------------------
  function compute_aux_potential(P, V, smooth) result(this)
  !
  ! Construct auxiliary potential from equi-potential contour(s).
  !
  ! **Parameters:**
  ! :P:        Polygonal representations of equi-potential contour(s)
  ! :V:        Value of potential on each contour
  ! :smooth:   Smoothing parameter for point charge potential
  !
  use moose_math, only: mdgesv
  type(polygon2d), intent(in) :: P(:)
  real(real64),    intent(in) :: V(size(P)), smooth
  type(aux_potential)         :: this

  type(rmesh)  :: support
  real(real64), allocatable :: Mtmp(:,:)
  real(real64) :: dx(2), t1, t2
  integer :: i, info, ip, is, j, n(size(P))


  call cpu_time(t1)
  ! initialize output
  call init_mfunc2d(this, 1)
  n = 0
  do ip=1,size(P)
     n(ip) = P(ip)%segments()
  enddo
  this%n = sum(n)
  allocate (this%q(this%n), source=0.d0)


  ! set target potential and reference points on centers of polygon segments
  allocate (this%x(2,this%n), this%w(this%n), source=0.d0)
  i = 0
  do ip=1,size(P)
     this%q(i+1:i+n(ip)) = V(ip)

     do is=0,n(ip)-1
        i = i + 1
        this%x(:,i) = (P(ip)%node(is) + P(ip)%node(is+1)) / 2
        this%w(i) = smooth * norm2(P(ip)%node(is+1) - P(ip)%node(is))
     enddo
  enddo


  ! calculate normalized potential at center of segments
  allocate (Mtmp(this%n,this%n), source=0.d0)
  ! diagonal element
  do i=1,this%n
     Mtmp(i,i) = -log(this%w(i))
  enddo
  ! off-diagonal elements
  do i=1,this%n
     do j=i+1,this%n
         dx = this%x(:,i) - this%x(:,j)
         Mtmp(i,j) = -log(sqrt(sum(dx**2)) + this%w(j))
         Mtmp(j,i) = -log(sqrt(sum(dx**2)) + this%w(i))
     enddo
  enddo
  call cpu_time(t2)
  print *, "computing matrix: ", t2 - t1, "s"


  ! solve for dummy charges
  call cpu_time(t1)
  call mdgesv(Mtmp, this%q, info)
  if (info /= 0) call ERROR("dgesv failed", "bspline2d_potential", info)
  deallocate (Mtmp)
  call cpu_time(t2)
  print *, "solving matrix: ", t2 - t1, "s"

  end function compute_aux_potential
  !-----------------------------------------------------------------------------


! type-bound procedues (aux_potential):
  !-----------------------------------------------------------------------------
  function aux_potential_eval(this, x) result(val)
  class(aux_potential), intent(in) :: this
  real(real64),         intent(in) :: x(this%ndim)
  real(real64)                     :: val

  real(real64) :: dx(2)
  integer :: k


  val = 0.d0
  do k=1,this%n
     dx = this%x(:,k) - x
     val = val - this%q(k) * log(sqrt(sum(dx**2)) + this%w(k))
  enddo

  end function aux_potential_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function aux_potential_deriv(this, x) result(D)
  class(aux_potential), intent(in) :: this
  real(real64),         intent(in) :: x(this%ndim)
  real(real64)                     :: D(this%ndim)

  real(real64) :: dx(2), r
  integer :: k


  D = 0.d0
  do k=1,this%n
     dx = this%x(:,k) - x
     r = sqrt(sum(dx**2))
     D = D + this%q(k) / (r + this%w(k)) / r * dx
  enddo

  end function aux_potential_deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function aux_potential_hessian(this, x) result(H)
  class(aux_potential), intent(in) :: this
  real(real64),         intent(in) :: x(this%ndim)
  real(real64)                     :: H(this%ndim, this%ndim)


  !print *, "aux_potential_hessian is not implemented"
  !stop
  H = 0.d0

  end function aux_potential_hessian
  !-----------------------------------------------------------------------------
! aux_potential ================================================================


! type-bound procedues (qwork):
  !-----------------------------------------------------------------------------
  function contour(this, dim, i) result(x)
  !
  ! pointer to mesh contour at index *i* along dimension *dim*
  !
  class(qwork),    intent(in) :: this
  integer,         intent(in) :: dim, i
  real(real64), pointer       :: x(:,:)


  nullify(x)
  if (dim < 1  .or.  dim > 2) return


  if (dim == 1) then
     x => this%x(i,:,:)
  elseif (dim == 2) then
     x => this%x(:,i,:)
  endif

  end function contour
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function contour_length(this, dim, i) result(l)
  !
  ! length of mesh contour at index *i* along dimension *dim*
  !
  use moose_math, only: diff
  class(qwork), intent(in) :: this
  integer,      intent(in) :: dim, i
  real(real64)             :: l


  l = sum(norm2(diff(this%contour(dim, i), dim=1), dim=2))

  end function contour_length
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function delta(this, dim, i, j)
  use moose_error
  use moose_utils, only: str
  class(qwork), intent(in) :: this
  integer,      intent(in) :: dim, i, j
  real(real64)             :: delta


  select case(dim)
  case(1)
     delta = norm2(this%x(i+1,j,:) - this%x(i,j,:))

  case(2)
     delta = norm2(this%x(i,j+1,:) - this%x(i,j,:))

  case default
     call ERROR("invalid dim = "//str(dim), "uqwork%delta")
  end select

  end function delta
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine insert_nodes(this, dim, i, n)
  !
  ! insert *n* nodes after position *i* along dimension *dim*
  !
  use moose_error
  class(qwork), intent(inout) :: this
  integer,      intent(in   ) :: dim, i, n

  real(real64), pointer :: x(:,:,:)
  integer :: nu, nv


  nu = this%n(1)
  nv = this%n(2)

  select case (dim)
  case(1)
     nu = nu + n
     allocate (x(0:nu-1, 0:nv-1, 2), source=0.d0)
     x(:i, :, :) = this%x(:i, :, :)
     x(i+n+1:, :, :) = this%x(i+1:, :, :)

  case(2)
     nv = nv + n
     allocate (x(0:nu-1, 0:nv-1, 2), source=0.d0)
     x(:, :i, :) = this%x(:, :i, :)
     x(:, i+n+1:, :) = this%x(:, i+1:, :)

  case default
     call ERROR("invalid idim in qwork%insert_nodes")
  end select
  this%x => x
  this%u(0:,0:) => this%x(:,:,1)
  this%v(0:,0:) => this%x(:,:,2)
  this%n = [nu, nv]
  call this%metadata%set("NODES", this%n)

  end subroutine insert_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine submesh(this, i1, i2, C2, Q, container, alpha, lambda, smooth, wmax, min_width, max_squeeze)
  !
  ! construct quasi-orthogonal sub-mesh between position *i1* and *i2* where
  ! where the geometry of the latter is determined by *C2*.
  !
  use moose_error
  use moose_utils,      only: str, user_option
  use moose_algorithms, only: reverse_array
  use moose_math
  use moose_mfunc
  use moose_geometry,   only: sample_mfunc2d, hypersurf2d
  use moose_fieldline
  class(qwork), intent(inout) :: this
  integer,      intent(in   ) :: i1, i2
  class(curve), intent(in   ) :: C2
  class(qfunc), intent(in   ), optional :: Q
  integer,      intent(in   ), optional :: container
  real(real64), intent(in   ), optional :: alpha, lambda, smooth, wmax, min_width, max_squeeze

  type(interp_qfunc) :: kappa_qfunc
  type(polygon2d) :: P(3)
  type(hypersurf2d), target :: boundary
  type(aux_potential) :: dummyU
  type(fieldline_driver) :: G
  real(real64) :: l1, l2, uP(3), ulevels(min(i1,i2):max(i1,i2)), epsabs, hstart, dmax, minr, xb(0:this%n(2)-1)
  integer :: i, idir, j, n, madd, mpoints


  ! 1. generate polygonal representation of contours
  ! 1.1. i1-th mesh contour
  P(1) = polygon2d(transpose(this%x(i1,:,:)));   l1 = P(1)%length()
  if (.not.P(1)%is_closed()) then
     print *, "first point: ", this%x(i1, 0, :)
     print *, "last point:  ", this%x(i1, this%n(2)-1, :)
     call ERROR("mesh contour is not closed at i = "//str(i1), "qwork%submesh")
  endif
  dmax = maxval(norm2(diff(this%x(i1,:,:), dim=1), dim=2))


  ! 1.2. geometry for i2-th mesh contour
  mpoints = points_per_segment * (this%n(2)-1) + 1
  kappa_qfunc = C2%curvature_qfunc(mpoints, alpha, lambda)
  P(2) = polygon2d(C2%discretization(mpoints, kappa_qfunc));   l2 = P(2)%length()
  boundary = hypersurf2d(P(2))


  ! 1.3. container contour
  n = 2
  if (present(container)) then
  if (container > 0) then
     P(3) = container_polygon(P(1:2), container)
     n = 3
  endif
  endif


  ! 2. construct auxiliary potential
  uP = [1.d0, 2.d0, 0.d0];   if (l1 < l2) uP = [2.d0, 1.d0, 0.d0]
  dummyU = aux_potential(P(1:n), uP(1:n), user_option(1.d0, smooth))
  if (debug_dummyU /= "") then
     do i=1,n
        call sample_mfunc2d(dummyU, trim(debug_dummyU)//"_poly"//str(i), P(i))
     enddo
     !call sample_mfunc2d(dummyU, trim(debug_dummyU), 512, 512)
     call sample_mfunc2d(dummyU, trim(debug_dummyU), 512, 512, 10.24d0, 15.26d0, -5.844d0, 5.844d0)
  endif


  ! 3. follow gradient lines in auxiliary potential
  ! 3.0. set potential levels for mesh contours
  ulevels = linspace(uP(1), uP(2), size(ulevels))
  if (present(Q)) then
     ulevels = uP(1) + (uP(2)-uP(1)) * Q%qquantiles(size(ulevels)-1)
  endif

  ! 3.1. set direction
  idir = 1
  if (i2 < i1) then
     ulevels = reverse_array(ulevels)
     idir = -1
  endif
  ! add margin at i2 (stop at P2 instead)
  ulevels(i2) = ulevels(i2) + 0.1d0 * (ulevels(i2) - ulevels(i1))

  ! 3.2. generate gradient lines
  epsabs = sqrt(product(dummyU%ub - dummyU%lb)) * 1.d-8
  hstart = 1.d-2 / abs(i2-i1)
  G = gradline_driver(dummyU, CARTESIAN2D, YMOD, hstart, "dopr5", epsabs, 0.d0, boundary=boundary)
  G%stepper%hmin = epsabs
  G%stepper%hmax = hstart
  do j=0,this%n(2)-1
     call set_gradline(j, i1, i2)
  enddo


  ! 4. add additional nodes if resolution is too coarse in some locations
  call add_additional_nodes(user_option(1.d0, wmax) * dmax)


  ! 5. homogenize node spacing
  minr = sqrt(1.d0 / this%n(2))
  if (present(min_width)) then
     if (min_width >= 0.d0) minr = min_width
  endif
  do i=i1,i2,idir
     if (minr > 0.d0) call adjust_node_spacing(this, i, 1, minr, user_option(0.8d0, max_squeeze))
  enddo

  contains
  !.............................................................................
  subroutine set_gradline(j, i1, i2)
  integer, intent(in) :: j, i1, i2

  real(real64) :: u, utmp(min(i1,i2):max(i1,i2)), x(2)
  integer :: i, ii, ierr


  x = this%x(i1, j, :)
  u = dummyU%eval(x)
  utmp(i1) = u

  call G%reset()
  do i=i1+idir,i2,idir
     ierr = G%evolve(u, ulevels(i), x)
     ! TODO: at i2, check for intersection with P2, save local coordinates
     if (ierr /= 0  .and.  ierr /= INTERSECT_BOUNDARY) then
        print *, "i1, i, j = ", i1, i, j
        print *, " i,    x(i),         y(i),         u(i),         utarget(i)"
        do ii=i1,i-idir,idir
           print 9001, ii, this%x(ii,j,:), utmp(ii), ulevels(ii)
        enddo
        call this%savetxt("ERROR_QMESH")
        call ERROR("gradline_driver%apply failed", error_code=ierr)
     endif
     if (ierr == INTERSECT_BOUNDARY) xb(j) = G%ub(1)
     this%x(i,j,:) = x
     utmp(i) = u
  enddo
 9001 format(i6,4(2x,e12.6))

  end subroutine set_gradline
  !.............................................................................
  subroutine add_additional_nodes(dmax)
  use moose_analysis, only: interp, pchip
  real(real64), intent(in) :: dmax

  type(interp) :: smap, kmap
  real(real64), allocatable :: pvalues(:), svalues(:)
  real(real64) :: d, x0(2), x1(2), t, ds, s0, s1
  integer :: j, jj, nP2, madd


  ! map curvature weighted parametrization to equal arc length parametrixation
  allocate (pvalues, source = kappa_qfunc%pvalues())
  allocate (svalues, source = linspace(0.d0, 1.d0, size(pvalues)))
  smap = pchip(pvalues, svalues)
  kmap = pchip(svalues, pvalues)
  nP2 = P(2)%segments()


  madd = 0
  do j=this%n(2)-2,0,-1
     ! find max. edge length
     d = 0.d0
     do i=i1,i2,idir
        d = max(d, this%delta(2, i, j))
     enddo

     ! arc length along C2
     s0 = smap%eval(xb(j) / nP2)
     s1 = smap%eval(xb(j+1) / nP2);   if (s1 < s0) s1 = s1 + 1.d0
     ds = s1 - s0
     d = max(d, l2 * ds)
     if (d <= dmax) cycle

     ! insert interpolated nodes
     madd = int(2*d/dmax)
     call this%insert_nodes(2, j, madd)
     do jj=j+1,j+madd
        this%x(:,jj,:) = this%x(:,j,:) + (this%x(:,j+madd+1,:) - this%x(:,j,:)) * (jj - j) / (madd + 1)
        t = mod(s0 + (s1 - s0) * (jj - j) / (madd + 1), 1.d0)
        this%x(i2,jj,:) = P(2)%interpolate(kmap%eval(t) * nP2)
     enddo
  enddo

  end subroutine add_additional_nodes
  !.............................................................................
  end subroutine submesh
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  function qmesh_distance_generator(P, m1, m2, hcontour, hsupport, margin) result(this)
  !
  ! construct qmesh with resolution (*m1*, *m2*) in annular domain between polygons *P(1)* and *P(2)*
  !
  use moose_error
  use moose_utils, only: str, ordinal
  use moose_geometry
  use moose_contours
  type(polygon2d), intent(in) :: P(2)
  integer,         intent(in) :: m1, m2
  real(real64),    intent(in), optional :: hcontour, hsupport, margin
  type(qmesh)                 :: this

  type(bspline2d) :: B(0:2)
  type(contour)   :: C
  type(rmesh)     :: mesh
  real(real64)    :: hcontour_, hsupport_, margin_, lb(2), ub(2), dbmax, x0(2), x1(2), x2(2), v(2), w1, w2
  integer :: i, istat, iu, sgn1, sgn2


  ! allocate new mesh and initialize lower boundary with P(1)
  this = qmesh(m1, m2)
  call set_qmesh_contour_from_polygon2d(this, 0, 1, P(1))


  ! compute distance functions
  margin_ = 1.d-2;   if (present(margin)) margin_ = margin
  call P(2)%get_bounding_box(lb, ub, margin_)
  dbmax = maxval(ub-lb)
  hsupport_ = dbmax / 128;   if (present(hsupport)) hsupport_ = hsupport
  mesh = rmesh_bbox(P, hsupport_, margin=margin_)
  B(1) = dummyU_polygon2d_distance(mesh, P(1))
  B(2) = dummyU_polygon2d_distance(mesh, P(2))
  B(0) = B(1)


  ! find reference points on P(1) and P(2)
  x1   = P(1)%node(0)
  v    = P(1)%secant_normal(0)
  sgn1 = int(sign(1.d0, P(1)%area()))
  sgn2 = int(sign(1.d0, P(2)%area()))
  x0   = x1-sgn1*dbmax*v
  if (.not.P(2)%intersects_segment(x1, x0, x2)) then
     print *, "x1   = ", x1
     print *, "xmax = ", x0
     open  (newunit=iu, file="XSECT_ERROR")
     write (iu, *) x1
     write (iu, *) x0
     close (iu)
     call ERROR("cannot find reference point on P(2), see XSECT_ERROR")
  endif


  ! generate qmesh from contours between P(1) and P(2)
  hcontour_ = hsupport_;   if (present(hcontour)) hcontour_ = hcontour
  do i=1,m1-1
     w2 = 1.d0 * i / (m1-1)
     w1 = 1.d0 - w2
     B(0)%bcoef = w1 * sgn1 * B(1)%bcoef  +  w2 * sgn2 * B(2)%bcoef

     ! find reference points on i-th contour between P(1) and P(2)
     x0 = x2;   if (i < m1-1) x0 = line_search(B(0), x1, x2, 0.d0, istat)
     if (istat > 0) then
        print *, "x1 = ", x1
        print *, "x2 = ", x2
        call sample_mfunc2d(B(0), "DUMMY_POTENTIAL", size(B(0)%uknot), size(B(0)%vknot))
        call ERROR("line search for reference point "//str(i)//" failed")
     endif

     ! generate i-th contour
     C = contour(B(0), x0, sgn1, hcontour_, istat)
     if (istat /= 0) call ERROR("tracing of "//ordinal(i)//" contour faild")
     this%x(i,:,:) = transpose(interp_contour_discretization(C, m2))
  enddo

  end function qmesh_distance_generator
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_qmesh_distance_contour_generator_curve(this, i0, i1, C1, eps, Rmax, hcontour, hsupport, bbox_margin)
  !
  ! i0:  mesh-index for initial contour
  ! i1:  mesh-index for final contour
  ! C1:  geometry for final mesh contour
  !
  ! **Optional parameters:**
  !
  ! :eps:           Max. displacement from arc segment for polygonal representation of *C1*
  ! :Rmax:          Max. curvature radius taken into account for discretization of *C1*
  ! :hcontour:      Spatial resolution for tracing of mesh contour
  ! :hsupport:      Spatial resolution for dummy potential support mesh
  ! :bbox_margin:   Relative margin for bounding box of dummy potential
  !
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i0, i1
  class(curve),    intent(in   ) :: C1
  real(real64),    intent(in   ), optional :: eps, Rmax, hcontour, hsupport, bbox_margin

  real(real64)    :: l0, eps_, hcontour_, Rmax_
  type(polygon2d) :: P(0:1)


  ! polygonal representation of initial contour
  P(0) = polygon2d_qmesh_contour(this, i0, 1);   l0 = P(0)%length()


  ! set default values for optional arguments
  hcontour_ = l0 / this%n(2);   if (present(hcontour)) hcontour_ = hcontour
  eps_ = hcontour_ / 100.d0;  if (present(eps)) eps_ = eps
  Rmax_ = hcontour_ * 10.d0;   if (present(Rmax)) Rmax_ = Rmax


  ! construct polygonal representation of C1 for implementation of mesh generator
  P(1) = polygon2d(C1%polygon(eps_, Rmax_))
  call aux_qmesh_distance_contour_generator_poly2d(this, i0, i1, P, hcontour, hsupport, bbox_margin)

  end subroutine aux_qmesh_distance_contour_generator_curve
  !-----------------------------------------------------------------------------
  subroutine aux_qmesh_distance_contour_generator_poly2d(this, i0, i1, P, hcontour, hsupport, bbox_margin)
  use moose_error
  use moose_utils, only: str, ordinal
  use moose_contours
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i0, i1
  type(polygon2d), intent(in   ) :: P(0:1)
  real(real64),    intent(in   ), optional :: hcontour, hsupport, bbox_margin

  type(bspline2d) :: U(0:1), Ui
  type(rmesh)     :: mesh
  type(contour)   :: C
  real(real64)    :: bbox_margin_, hcontour_
  real(real64)    :: dmax, l0, v(2), w(0:1), x0(2), x1(2), xi(2)
  integer :: i, idir, istat, iu, sgn(0:1)


  ! set default values for optional arguments
  l0 = P(0)%length()
  bbox_margin_ = 1.d-2;   if (present(bbox_margin)) bbox_margin_ = bbox_margin
  hcontour_ = l0 / this%n(2);   if (present(hcontour)) hcontour_ = hcontour


  ! polygonal representation of final mesh contour
  mesh = rmesh_bbox(P, hsupport, margin=bbox_margin_)
  dmax = max(mesh%u(size(mesh%u)-1) - mesh%u(0), mesh%u(size(mesh%v)-1) - mesh%v(0))


  ! construct dummy potential from distance function
  do i=0,1
     sgn(i) = int(sign(1.d0, P(i)%area()))
     U(i) = dummyU_polygon2d_distance(mesh, P(i))
     U(i)%bcoef = U(i)%bcoef * sgn(i)
  enddo
  Ui = U(0)


  ! find reference points on P(0) and P(1)
  x0   = P(0)%node(0)
  v    = P(0)%secant_normal(0)
  idir = 1;   if (i1 < i0) idir = -1
  xi   = x0 - sgn(0) * idir * dmax * v
  if (.not.P(1)%intersects_segment(x0, xi, x1)) then
     call this%savetxt("ERROR_MESH")
     print *, "x0   = ", x0
     print *, "xmax = ", xi
     open  (newunit=iu, file="XSECT_ERROR")
     write (iu, *) x0
     write (iu, *) xi
     close (iu)
     call ERROR("cannot find reference point on P(1), see XSECT_ERROR")
  endif


  ! generate qmesh from contours between P(0) and P(1)
  do i=i0+idir,i1,idir
     w(1) = 1.d0 * (i-i0) / (i1-i0)
     w(0) = 1.d0 - w(1)

     ! combine U(0) and U(1) for i-th mesh contour
     Ui%bcoef = w(0) * U(0)%bcoef  +  w(1) * U(1)%bcoef

     ! find reference points on i-th contour between P(0) and P(1)
     xi = x1;   if (i /= i1) xi = line_search(Ui, x0, x1, 0.d0, istat)
     if (istat > 0) then
        print *, "i, w0, w1 = ", i, w
        print *, "x0 = ", x0
        print *, "x1 = ", x1
        call sample_mfunc2d(Ui, "DUMMY_POTENTIAL", size(Ui%uknot), size(Ui%vknot))
        call ERROR("line search for reference point "//str(i)//" failed")
     endif

     ! generate i-th contour
     C = contour(Ui, xi, sgn(0), hcontour_, istat)
     if (istat /= 0) then
        call sample_mfunc2d(Ui, "DUMMY_POTENTIAL", size(Ui%uknot), size(Ui%vknot))
        call C%savetxt("CONTOUR_ERROR")
        call ERROR("tracing of "//ordinal(i)//" contour failed")
     endif
     this%x(i,:,:) = transpose(interp_contour_discretization(C, this%n(2)))
  enddo

  end subroutine aux_qmesh_distance_contour_generator_poly2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_qmesh_distance_gradient_generator(this, i0, i1, C1, eps, Rmax, hcontour, hsupport, bbox_margin, minr, alpha)
  !
  ! i0:  mesh-index for initial contour
  ! i1:  mesh-index for final contour
  ! C1:  geometry for final mesh contour
  !
  ! **Optional parameters:**
  !
  ! :eps:           Max. displacement from arc segment for polygonal representation of *C1*
  ! :Rmax:          Max. curvature radius taken into account for discretization of *C1*
  ! :hcontour:      Spatial resolution for tracing of mesh contour
  ! :hsupport:      Spatial resolution for dummy potential support mesh
  ! :bbox_margin:   Relative margin for bounding box of dummy potential
  !
  use moose_error
  use moose_utils, only: str, ordinal
  use moose_rlist
  use moose_math
  use moose_contours
  use moose_fieldline
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i0, i1
  class(curve),    intent(in   ) :: C1
  real(real64),    intent(in   ), optional :: eps, Rmax, hcontour, hsupport, bbox_margin, minr, alpha

  real(real64), parameter :: hstart = 1.d-5

  type(polygon2d) :: P(0:1)
  type(bspline2d) :: U(0:1), Ui
  type(rmesh)     :: mesh
  type(rlist)     :: T
  type(fieldline_driver) :: G

  real(real64)    :: bbox_margin_, eps_, hcontour_, Rmax_, minr_, alpha_
  real(real64)    :: dmax, l0, w(0:1), x0(2), tmp(3)
  integer :: i, idir, istat, iu, j, sgn(0:1)


  ! polygonal representation of initial contour
  P(0) = polygon2d_qmesh_contour(this, i0, 1);   l0 = P(0)%length()


  ! set default values for optional arguments
  bbox_margin_ = 1.d-2;   if (present(bbox_margin)) bbox_margin_ = bbox_margin
  hcontour_ = l0 / this%n(2);   if (present(hcontour)) hcontour_ = hcontour
  eps_ = hcontour_ / 100.d0;  if (present(eps)) eps_ = eps
  Rmax_ = hcontour_ * 10.d0;   if (present(Rmax)) Rmax_ = Rmax
  minr_  = 1.d0 / this%n(2);   if (present(minr)) minr_ = minr
  alpha_ = 0.8d0;   if (present(alpha)) alpha_ = alpha


  ! polygonal representation of final mesh contour
  P(1) = polygon2d(C1%polygon(eps_, Rmax_))
  mesh = rmesh_bbox(P, hsupport, margin=bbox_margin_)
  dmax = max(mesh%u(size(mesh%u)-1) - mesh%u(0), mesh%u(size(mesh%v)-1) - mesh%v(0))


  ! construct dummy potential from distance function
  do i=0,1
     sgn(i) = int(sign(1.d0, P(i)%area()))
     U(i) = dummyU_polygon2d_distance(mesh, P(i))
     U(i)%bcoef = U(i)%bcoef * sgn(i)
  enddo
  Ui = U(0)


  ! generate qmesh from contours between P(0) and P(1)
  idir = 1;   if (i1 < i0) idir = -1
  G = gradline_driver(Ui, CARTESIAN2D, YMOD, hstart, "dopr5", 1.d-8, 0.d0)
  do i=i0+idir,i1,idir
     w(1) = 1.d0 * (i-i0) / (i1-i0)
     w(0) = 1.d0 - w(1)

     ! combine U(0) and U(1) for i-th mesh contour
     Ui%bcoef = w(0) * U(0)%bcoef  +  w(1) * U(1)%bcoef

     do j=0,this%n(2)-1
        x0 = this%node(i-idir,j)

        call G%reset()
        T = G%trace(Ui%eval(x0), x0, 0.d0)

        tmp = T%element(-1);   this%x(i,j,:) = tmp(1:2)
     enddo

     ! expand nodes which are too close
     call adjust_node_spacing(this, i, 1, minr_, alpha_)
  enddo



  end subroutine aux_qmesh_distance_gradient_generator
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function interp_contour_discretization(C, m, Q) result(x)
  !
  ! generate discretization of interpolated contour *C*
  !
  use moose_analysis, only: qfunc
  use moose_geometry, only: interp_curve, contour
  class(contour), intent(in) :: C
  integer,        intent(in) :: m
  class(qfunc),   intent(in), optional :: Q
  real(real64)               :: x(2, 0:m-1)

  type(interp_curve) :: tmp


  tmp = C%interp()
  x   = tmp%discretization(m, Q)
  call tmp%free()

  end function interp_contour_discretization
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function dummyU_polygon2d_distance(mesh, P) result(this)
  !
  ! construct dummy potential on *mesh* from distance to polygon *P*
  !
  type(rmesh),     intent(in) :: mesh
  type(polygon2d), intent(in) :: P
  type(bspline2d)             :: this

  real(real64), allocatable :: d(:,:)
  real(real64) :: x(2)
  integer      :: i, j, n(2)


  n = mesh%n
  allocate (d(0:n(1)-1, 0:n(2)-1), source=0.d0)
  do i=0,n(1)-1
  do j=0,n(2)-1
     x = [mesh%u(i), mesh%v(j)]
     d(i,j) = P%get_distance(x)
  enddo
  enddo
  this = bspline2d(mesh%u, mesh%v, d)

  end function dummyU_polygon2d_distance
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_qmesh_contour_from_polygon2d(this, i, idim, P, Q)
  !
  ! set *i*-th mesh contour in dimention *idim* from polygon *P*
  !
  use moose_geometry, only: bspline_polygon
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i, idim
  type(polygon2d), intent(in   ) :: P
  class(qfunc),    intent(in   ), optional :: Q


  call set_qmesh_contour(this, i, idim, bspline_polygon(P), Q)

  end subroutine set_qmesh_contour_from_polygon2d
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine adjust_node_spacing(this, i, idim, minr, alpha)
  !
  ! adjust nodes on *i*-th mesh contour if they are very close to each other
  !
  ! **Parameters:**
  !
  ! :minr:      Minimal distance between nodes relative to mean cell width (1: equidistant)
  !
  ! :alpha:     Allow segments to be squeezed up to this fraction (0: stiff, 1: limited by minr)
  !
  use moose_error
  use moose_utils, only: str
  use moose_grids, only: cgrid
  use moose_geometry
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: i, idim
  real(real64), intent(in   ) :: minr, alpha

  type(polygon2d)     :: P
  type(bspline_curve) :: C
  type(cgrid)         :: tmp
  real(real64), allocatable :: t(:)
  real(real64) :: min_dt
  integer :: ierr, iu, j, n


  ! initialize polygonal representation of mesh contour and spacings array
  P = polygon2d_qmesh_contour(this, i, idim)
  n = P%segments()
  allocate (t(0:n), source=P%accumulated_lengths())


  ! compute adjusted node spacings
  min_dt = t(n) / n * minr
  call aux_adjust_node_spacing(t, min_dt, alpha, ierr)
  if (ierr /= 0) then
     print *, "min_dt = ", min_dt
     print *, "alpha  = ", alpha
     tmp = P%cgrid()
     call tmp%savetxt("ERROR_NODE_ADJUSTMENT")
     call ERROR("node adjustment failed on mesh contour "//str(i), error_code=ierr)
  endif


  ! update mesh contour based on new spacings
  C = bspline_polygon(P)
  do j=0,n
     call P%set_node(j, C%eval(t(j)))
  enddo
  call set_qmesh_contour(this, i, idim, P)
  call C%free()
  call P%free()

  end subroutine adjust_node_spacing
  !-----------------------------------------------------------------------------
  subroutine aux_adjust_node_spacing(t, min_dt, alpha, ierr)
  !
  ! adjust intervals between nodes *t* such that each segment is at last *min_dt*,
  ! but do not compress other segments by more than a fator of *alpha*
  !
  ! error codes:
  !    1   estimated adjustment exceeds available space on left hand side
  !    2   estimated adjustment exceeds available space on right hand side
  !    3   no further squeezing of intervals possible
  !    9   min_dt exceeds t-range / size(t)
  !
  use moose_error
  real(real64), intent(inout) :: t(0:)
  real(real64), intent(in   ) :: min_dt, alpha
  integer,      intent(  out) :: ierr

  real(real64) :: d, dab, h, t0
  integer :: i, ia(0:size(t)-1), ib(0:size(t)-1), ik, ik1, im1, ip1, k, n, nab


  if (min_dt > abs(t(ubound(t,1)) - t(0)) / (size(t)-1)) then
     ierr = 9
     return
  endif


  ierr = 0
  ! 1. identify index range [ia, ib] for which node adjustment is required
  n  = size(t)-1
  t0 = t(0);   t = t-t0;   d = t(n)
  i  = 1
  k  = 0;   ia = 0
  do
     if (i == ia(0)) exit
     im1 = modulo(i-1,n);   ip1 = modulo(i+1,n)

     ! identify first node for left-shift
     if (dt(i,im1) >= min_dt  .and.  dt(ip1,i) < min_dt) then
        ia(k) = i
        do
           i = modulo(i+1,n)
           if (i == ia(k)) exit
           im1 = modulo(i-1,n);   ip1 = modulo(i+1,n)

           ! identify last node for right-shift
           if (dt(ip1,i) >= min_dt) then
              ib(k) = i
              exit
           endif
        enddo
        !print *, "ia, ib = ", ia(k), ib(k)
        k = k + 1
     endif
     i = modulo(i+1,n)
  enddo


  ! 2. merge adjacent index ranges, if necessary
  ik = 0
  do
     if (ik == k) exit

     ik1 = modulo(ik+1,k)
     h = (required_adjustment(ik) + required_adjustment(ik1)) / 2
     !print *, ik, ib(ik), ia(ik1), h, squeezable_space(ib(ik), ia(ik1), 1)

     ! not enough space
     if (h > squeezable_space(ib(ik), ia(ik1), 1)) then
        !print *, "merging ", ik, " and ", ik1
        if (ik1 == 0) then
           ia(0) = ia(ik)
        else
           ib(ik:k-2) = ib(ik+1:k-1)
           if (ik < k-2) ia(ik+1:k-2) = ia(ik+2:k-1)
        endif
        k = k - 1
     else
        ik = ik + 1
     endif
  enddo


  ! 3. adjust nodes
  do ik=0,k-1
     ! width of adjusted domain
     h = required_adjustment(ik)
     !print *, ik, t(ia(ik)), t(ib(ik)), t(ia(modulo(ik+1,k))), t(ib(modulo(ik-1,k)))
     !print *, ik, "h = ", h, ib(ik) - ia(ik), dsmin, dt(ib(ik), ia(ik))
     !print *, ik, ia(ik), ib(ik), t(ia(ik)), t(ib(ik))
     !print *, h, t(ia(ik)) - h/2, t(ib(ik)) + h/2

     ! verif that the is enough space to the left
     !print *, i, "left:  ", dt(ia(i), ib(modulo(i-1,k)))
     if (dt(ia(ik), ib(modulo(ik-1,k))) < h/2) then
        ierr = 1
        return
        !call ERROR("not enough space for adjustment (left)")
     endif
     call squeeze_nodes(ia(ik), ib(modulo(ik-1,k)), -1, ierr);   if (ierr /= 0) return


     ! verif that the is enough space to the right
     !print *, i, "right: ", dt(ia(modulo(i+1,k)), ib(i))
     if (dt(ia(modulo(ik+1,k)), ib(ik)) < h/2) then
        ierr = 2
        return
        ! call ERROR("not enough space for adjustment (right)")
     endif
     call squeeze_nodes(ib(ik), ia(modulo(ik+1,k)), 1, ierr);   if (ierr /= 0) return


     ! update nodes [ia, ib]
     dab = dt(ib(ik), ia(ik))
     nab = di(ib(ik), ia(ik))
     !print *, ia(ik), ib(ik), dab, nab
     do i=ia(ik)+1,ia(ik)+nab-1
        t(modulo(i,n)) = modulo(t(ia(ik)) + dab * (i-ia(ik)) / nab, d)
        !print *, t(modulo(i,n))
     enddo
  enddo


  ! 4. update t(n) and restore offset t0
  t(n) = t(0) + d;   t = modulo(t + t0, d)

  contains
  !.............................................................................
  function di(i2, i1)
  integer, intent(in) :: i1, i2
  integer             :: di


  di = modulo(i2 - i1, n)

  end function di
  !.............................................................................
  function dt(i2, i1)
  integer, intent(in) :: i1, i2
  real(real64)        :: dt


  dt = dtmod(t(i2), t(i1), 1)

  end function dt
  !.............................................................................
  function dtmod(t2, t1, idir)
  real(real64), intent(in) :: t2, t1
  integer,      intent(in) :: idir
  real(real64)             :: dtmod


  dtmod = modulo(idir*(t2-t1), d)

  end function dtmod
  !.............................................................................
  function required_adjustment(ik) result(h)
  integer, intent(in) :: ik
  real(real64)        :: h


  h = di(ib(ik),ia(ik)) * min_dt  -  dt(ib(ik), ia(ik))

  end function required_adjustment
  !.............................................................................
  function squeezable_space(i1, i2, idir)
  integer, intent(in) :: i1, i2, idir
  real(real64)        :: squeezable_space

  real(real64) :: dti, ti1
  integer :: i


  squeezable_space = 0.d0
  i = i1
  do
     if (i == i2) exit

     ti1 = t(modulo(i+idir,n))
     dti = abs(ti1 - t(i))

     squeezable_space = squeezable_space + min(dti-min_dt, (1-alpha)*dti)
     i = modulo(i+idir,n)
  enddo

  end function squeezable_space
  !.............................................................................
  subroutine squeeze_nodes(ia, ib, idir, ierr)
  integer, intent(in   ) :: ia, ib, idir
  integer, intent(  out) :: ierr

  real(real64) :: dh, dti, ti1
  integer :: i, i1, idir1


  i = ia;   dh = h/2
  !print *, "squeezing nodes ", ia, ib, idir, dh
  do
     ! get next node and (initial) distance to it, and push this node
     i1   = modulo(i+idir,n);   ti1 = t(i1)
     dti  = dtmod(ti1,t(i),idir)
     t(i) = modulo(t(i) + idir*dh, d)
     !print *, i, dh, t(i), dti

     ! remaining adjustment can be absorbed by this segment
     if (dh < dti*alpha  .and.  dtmod(ti1,t(i),idir) >= min_dt) exit

     ! this is the last segment - but adjustment is not completed yet
     if (i == ib-idir) then
        ierr = 3
        return
        ! call ERROR("no further adjustment possible")
     endif

     ! partial adjustment, move on to next segment
     idir1 = 1;   if (dh > dti) idir1 = -1
     dh = min_dt - idir1*dtmod(ti1,t(i),idir*idir1)
     i  = modulo(i+idir,n)
  enddo

  end subroutine squeeze_nodes
  !.............................................................................
  end subroutine aux_adjust_node_spacing
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rmesh_bbox(P, h, n, margin) result(this)
  !
  ! Construct :class:`Rmesh` with spatial resolution *h* in each direction
  ! across bounding box for set of polygons *P*
  !
  type(polygon2d), intent(in) :: P(:)
  real(real64),    intent(in), optional :: h, margin
  integer,         intent(in), optional :: n(2)
  type(rmesh)                 :: this

  real(real64) :: lb(2), ub(2), xmin, xmax, ymin, ymax
  integer :: i, n1, n2


  if (size(P) == 0) return


  call P(1)%get_bounding_box(lb, ub, margin)
  xmin = lb(1);   xmax = ub(1)
  ymin = lb(2);   ymax = ub(2)
  do i=2,size(P)
     call P(i)%get_bounding_box(lb, ub, margin)
     xmin = min(xmin, lb(1));   xmax = max(xmax, ub(1))
     ymin = min(ymin, lb(2));   ymax = max(ymax, ub(2))
  enddo


  ! explicit mesh resolution is given
  if (present(n)) then
     n1 = n(1)
     n2 = n(2)

  ! cell spacing parameter is given
  elseif (present(h)) then
     n1 = (xmax-xmin) / h + 1
     n2 = (ymax-ymin) / h + 1

  ! fallback
  else
     n1 = P(1)%nnodes()
     do i=2,size(P)
        n1 = n1 * P(1)%nnodes()
     enddo
     n1 = n1**(1.d0/size(P))
     n2 = n1
  endif

  this = rmesh(xmin, xmax, n1, ymin, ymax, n2)

  end function rmesh_bbox
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function bspline2d_potential(P, V, smooth) result(this)
  !
  ! Construct B-Spline potential from equi-potential contour(s).
  !
  ! **Parameters:**
  ! :P:        Polygonal representations of equi-potential contour(s)
  ! :V:        Value of potential on each contour
  ! :smooth:   Smoothing parameter for point charge potential
  !
  use moose_math, only: mdgesv
  type(polygon2d), intent(in) :: P(:)
  real(real64),    intent(in) :: V(size(P)), smooth
  type(bspline2d)             :: this

  type(rmesh)  :: support
  real(real64), allocatable :: x(:,:), Mtmp(:,:), q(:)
  real(real64) :: dx(2)
  integer :: i, info, ip, is, j, k, n(size(P)), nsum


  ! initialize output
  n = 0
  do ip=1,size(P)
     n(ip) = P(ip)%segments()
  enddo
  nsum = sum(n)
  allocate (q(nsum), source=0.d0)


  ! set target potential and reference points on centers of polygon segments
  allocate (x(2,nsum), source=0.d0)
  i = 0
  do ip=1,size(P)
     q(i+1:i+n(ip)) = V(ip)

     do is=0,n(ip)-1
        i = i + 1
        x(:,i) = (P(ip)%node(is) + P(ip)%node(is+1)) / 2
     enddo
  enddo


  ! calculate normalized potential at center of segments
  allocate (Mtmp(nsum,nsum), source=0.d0)
  ! diagonal element
  do i=1,nsum
     Mtmp(i,i) = -log(smooth)
  enddo
  ! off-diagonal elements
  do i=1,nsum
     do j=i+1,nsum
         dx = x(:,i) - x(:,j)
         Mtmp(i,j) = -log(sqrt(sum(dx**2)) + smooth)
         Mtmp(j,i) = Mtmp(i,j)
     enddo
  enddo


  ! solve for dummy charges
  call mdgesv(Mtmp, q, info)
  if (info /= 0) call ERROR("dgesv failed", "bspline2d_potential", info)
  deallocate (Mtmp)


  ! construct B-spline representation
  support = rmesh_bbox(P, smooth, margin=bbox_margin)
  allocate (Mtmp(0:support%n(1)-1, 0:support%n(2)-1), source=0.d0)
  do i=0,support%n(1)-1
     do j=0,support%n(2)-1
        do k=1,nsum
           dx = x(:,k) - support%node(i,j)
           Mtmp(i,j) = Mtmp(i,j) - q(k) * log(sqrt(sum(dx**2)) + smooth)
        enddo
     enddo
  enddo
  this = bspline2d(support%u, support%v, Mtmp)
  deallocate (x, Mtmp)
  call support%free()

  end function bspline2d_potential
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function container_polygon(P, n, margin) result(this)
  use moose_error
  use moose_math
  use moose_utils, only: str
  use moose_geometry, only: shifted_polygon2d, bspline_curve, bspline_polygon
  class(polygon2d), intent(in) :: P(:)
  integer,          intent(in) :: n
  real(real64),     intent(in), optional :: margin
  type(polygon2d)              :: this

  type(bspline_curve) :: C
  real(real64)        :: area, ds, sgn, x(2)
  integer             :: i, m


  ! check if all polygons are closed
  m = size(P)
  do i=1,m
     if (.not.P(i)%is_closed()) then
        call P(i)%savetxt("ERROR_CONTAINER_POLYGON")
        call ERROR("container_polygon called with open polygon "//str(i))
     endif
  enddo


  ! find outermost polygon (NOTE: polygons must not be intersecting each other)
  if (m <= 0) then
     call ERROR("at least 1 polygon is required for container_polygon")

  elseif (m == 1) then
     i = 1

  elseif (m == 2) then
     x = P(1)%node(0)
     if (P(2)%winding_number(x) == 0) then
        i = 1
     else
        i = 2
     endif

  else
     call ERROR("container_polygon not implemented for size(P) > 2")
  endif


  ! set parameters
  area = P(i)%area()
  sgn  = int(sign(1.d0, area))
  ds   = sgn * sqrt(abs(area))
  if (present(margin)) then
     ds = ds * margin
  else
     ds = ds * 0.4d0
  endif


  ! construct container polygon
  C = bspline_polygon(shifted_polygon2d(P(i), -ds))
  this = polygon2d(C%discretization(n))

  end function container_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_submesh(this, i1, i2, C2, Q, container, minr, alpha)
  !
  ! Construct quasi-orthogonal sub-mesh between the *i1*-th and *i2*-th mesh
  ! edges where the geometry of the latter is determined by *C2*.
  !
  use moose_error
  use moose_math
  use moose_algorithms, only: reverse_array
  use moose_mfunc
  use moose_geometry, only: sample_mfunc2d
  use moose_utils, only: str
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: i1, i2
  class(curve), intent(in   ) :: C2
  class(qfunc), intent(in   ), optional :: Q
  integer,      intent(in   ), optional :: container
  real(real64), intent(in   ), optional :: minr, alpha

  type(polygon2d) :: P(0:3)
  type(bspline2d) :: dummyU
  real(real64)    :: smooth, l1, l2, uP(3), ulevels(min(i1,i2):max(i1,i2))
  integer :: i, n, mpoints


  ! polygonal representation of i1-th mesh surface
  P(1) = polygon2d(transpose(this%x(i1,:,:)));   l1 = P(1)%length()
  if (.not.P(1)%is_closed()) then
     print *, "first point: ", this%x(i1,0,:)
     print *, "last point:  ", this%x(i1,this%n(2)-1,:)
     call ERROR("initial mesh contour is not closed")
  endif


  ! polygonal representation of C2
  mpoints = points_per_segment * (this%n(2)-1) + 1
  P(2) = polygon2d(C2%arclength_discretization(mpoints));   l2 = P(2)%length()


  n = 2
  ! add container contour
  if (present(container)) then
  if (container > 0) then
     P(3) = container_polygon(P(1:2), container)
     n = 3
  endif
  endif


  ! construct dummy potential
  uP = [1.d0, 2.d0, 0.d0];   if (l1 < l2) uP = [2.d0, 1.d0, 0.d0]
  smooth = l1 / this%n(2)
  dummyU = bspline2d_potential(P(1:n), uP(1:n), smooth)
  if (debug_dummyU /= "") then
     do i=1,n
        call sample_mfunc2d(dummyU, trim(debug_dummyU)//"_poly"//str(i), P(i))
     enddo
     call sample_mfunc2d(dummyU, trim(debug_dummyU), 512, 512)
  endif


  ! follow gradient lines in dummy potential
  ulevels = linspace(uP(1), uP(2), size(ulevels))
  if (present(Q)) then
     ulevels = uP(1) + (uP(2)-uP(1)) * Q%qquantiles(size(ulevels)-1)
  endif
  if (i2 < i1) ulevels = reverse_array(ulevels)
  call aux_construct_submesh(this, i1, i2, ulevels, dummyU, minr, alpha)

  end subroutine construct_submesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_construct_submesh(this, i1, i2, ulevels, dummyU, minr, alpha)
  use moose_error
  use moose_math
  use moose_fieldline
  class(qmesh),          intent(inout) :: this
  integer,               intent(in   ) :: i1, i2
  real(real64),          intent(in   ) :: ulevels(min(i1,i2):max(i1,i2))
  class(scalar_mfunc2d), intent(in   ) :: dummyU
  real(real64),          intent(in   ), optional :: minr, alpha

  type(fieldline_driver) :: G
  real(real64) :: epsabs, hstart, u, utmp(min(i1,i2):max(i1,i2)), x(2), minr_, alpha_
  integer :: i, idir, ierr, ii, j


  ! follow gradient lines in dummy potential
  idir = 1;   if (i2 < i1) idir = -1
  epsabs = sqrt(product(dummyU%ub - dummyU%lb)) * 1.d-8
  hstart = 1.d-1/abs(i2-i1)
  G = gradline_driver(dummyU, CARTESIAN2D, YMOD, hstart, "dopr5", epsabs, 0.d0)
  G%stepper%hmin = epsabs
  do j=0,this%n(2)-1
     x = this%x(i1,j,:)
     u = dummyU%eval(x)
     utmp(i1) = u

     call G%reset()
     do i=i1+idir,i2,idir
        ierr = G%evolve(u, ulevels(i), x)
        if (ierr /= 0) then
           print *, "i1, i, j = ", i1, i, j
           print *, " i,    x(i),         y(i),         u(i),         utarget(i)"
           do ii=i1,i-idir,idir
              print 9001, ii, this%x(ii,j,:), utmp(ii), ulevels(ii)
           enddo
           call this%savetxt("ERROR_QMESH")
           call ERROR("gradline_driver%apply failed", error_code=ierr)
        endif
        this%x(i,j,:) = x
        utmp(i) = u
     enddo
  enddo
 9001 format(i6,4(2x,e12.6))


  ! expand nodes which are too close
  minr_  = sqrt(1.d0 / this%n(2))
  if (present(minr)) then
     if (minr >= 0.d0) minr_ = minr
  endif
  alpha_ = 0.8d0;   if (present(alpha)) alpha_ = alpha
  do i=i1+idir,i2,idir
     if (minr_ > 0.d0) call adjust_node_spacing(this, i, 1, minr_, alpha_)
  enddo

  end subroutine aux_construct_submesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine interpolate_submesh(this, i1, i2, mdim, Q)
  !
  ! interpolate mesh between i1-th and i2-th surface in mdim-th dimension with
  ! optional node spacing function Q
  !
  use moose_error
  use moose_math
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: i1, i2, mdim
  class(qfunc), intent(in   ), optional :: Q

  real(real64) :: xi(min(i1,i2):max(i1,i2)), x1(2), x2(2)
  integer :: i, i12, i21, j


  i12 = min(i1,i2)
  i21 = max(i1,i2)
  if (present(Q)) then
     xi = Q%qquantiles(abs(i2-i1))
  else
     xi = linspace(0.d0, 1.d0, abs(i2-i1)+1)
  endif


  if (mdim == 1) then
     do j=0,this%n(2)-1
        x1 = this%x(i1,j,:)
        x2 = this%x(i2,j,:)
        do i=i12+1,i21-1
           this%x(i,j,:) = x1 + xi(i) * (x2-x1)
        enddo
     enddo

  elseif (mdim == 2) then
     do j=0,this%n(1)-1
        x1 = this%x(j,i1,:)
        x2 = this%x(j,i2,:)
        do i=i12+1,i21-1
           this%x(j,i,:) = x1 + xi(i) * (x2-x1)
        enddo
     enddo
  else
     call ERROR("mdim = 1 or 2 required", "interpolate_submesh")
  endif

  end subroutine interpolate_submesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_wfunc(t, a, b)
  use moose_error
  real(real64), intent(inout) :: t(:)
  real(real64), intent(in   ) :: a, b

  logical :: wrapped
  integer :: i, itype, n


  n = size(t)
  wrapped = .false.
  do i=2,n-1
     if (t(i) - t(i-1) < 0.d0) then
         if (wrapped) call ERROR("non-monotonic t detected")
         t(i:) = t(i:) + b - a
         wrapped = .true.
     endif
  enddo


  if (t(1) == a) then
     t(n) = b
  else
     t(n) = t(1) + b - a
  endif

  end subroutine aux_wfunc
  !-----------------------------------------------------------------------------
  function make_wfunc(t, a, b) result(this)
  use moose_analysis
  real(real64), intent(in) :: t(:)
  real(real64), intent(in) :: a, b
  type(interp)             :: this

  real(real64) :: tt(size(t))
  integer :: n


  tt = t
  call aux_wfunc(tt, a, b)


  n = size(t)
  this = pchip(linspace(0.d0, 1.d0, n), tt)

  end function make_wfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function normal_ray_wfunc(C, m, P, idir) result(wfunc)
  !
  ! follow normal rays from *C* to *P*, and construct normalized weight
  ! function for *m* blocks.
  !
  ! idir = -1: left hand side of C
  !         1: right hand side of C
  !
  use moose_math,     only: linspace
  use moose_analysis, only: interp, pchip
  use moose_geometry, only: XSECT_RAY
  class(curve),     intent(in) :: C
  integer,          intent(in) :: m
  class(polygon2d), intent(in) :: P
  integer,          intent(in) :: idir
  type(interp)                 :: wfunc

  logical      :: periodic
  real(real64) :: r, s, t, w(0:m), d(2,0:1), x(2), xsect(2), u(2)
  integer      :: i, i0, iu, n


  periodic = C%is_closed
  if (periodic) then
     if (.not.P%is_closed()) then
        call C%savetxt("ERROR_C")
        call P%savetxt("ERROR_P")
        call ERROR("P must be closed in normal_ray_wfunc if C is closed")
     endif
     i0   = 0
  else
     if (P%is_closed()) call ERROR("P must be open in normal_ray_wfunc if C is open")
     i0   = 1
     w(0) = 0.d0
     w(m) = 1.d0
  endif


  ! find intersections of normal rays with P
  do i=i0,m-1
     t = C%a + (C%b-C%a) * i / m
     d = C%deriv(t,1)
     x = d(:,0)
     u = idir * [d(2,1), -d(1,1)]
     call P%intersect(x, u, XSECT_RAY, xsect, r, s, n)
     if (n < 0) then
        open  (newunit=iu, file="ERROR_NORMAL_RAY")
        write (iu, *) x
        write (iu, *) x+u
        close (iu)
        call ERROR_("no intersection found for i-th ray in normal_ray_wfunc failed")
     endif
     w(i) = s + n
  enddo
  w = w / P%segments()
  if (abs(w(0)) < 1.d-99) w(0) = 0.d0
  if (periodic) w(m) = w(0) + 1.d0


  ! construct mapping function
  do i=1,m
     do
        if (w(i) - w(i-1) > 0.d0) exit
        if (w(i) - w(i-1) > -0.3d0) then
           call ERROR_("x-like block or invalid orientation in normal_ray_wfunc")
        endif
        w(i) = w(i) + 1.d0
     enddo
  enddo
  wfunc = pchip(linspace(0.d0, 1.d0, m+1), w)

  contains
  !.............................................................................
  subroutine ERROR_(msg)
  character(len=*), intent(in) :: msg


  print *, "i = ", i
  print *, "m = ", m
  call C%savetxt("ERROR_NORMAL_RAY_WFUNC_C")
  call P%savetxt("ERROR_NORMAL_RAY_WFUNC_P")
  call ERROR(msg)

  end subroutine ERROR_
  !.............................................................................
  end function normal_ray_wfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine normal_ray_interpolate(this, i1, i2, C2, m, direction, Q)
  use moose_error
  use moose_analysis
  use moose_geometry
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: i1, i2, m
  class(curve), intent(in   ) :: C2
  integer,      intent(in   ), optional :: direction
  class(qfunc), intent(in   ), optional :: Q

  type(polygon2d)    :: P1, P2
  type(interp_curve) :: C1
  type(interp) :: wfunc
  real(real64) :: s
  integer      :: idir, j, n


  ! initial contour at i1
  P1 = polygon2d_qmesh_contour(this, i1, 1)
  C1 = interp_curve(linspace(0.d0, 1.d0, P1%nnodes()), P1%nodes())
  ! final contour at i2
  P2 = polygon2d(C2%polygon(256))


  ! construct mapping function and interpolate nodes in blocks
  if (present(direction)) then
     idir = sign(1,direction)
  else
     idir = P1%orientation()
     if (C1%is_closed) then
        if (P1%winding_number(P2%node(0)) /= 0) idir = -idir
     else
        call ERROR("direction must be given for open mesh contours")
     endif
  endif
  n = this%n(2)
  wfunc = normal_ray_wfunc(C1, m, P2, idir)
  do j=0,n-1
     s = wfunc%eval(j,n-1);   if (C1%is_closed) s = mod(s,1.d0)
     this%x(i2,j,:) = C2%eval(C2%a + (C2%b-C2%a) * s)
  enddo
  call interpolate_submesh(this, i1, i2, 1, Q)

  end subroutine normal_ray_interpolate
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine footpoints_interpolate(this, i1, i2, C2, m, Q)
  use moose_error
  use moose_analysis
  use moose_geometry
  class(qmesh), intent(inout) :: this
  integer,      intent(in   ) :: i1, i2, m
  class(curve), intent(in   ) :: C2
  class(qfunc), intent(in   ), optional :: Q

  type(polygon2d)     :: P1
  type(bspline_curve) :: C1
  type(interp)        :: wfunc
  real(real64), allocatable :: x(:,:), P(:,:), t(:), e(:)
  real(real64), pointer :: xtmp(:,:)
  integer :: ierr, j, n


  n = this%n(2)


  ! 1. reference points on initial mesh contour
  ! all nodes
  if (m == -1) then
     allocate (x(2,0:n-1), P(2,0:n-1), t(0:n-1), e(0:n-1))
     xtmp => qmesh_contour(this, i1, 1);   x = transpose(xtmp)

  ! m blocks
  else
     P1 = polygon2d_qmesh_contour(this, i1, 1)
     C1 = bspline_polygon(P1)

     allocate (x(2,0:m), P(2,0:m), t(0:m), e(0:m))
     do j=0,m
        x(:,j) = C1%eval(j,m)
     enddo
  endif


  ! 2. find footpoints on C2
  call C2%find_footpoints(size(x,2), x, t, P, e, ierr)
  if (ierr /= 0) call ERROR("footpoint construction failed for mesh nodes")


  ! 3. set mesh contour at i2
  if (m == -1) then
     do j=0,n-1
        this%x(i2,j,:) = P(:,j)
     enddo

  else
     wfunc = make_wfunc(t, C2%a, C2%b)
     do j=0,n-1
        this%x(i2,j,:) = C2%eval(C2%inbounds(wfunc%eval(j,n-1)))
     enddo
  endif


  ! 4. interpolate mesh between i1 and i2
  call interpolate_submesh(this, i1, i2, 1, Q)

  end subroutine footpoints_interpolate
  !-----------------------------------------------------------------------------


! auxiliary procedures for mesh contours:
  !-----------------------------------------------------------------------------
  function qmesh_contour(this, i, mdim) result(x)
  !
  ! mesh contour at index *i* along mesh dimensions *mdim* as array pointer
  !
  class(qmesh),    intent(in) :: this
  integer,         intent(in) :: i, mdim
  real(real64), pointer       :: x(:,:)


  nullify(x)
  if (mdim < 1  .or.  mdim > 2) return


  if (mdim == 1) then
     x => this%x(i,:,:)
  elseif (mdim == 2) then
     x => this%x(:,i,:)
  endif

  end function qmesh_contour
  !-----------------------------------------------------------------------------
  function polygon2d_qmesh_contour(this, i, mdim) result(P)
  !
  ! polygon2d representation of qmesh_contour
  !
  use moose_error
  class(qmesh),    intent(in) :: this
  integer,         intent(in) :: i, mdim
  type(polygon2d)             :: P

  real(real64), pointer :: xtmp(:,:)


  if (mdim < 1  .or.  mdim > 2) call ERROR("mdim = 1 or 2 required", "polygon2d_qmesh_contour")
  ! P = polygon2d(transpose(qmesh_contour(this,i,mdim)))
  xtmp => qmesh_contour(this,i,mdim);   P = polygon2d(transpose(xtmp))

  end function polygon2d_qmesh_contour
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_qmesh_contour_array(this, i, mdim, x)
  !
  ! set mesh contour at index *i* along mesh dimension *mdim* from array *x*
  !
  use moose_error
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i, mdim
  real(real64),    intent(in   ) :: x(this%n(3-mdim),2)

  real(real64), pointer :: xtmp(:,:)


  if (mdim < 1  .or.  mdim > 2) call ERROR("mdim = 1 or 2 required", "set_qmesh_contour")
  !qmesh_contour(this, i, mdim) = x
  xtmp => qmesh_contour(this, i, mdim);   xtmp = x

  end subroutine set_qmesh_contour_array
  !-----------------------------------------------------------------------------
  subroutine set_qmesh_contour_polygon2d(this, i, mdim, P)
  !
  ! set mesh contour from polygon *P*
  !
  use moose_error
  class(qmesh),    intent(inout) :: this
  integer,         intent(in   ) :: i, mdim
  type(polygon2d)                :: P


  if (mdim < 1  .or.  mdim > 2) call ERROR("idim = 1 or 2 required", "set_qmesh_contour")
  if (this%n(3-mdim) /= P%nnodes()) call ERROR("number of nodes in polygon is incompatible with mesh")
  call set_qmesh_contour(this, i, mdim, P%nodes())

  end subroutine set_qmesh_contour_polygon2d
  !-----------------------------------------------------------------------------
  subroutine set_qmesh_contour_curve(this, i, mdim, C, Q, parametrization)
  !
  ! set mesh contour from curve *C*
  !
  ! optional:
  !   parametrization    can be set to "arclength" for arc length based node spacings
  !   Q                  quantile function for sampling nodes along *C*
  !
  use moose_error
  use moose_quantiles
  class(qmesh),     intent(inout) :: this
  integer,          intent(in   ) :: i, mdim
  class(curve),     intent(in   ) :: C
  class(qfunc),     intent(in   ), optional :: Q
  character(len=*), intent(in   ), optional :: parametrization

  character(len=128) :: P
  integer :: n


  if (mdim < 1  .or.  mdim > 2) call ERROR("idim = 1 or 2 required", "set_qmesh_contour")
  n = this%n(3-mdim)


  P = "intrinsic";   if (present(parametrization)) P = parametrization
  select case(P)
  case("intrinsic")
     call set_qmesh_contour(this, i, mdim, transpose(C%discretization(n, Q)))

  case("arclength")
     call set_qmesh_contour(this, i, mdim, transpose(C%arclength_discretization(n, Q)))

  case default
     call ERROR("invalid parametrization '"//trim(P)//"'")
  end select

  end subroutine set_qmesh_contour_curve
  !-----------------------------------------------------------------------------


! high-level qmesh constructors:
  !-----------------------------------------------------------------------------
  function quasi_orthogonal_qmesh(C1, C2, nlevels, m, idir, Qlevels, Q) result(this)
  use moose_error
  use moose_math
  use moose_fieldline
  class(curve), intent(in) :: C1, C2
  integer,      intent(in) :: nlevels, m, idir
  class(qfunc), intent(in), optional :: Qlevels, Q
  type(qmesh)              :: this


  this = qmesh(nlevels, m)
  ! start from C1
  if (idir == -1) then
     call set_qmesh_contour(this, 0, 1, C1, Q, "arclength")
     call construct_submesh(this, 0, nlevels-1, C2, Qlevels, 16)

  ! start from C2
  elseif (idir == 1) then
     call set_qmesh_contour(this, nlevels-1, 1, C2, Q, "arclength")
     call construct_submesh(this, nlevels-1, 0, C1, Qlevels, 16)

  else
     call ERROR("invalid dir in quasi_orthogonal_qmesh")
  endif

  end function quasi_orthogonal_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function nray_blocks_qmesh(C1, C2, nu, nv, nblocks, Qn, Qm) result(this)
  class(curve), intent(in) :: C1, C2
  integer,      intent(in) :: nu, nv, nblocks
  class(qfunc), intent(in), optional :: Qn, Qm
  type(qmesh)              :: this


  this = qmesh(nu, nv)
  call set_qmesh_contour(this, 0, 1, C1, Qm, "arclength")
  call normal_ray_interpolate(this, 0, nu-1, C2, nblocks, Q=Qn)

  end function nray_blocks_qmesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function footpoints_qmesh(C1, C2, nu, nv, nblocks, Qu, Qv) result(this)
  class(curve), intent(in) :: C1, C2
  integer,      intent(in) :: nu, nv, nblocks
  class(qfunc), intent(in), optional :: Qu, Qv
  type(qmesh)              :: this


  this = qmesh(nu, nv)
  call set_qmesh_contour(this, 0, 1, C1, Qv, "arclength")
  call footpoints_interpolate(this, 0, nu-1, C2, nblocks, Qu)

  end function footpoints_qmesh
  !-----------------------------------------------------------------------------

end module moose_qmesh_generator
