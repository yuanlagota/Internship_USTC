!===============================================================================
! Multivariate, real-valued scalar and vector functions
!===============================================================================
module moose_mfunc
  use iso_fortran_env
  implicit none
  private


  ! base class for (multivariate) functions ....................................
  ! This is rather the definition of the domain than the map itself.
  ! Definition of the interface to perform the mapping is deferred for
  ! scalar and vector implementations.
  type, abstract :: mfunc
     ! dimension of domain and codomain
     integer :: ndim, mdim

     ! bounding box for domain of definition
     real(real64), allocatable :: lb(:), ub(:)
     ! flag for periodic boundaries
     logical, allocatable :: periodic(:)

     contains
     ! broadcast mfunc to all mpi processes
     procedure :: broadcast
     procedure :: mfunc_broadcast => broadcast

     ! finalize mfunc
     procedure :: free
     procedure :: mfunc_free => free

     ! check if point is outside of bounding box
     procedure :: strict_out_of_bounds
     procedure :: out_of_bounds, mfunc_out_of_bounds => out_of_bounds
     procedure :: bounded_domain
  end type mfunc
  ! type mfunc .................................................................



  ! bass class for scalar functions ............................................
  type, extends(mfunc), abstract, public :: scalar_mfunc
     contains
     ! return function value
     procedure(eval),    deferred :: eval

     ! return first-order partial derivatives
     procedure(deriv),   deferred :: deriv

     ! compute both function value and first-order derivatives
     ! NOTE: this procedure is provided for convenience; implementations of
     ! scalar_mfunc may overload this procedure for speedup
     procedure :: fdf

     ! return Hessian matrix
     procedure(hessian), deferred :: hessian
  end type scalar_mfunc


  abstract interface
     ! return function value at x
     function eval(this, x) result(v)
     use iso_fortran_env
     import scalar_mfunc
     class(scalar_mfunc), intent(in) :: this
     real(real64),        intent(in) :: x(this%ndim)
     real(real64)                    :: v
     end function eval


     ! return first-order partial derivatives at x
     function deriv(this, x)
     use iso_fortran_env
     import scalar_mfunc
     class(scalar_mfunc), intent(in) :: this
     real(real64),        intent(in) :: x(this%ndim)
     real(real64)                    :: deriv(this%ndim)
     end function deriv


     ! return Hessian matrix at x
     function hessian(this, x)
     use iso_fortran_env
     import scalar_mfunc
     class(scalar_mfunc), intent(in) :: this
     real(real64),        intent(in) :: x(this%ndim)
     real(real64)                    :: hessian(this%ndim, this%ndim)
     end function hessian
  end interface
  ! type scalar_mfunc ..........................................................



  ! bass class for vector functions ............................................
  type, extends(mfunc), abstract, public :: vector_mfunc
     contains
     ! return function value
     procedure(vector_eval), deferred :: eval

     ! return Jacobian
     procedure(jac), deferred :: jac

     ! sample vector function on domain
     procedure :: sample => sample_vector_mfunc
  end type vector_mfunc


  abstract interface
     ! return function value at x
     function vector_eval(this, x) result(v)
     use iso_fortran_env
     import vector_mfunc
     class(vector_mfunc), intent(in) :: this
     real(real64),        intent(in) :: x(this%ndim)
     real(real64)                    :: v(this%mdim)
     end function vector_eval


     ! return Jacobian at x
     function jac(this, x)
     use iso_fortran_env
     import vector_mfunc
     class(vector_mfunc), intent(in) :: this
     real(real64),        intent(in) :: x(this%ndim)
     real(real64)                    :: jac(this%mdim, this%ndim)
     end function jac
  end interface
  ! type vector_mfunc ..........................................................



  ! extended base classes for 2D and 3D functions ..............................
  type, extends(scalar_mfunc), abstract, public :: scalar_mfunc2d
  end type scalar_mfunc2d
  type, extends(scalar_mfunc), abstract, public :: scalar_mfunc3d
  end type scalar_mfunc3d
  type, extends(vector_mfunc), abstract, public :: vector_mfunc2d
  end type vector_mfunc2d
  type, extends(vector_mfunc), abstract, public :: vector_mfunc3d
  end type vector_mfunc3d
  !.............................................................................



  interface sample_mfunc2d
     procedure :: sample_mfunc2d_mesh
  end interface sample_mfunc2d



  public :: &
     init_mfunc, &
     init_mfunc2d, &
     init_mfunc3d, &
     init_scalar_mfunc2d, &
     init_scalar_mfunc3d, &
     assert_mfunc2d, &
     sample_mfunc2d, aux_sample_mfunc2d, &
     find_nearest_critical_point, &
     find_critical_points, &
     find_critical_points2d, &
     find_critical_points3d, &
     line_search, root_finder


  contains
  !-----------------------------------------------------------------------------


! class mfunc ==================================================================
! constructor procedures:
  !-----------------------------------------------------------------------------
  subroutine init_mfunc(this, ndim, mdim, lb, ub, periodic)
  class(mfunc), intent(inout) :: this
  integer,      intent(in)    :: ndim, mdim
  real(real64), intent(in), optional :: lb(ndim), ub(ndim)
  logical,      intent(in), optional :: periodic(ndim)


  ! set dimension of domain and codomain
  this%ndim = ndim
  this%mdim = mdim


  ! set lower boundaries
  if (present(lb)) then
     allocate (this%lb(ndim), source=lb)
  else
     allocate (this%lb(ndim), source=-huge(1.d0))
  endif


  ! set upper boundaries
  if (present(ub)) then
     allocate (this%ub(ndim), source=ub)
  else
     allocate (this%ub(ndim), source=huge(1.d0))
  endif


  ! set boundary types
  if (present(periodic)) then
     allocate (this%periodic(ndim), source=periodic)
  else
     allocate (this%periodic(ndim), source=.false.)
  endif

  end subroutine init_mfunc
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_mfunc2d(this, mdim, lb, ub, periodic)
  class(mfunc), intent(inout) :: this
  integer,      intent(in)    :: mdim
  real(real64), intent(in), optional :: lb(2), ub(2)
  logical,      intent(in), optional :: periodic(2)


  call init_mfunc(this, 2, mdim, lb, ub, periodic)

  end subroutine init_mfunc2d
  !-----------------------------------------------------------------------------
  subroutine init_scalar_mfunc2d(this, lb, ub, periodic)
  class(scalar_mfunc2d), intent(inout) :: this
  real(real64),          intent(in), optional :: lb(2), ub(2)
  logical,               intent(in), optional :: periodic(2)


  call init_mfunc(this, 2, 1, lb, ub, periodic)

  end subroutine init_scalar_mfunc2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_mfunc3d(this, mdim, lb, ub, periodic)
  class(mfunc), intent(inout) :: this
  integer,      intent(in)    :: mdim
  real(real64), intent(in), optional :: lb(3), ub(3)
  logical,      intent(in), optional :: periodic(3)


  call init_mfunc(this, 3, mdim, lb, ub, periodic)

  end subroutine init_mfunc3d
  !-----------------------------------------------------------------------------
  subroutine init_scalar_mfunc3d(this, lb, ub, periodic)
  class(scalar_mfunc3d), intent(inout) :: this
  real(real64),          intent(in), optional :: lb(3), ub(3)
  logical,               intent(in), optional :: periodic(3)


  call init_mfunc(this, 3, 1, lb, ub, periodic)

  end subroutine init_scalar_mfunc3d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(mfunc), intent(inout) :: this


  call proc(0)%broadcast(this%ndim)
  call proc(0)%broadcast(this%mdim)
  call proc(0)%broadcast_allocatable(this%lb)
  call proc(0)%broadcast_allocatable(this%ub)
  call proc(0)%broadcast_allocatable(this%periodic)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(mfunc), intent(inout) :: this


  deallocate(this%lb, this%ub, this%periodic)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function strict_out_of_bounds(this, x) result(out_of_bounds)
  !
  ! return true if point is outside of bounding box, including periodic dimensions
  !
  class(mfunc), intent(in) :: this
  real(real64), intent(in) :: x(this%ndim)
  logical                  :: out_of_bounds

  integer :: i


  out_of_bounds = .false.
  do i=1,this%ndim
     if (x(i) < this%lb(i)  .or.  x(i) > this%ub(i)) then
        out_of_bounds = .true.
        return
     endif
  enddo

  end function strict_out_of_bounds
  !-----------------------------------------------------------------------------
  function out_of_bounds(this, x)
  !
  ! return true if point is outside of bounding box, excluding periodic dimensions
  !
  class(mfunc), intent(in) :: this
  real(real64), intent(in) :: x(this%ndim)
  logical                  :: out_of_bounds

  integer :: i


  out_of_bounds = .false.
  do i=1,this%ndim
     if (this%periodic(i)) cycle

     if (x(i) < this%lb(i)  .or.  x(i) > this%ub(i)) then
        out_of_bounds = .true.
        return
     endif
  enddo

  end function out_of_bounds
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function bounded_domain(this)
  class(mfunc), intent(in) :: this
  logical                  :: bounded_domain


  bounded_domain = .true.
  if (any(this%ub == huge(1.d0))  .or.  any(this%lb == -huge(1.d0))) then
     bounded_domain = .false.
  endif

  end function bounded_domain
  !-----------------------------------------------------------------------------
! class mfunc ==================================================================



! class scalar_mfunc ===========================================================
! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine fdf(this, x, f, df)
  class(scalar_mfunc), intent(in)  :: this
  real(real64),        intent(in)  :: x(this%ndim)
  real(real64),        intent(out) :: f, df(this%ndim)


  f  = this%eval(x)
  df = this%deriv(x)

  end subroutine fdf
  !-----------------------------------------------------------------------------
! class scalar_mfunc ===========================================================



! class vector_mfunc ===========================================================
  !-----------------------------------------------------------------------------
  subroutine sample_vector_mfunc(this, basename, n, lb, ub)
  use moose_dataset
  use moose_grids
  class(vector_mfunc), intent(in) :: this
  character(len=*),    intent(in) :: basename
  integer,             intent(in), optional :: n(this%ndim)
  real(real64),        intent(in), optional :: lb(this%ndim), ub(this%ndim)

  character(len=len(basename)+5)  :: filename
  class(grid), allocatable :: G
  real(real64), pointer :: element(:)
  type(dataset) :: D
  real(real64)  :: x(this%ndim), x1(this%ndim), x2(this%ndim), v(this%mdim)
  integer       :: i, nn(this%ndim)


  ! set domain
  x1 = this%lb;   if (present(lb)) x1 = lb
  x2 = this%ub;   if (present(ub)) x2 = ub
  nn = 128;       if (present(n))  nn = n
  select case(this%ndim)
  case(2)
     allocate (G, source=rmesh(x1(1), x2(1), nn(1), x1(2), x2(2), nn(2)))

  case(3)
     allocate (G, source=cmesh(nn(1), nn(2), nn(3), x1(1), x2(2), x1(2), x2(2), x1(3), x2(3)))

  case default
     return
  end select
  write (filename, 1001) basename
 1001 format(a,".grid")
  call G%savetxt(filename)


  ! evaluate vector field on grid
  D = dataset(this%mdim, G%nnodes(), geometry=filename)
  do i=0,G%nnodes()-1
     x = G%node(i)
     v = this%eval(x)
     element => D%element(i)
     element = v
  enddo
  write (filename, 2001) basename
 2001 format(a,".dat")
  call D%savetxt(filename)

  end subroutine sample_vector_mfunc
  !-----------------------------------------------------------------------------
! class vector_mfunc ===========================================================



! module procedures:
! sample =======================================================================
  !---------------------------------------------------------------------
  subroutine assert_mfunc2d(F, procname)
  use moose_error
  class(mfunc),     intent(in) :: F
  character(len=*), intent(in) :: procname


  if (F%ndim /= 2) then
     call ERROR("mfunc with ndim == 2 required in procedure '"//procname//"'")
  endif

  end subroutine assert_mfunc2d
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine sample_mfunc2d_mesh(F, basename, nx, ny, x1, x2, y1, y2)
  use moose_grids, only: rmesh
  class(scalar_mfunc),           intent(in) :: F
  character(len=*),              intent(in) :: basename
  integer,                       intent(in) :: nx, ny
  real(real64),        optional, intent(in) :: x1, x2, y1, y2

  type(rmesh)   :: M
  real(real64)  :: lb(2), ub(2)


  call assert_mfunc2d(F, "sample_mfunc2d")


  ! set sample domain
  lb = F%lb
  ub = F%ub
  if (present(x1)) lb(1) = x1
  if (present(x2)) ub(1) = x2
  if (present(y1)) lb(2) = y1
  if (present(y2)) ub(2) = y2
  ! check domain boundaries
  if (any(lb == -huge(1.d0))  .or.  any(ub == huge(1.d0))) then
     write (6, *) "ERROR: finite domain boundaries must be set in sample_field2D"
     stop
  endif

  M = rmesh(lb(1), ub(1), nx, lb(2), ub(2), ny)
  call aux_sample_mfunc2d(F, basename, M)

  end subroutine sample_mfunc2d_mesh
  !.............................................................................
  subroutine aux_sample_mfunc2d(F, basename, G)
  use moose_grids,   only: grid
  use moose_dataset, only: dataset, metadata
  class(scalar_mfunc), intent(in) :: F
  character(len=*),    intent(in) :: basename
  class(grid),         intent(in) :: G


  character(len=len(basename)+5) :: filename
  real(real64), pointer :: element(:)
  type(dataset) :: D
  real(real64)  :: x(2), Fval, dF1(2), Hessian(2,2)
  integer       :: i


  write (filename, 1001) basename
 1001 format(a,".grid")
  call G%savetxt(filename)


  D = dataset(6, G%nnodes(), geometry=filename)
  call D%set_metadata(1, "value", "Field value")
  call D%set_metadata(2, "d1",    "Derivative in first direction")
  call D%set_metadata(3, "d2",    "Derivative in second direction")
  call D%set_metadata(4, "H11",   "1-1 element of Hessian matrix")
  call D%set_metadata(5, "H12",   "1-2 element of Hessian matrix")
  call D%set_metadata(6, "H22",   "2-2 element of Hessian matrix")
  do i=0,G%nnodes()-1
     x = G%node(i)

     Fval    = F%eval(x)
     dF1     = F%deriv(x)
     Hessian = F%hessian(x)

     element => D%element(i)
     element = [Fval, dF1(1), dF1(2), Hessian(1,1), Hessian(1,2), Hessian(2,2)]
  enddo
  write (filename, 2001) basename
 2001 format(a,".dat")
  call D%savetxt(filename)

  end subroutine aux_sample_mfunc2d
  !---------------------------------------------------------------------
! sample =======================================================================



! critical points ==============================================================
  !-----------------------------------------------------------------------------
  subroutine find_nearest_critical_point(this, x, x0, ierr, method)
  !
  ! find nearest critical point of scalar function
  !
  class(scalar_mfunc), intent(in)  :: this
  real(real64),        intent(in)  :: x(this%ndim)
  real(real64),        intent(out) :: x0(this%ndim)
  integer,             intent(out) :: ierr
  integer,             intent(in), optional :: method


  call newton_method(this, x, 8.d-1, 1.d-10, 100, x0, ierr)

  end subroutine find_nearest_critical_point
  !-----------------------------------------------------------------------------
  subroutine newton_method(this, x, lambda, delta, nmax, x0, ierr)
  use moose_linalg, only: inverse => inverse_2d ! @todo: implementation for ndim > 2
  class(scalar_mfunc), intent(in)  :: this
  real(real64),        intent(in)  :: x(this%ndim), lambda, delta
  integer,             intent(in)  :: nmax
  real(real64),        intent(out) :: x0(this%ndim)
  integer,             intent(out) :: ierr

  real(real64) :: dx(this%ndim), dxmod
  real(real64) :: df(this%ndim), H(this%ndim, this%ndim), H1(this%ndim, this%ndim)
  integer :: i, j


  ierr = 0
  x0   = x
  do i=1,nmax
     ! check boundaries of domain
     if (this%out_of_bounds(x0)) then
        ierr = 1
        return
     endif

     ! calculate Gradient and Hessian
     df = this%deriv(x0)
     H  = this%Hessian(x0)
     call inverse(H, H1, ierr)
     if (ierr /= 0) then
        ierr = 2
        return
     endif

     ! calculate and apply step
     do j=1,this%ndim
        dx(j) = lambda * sum(H1(j,:) * df)
     enddo
     x0 = x0 - dx

     ! check convergence
     dxmod = sqrt(sum(dx**2))
     if (dxmod < delta) return
  enddo
  ierr = 3

  end subroutine newton_method
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function sample_critical_points(Phi, G) result(X)
  !
  ! return list of critical points for scalar multi-variate function Phi
  ! G: grid with sample points
  !
  use moose_rlist
  use moose_grids
  class(scalar_mfunc), intent(in) :: Phi
  class(grid),         intent(in) :: G
  type(rlist)                     :: X

  real(real64), parameter :: min2 = 1.d-12

  real(real64) :: y(Phi%ndim), y0(Phi%ndim), d
  integer :: i, istat, j


  X = rlist(Phi%ndim)
  if (G%ndim /= Phi%ndim) return

  grid_loop: do i=0,G%nnodes()-1
     y = G%node(i)
     if (Phi%out_of_bounds(y)) cycle

     ! find closest critical point to grid node
     call find_nearest_critical_point(Phi, y, y0, istat)
     if (istat /= 0) cycle


     ! check if this critical point is already known
     do j=0,X%nelements()-1
        d = sum((X%element(j) - y0)**2)
        if (d < min2) cycle grid_loop
     enddo


     ! save this critical point
     call X%append(y0)
  enddo grid_loop

  end function sample_critical_points
  !-----------------------------------------------------------------------------
  function find_critical_points(Phi, nsample, xrange) result(X)
  !
  ! return list of critical points for scalar multi-variate function Phi
  !
  use moose_rlist
  use moose_grids
  class(scalar_mfunc), intent(in) :: Phi
  integer,             intent(in), optional :: nsample(Phi%ndim)
  real(real64),        intent(in), optional :: xrange(2,Phi%ndim)
  type(rlist)                     :: X

  class(grid), allocatable :: G
  real(real64) :: lb(Phi%ndim), ub(Phi%ndim)
  integer :: n(Phi%ndim)


  ! set bounding box for scan domain
  lb = Phi%lb;   if (present(xrange)) lb = xrange(1,:)
  ub = Phi%ub;   if (present(xrange)) ub = xrange(2,:)


  n = 32;   if (present(nsample)) n = nsample
  select case(Phi%ndim)
  case(2)
     allocate (G, source=rmesh(lb(1), ub(1), n(1), lb(2), ub(2), n(2)))

  case(3)
     allocate (G, source=cmesh(n(1), n(2), n(3), lb(1), ub(1), lb(2), ub(2), lb(3), ub(3)))

  case default
     X = rlist(Phi%ndim)
     return
  end select
  X = sample_critical_points(Phi, G)


  end function find_critical_points
  !-----------------------------------------------------------------------------
  function find_critical_points2d(Phi, nsample, xrange) result(X)
  !
  ! wrapper for return type rlist2
  !
  use moose_rlist
  class(scalar_mfunc), intent(in) :: Phi
  integer,             intent(in), optional :: nsample(Phi%ndim)
  real(real64),        intent(in), optional :: xrange(2,Phi%ndim)
  type(rlist2)                    :: X


  if (Phi%ndim /= 2) then
     write (6, 9000) Phi%ndim;   stop
  endif
 9000 format("ERROR: find_critical_points2d called with Phi%ndim = ",i0," /= 2!")

  X%rlist = find_critical_points(Phi, nsample, xrange)

  end function find_critical_points2d
  !-----------------------------------------------------------------------------
  function find_critical_points3d(Phi, nsample, xrange) result(X)
  !
  ! wrapper for return type rlist3
  !
  use moose_rlist
  class(scalar_mfunc), intent(in) :: Phi
  integer,             intent(in), optional :: nsample(Phi%ndim)
  real(real64),        intent(in), optional :: xrange(2,Phi%ndim)
  type(rlist3)                    :: X


  if (Phi%ndim /= 3) then
     write (6, 9000) Phi%ndim;   stop
  endif
 9000 format("ERROR: find_critical_points2d called with Phi%ndim = ",i0," /= 3!")

  X%rlist = find_critical_points(Phi, nsample, xrange)

  end function find_critical_points3d
  !-----------------------------------------------------------------------------
! critical points ==============================================================



! root-finding =================================================================
  !-----------------------------------------------------------------------------
  function line_search(Phi, xa, xb, T, istat, nmax, delta) result(x0)
  !
  ! return x0 on line between xa and xb with Phi(x0) = T, or istat > 0
  !
  ! algorithm: bisection method
  !
  class(scalar_mfunc), intent(in)  :: Phi
  real(real64),        intent(in)  :: xa(Phi%ndim), xb(Phi%ndim), T
  integer,             intent(out) :: istat
  integer,             intent(in), optional :: nmax
  real(real64),        intent(in), optional :: delta
  real(real64)                     :: x0(Phi%ndim)

  real(real64), parameter :: default_delta = 1.d-10
  integer,      parameter :: default_nmax  = 256

  real(real64) :: Phia, Phib, Phic, d, a(Phi%ndim), b(Phi%ndim), c(Phi%ndim)
  integer      :: i, n


  istat = 0
  a     = xa
  b     = xb
  Phia  = Phi%eval(a)
  Phib  = Phi%eval(b)
  if ((Phia-T) * (Phib-T) >= 0) then
     istat = 1
     return
  endif


  n = default_nmax;          if (present(nmax))  n = nmax
  d = 4*default_delta**2 ;   if (present(delta)) d = 4*delta**2
  do i=1,n
     ! new midpoint
     c    = (a + b) / 2
     Phic = Phi%eval(c)

     ! (good enough) solution found?
     if (Phic == T  .or.  sum((a-b)**2) < d) then
        x0 = c
        return
     endif

     ! select new subinterval
     if ((Phia-T) * (Phic-T) > 0.d0) then
        a    = c
        Phia = Phic
     else
        b    = c
        Phib = Phic
     endif
  enddo
  istat = 2

  end function line_search
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine root_finder(this, x0, T, v, istat, epsabs, damping, periodic)
  !
  ! find nearest location at target value in direction v
  ! Input:
  !    x0        start/reference location
  !    T         target value
  !    v         search direction
  !    periodic  treat boundaries as periodic (default = .false.)
  !
  ! Output:
  !    x0      location with this%eval(x0) = T
  !    istat   = 0: success
  !              1: max. # iterations reached without convergence
  !              2: df/ds = 0 encountered along v
  !              3: start location out of bounds
  !              4: at boundary
  !
  class(scalar_mfunc), intent(in)    :: this
  real(real64),        intent(inout) :: x0(this%ndim)
  real(real64),        intent(in)    :: T, v(this%ndim)
  integer,             intent(out)   :: istat
  real(real64),        intent(in), optional :: epsabs, damping
  logical,             intent(in), optional :: periodic(this%ndim)


  integer, parameter :: nmax = 256

  logical      :: lp(this%ndim)
  real(real64) :: a, eps, ds, x(this%ndim), xl(this%ndim), f, dfds, fT
  integer :: i


  ! set defaults for numerical parameters (optional input)
  lp  = .false.;   if (present(periodic)) lp  = periodic  ! periodic boundaries
  eps = 1.d-12;    if (present(epsabs))   eps = epsabs    ! required accuracy of root
  a   = 0.9d0;     if (present(damping))  a   = damping


  ! start location out of bounds
  if (this%out_of_bounds(x0)) then
     istat = 3
     return
  endif


  istat = 1
  x     = x0
  xl    = x
  fT    = T
  do i=1,nmax
     ! evaluate f and df/ds at x
     f    = this%eval(x) - fT
     dfds = sum(this%deriv(x) * v)


     ! quit if df/ds = 0
     if (dfds == 0.d0) then
        istat = 2
        return
     endif


     ! calculate increment
     ds = -f / dfds * a
     x  = xl + ds * v


     ! check if location is out of bounds
     if (this%out_of_bounds(x)) then
        call boundary_point(istat)
        if (istat == 4) return
        cycle
     endif


     ! convergence check
     if (abs(ds) < eps) then
        x0    = x
        istat = 0
        return
     endif


     ! prepare next step
     xl = x
  enddo

  contains
  !.............................................................................
  subroutine boundary_point(istat)
  integer, intent(out) :: istat

  real(real64) :: xm(this%ndim), xmjb, f, fm
  integer :: j, jb


  ! calculate point on boundary
  do j=1,this%ndim
     ! at upper boundary
     if (this%ub(j) - xl(j)  <  ds * v(j)) then
        ds    = (this%ub(j) - xl(j)) / v(j)

        ! mirror point for periodic boundaries
        jb    = j
        xmjb  = this%lb(j)
     endif

     ! at lower boundary
     if (xl(j) - this%lb(j) <  -ds * v(j)) then
        ds    = (this%lb(j) - xl(j)) / v(j)

        ! mirror point for periodic boundaries
        jb    = j
        xmjb  = this%ub(j)
     endif
  enddo
  x0    = xl + ds * v
  istat = 4

  ! should this be treated as periodic boundary?
  if (lp(jb)) then
     ! evaluate function on boundary
     f      = this%eval(x0)

     ! evaluate function at mirror point on boundary
     xm     = x0
     xm(jb) = xmjb
     fm     = this%eval(xm)


     ! set mirror point as new point, and update target value
     x      = xm
     fT     = fT + fm - f

     ! continue from here
     xl     = x
     istat  = 1
     !write (6, 4002) x0(2), xm(2), fT
 4002 format("reached boundary at ",f10.5," continuing at ",f10.5," with new target ",f10.5)
  endif

  end subroutine boundary_point
  !.............................................................................
  end subroutine root_finder
  !-----------------------------------------------------------------------------
! root-finding =================================================================


end module moose_mfunc
