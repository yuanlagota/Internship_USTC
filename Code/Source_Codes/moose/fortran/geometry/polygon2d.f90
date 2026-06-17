module moose_polygon2d
  use iso_fortran_env
  use moose_polygon
  implicit none
  private


  integer, parameter, public :: &
     XSECT_SEGMENT = 1, &
     XSECT_RAY     = 2



  ! Polygon in 2D ..............................................................
  type, extends(polygon), public :: polygon2d
     ! bounding boxes for polygon segments
     real(real64), dimension(:,:), allocatable :: rbox, zbox

     contains
     procedure :: set_node
     procedure :: normal
     procedure :: secant_normal

     procedure :: aux_distance, get_distance
     procedure :: orientation
     procedure :: winding_number => polygon2d_winding_number
     generic   :: intersect => xsect, intersect_polygon2d
     procedure :: xsect, intersect_polygon2d
     generic   :: intersects => intersects_segment, intersects_polygon2d
     procedure :: intersects_segment, intersects_polygon2d
     procedure :: area
     procedure :: shift

     procedure :: broadcast
     procedure :: free
     procedure :: update, reverse
  end type polygon2d


  interface polygon2d
     procedure :: new
     procedure :: convert_polygon
     procedure :: construct
     procedure :: construct_from_list
     procedure :: load
  end interface
  ! polygon2d ..................................................................



  interface sample_mfunc2d
     procedure :: sample_mfunc2d_polygon
  end interface sample_mfunc2d



  public :: &
     shifted_polygon2d, winding_number, sample_mfunc2d


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(n) result(this)
  !
  ! allocate memory for polygon with n segments
  !
  integer, intent(in) :: n
  type(polygon2d)     :: this


  this = convert_polygon(polygon(n, 2))

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function convert_polygon(P) result(this)
  !
  ! conversion of polygon to polygon2d
  !
  use moose_error
  use moose_txtio
  class(polygon), intent(in) :: P
  type(polygon2d)            :: this


  if (P%ndim /= 2) call ERROR("polygon must be 2-D for conversion to polygon2d")

  call init_txtio(this, "polygon2d")
  this%polygon%implementation = P%implementation
  this%ndim = 2
  call aux_init_polygon2d(this)
  call this%update()

  end function convert_polygon
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct(nodes) result(this)
  !
  ! construct 2D polygon from array of nodes
  !
  real(real64), dimension(:,:) :: nodes
  type(polygon2d)              :: this


  if (size(nodes,1) /= 2) then
     write (6, 9000) size(nodes,1);   stop
  endif
 9000 format("error in polygon2d: dimension of nodes is ",i0," /= 2!")


  this = convert_polygon(polygon(nodes))

  end function construct
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function construct_from_list(L) result(this)
  !
  ! construct 2D polygon from list of 2D nodes
  !
  use moose_rlist
  class(rlist2), intent(in) :: L
  type(polygon2d)           :: this


  this = convert_polygon(polygon(L%rlist))

  end function construct_from_list
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename, scale) result(P)
  !
  ! load 2D polygon from text file
  !
  use moose_txtio
  use moose_table
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale
  type(polygon2d)              :: P


  call init_txtio(P, "polygon2d")
  P%ndim = 2
  P%implementation = table(filename, columns=2, transposed=.true., lbounds=[1,0])
  if (present(scale)) P%implementation%values = scale * P%implementation%values
  call P%remove_duplicate_nodes()
  call aux_init_polygon2d(P)
  call P%update()

  end function load
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_polygon2d(P)
  class(polygon2d), intent(inout) :: P

  integer :: n


  n = P%segments()
  allocate (P%rbox(2,n), P%zbox(2,n), source=0.d0)

  end subroutine aux_init_polygon2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function shifted_polygon2d(P, ds) result(Pshifted)
  !
  ! Construct shifted polygon from P.
  !
  class(polygon2d), intent(in) :: P
  real(real64),     intent(in) :: ds
  type(polygon2d)              :: Pshifted

  integer :: ierr


  Pshifted = P
  call Pshifted%shift(ds, ierr)
  if (ierr /= 0) then
     print 9000
     stop
  endif
 9000 format("ERROR in shifted_polygon2d: ierr = ",i0)

  end function shifted_polygon2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  pure subroutine set_node(this, i, x)
  class(polygon2d), intent(inout) :: this
  integer,          intent(in   ) :: i
  real(real64),     intent(in   ) :: x(this%ndim)

  real(real64) :: x1(2), x2(2)


  call this%polygon%set_node(i, x)
  if (i > 0) then
     x1 = this%node(i-1)
     this%rbox(1,i) = min(x1(1), x(1))
     this%rbox(2,i) = max(x1(1), x(1))
     this%zbox(1,i) = min(x1(2), x(2))
     this%zbox(2,i) = max(x1(2), x(2))
  endif
  if (i < this%segments()) then
     x2 = this%node(i+1)
     this%rbox(1,i+1) = min(x(1), x2(1))
     this%rbox(2,i+1) = max(x(1), x2(1))
     this%zbox(1,i+1) = min(x(2), x2(2))
     this%zbox(2,i+1) = max(x(2), x2(2))
  endif

  end subroutine set_node
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function normal(this, i) result(u)
  !
  ! return normal vector for segment i
  !
  class(polygon2d), intent(in) :: this
  integer, intent(in) :: i

  real(real64) :: v(2), u(2)


  v    = this%node(i+1) - this%node(i)
  u(1) = - v(2)
  u(2) =   v(1)
  u    = u / sqrt(sum(u**2))

  end function normal
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function secant_normal(this, i) result(u)
  !
  ! normal vector for secant through nodes i-1 and i+1
  !
  class(polygon2d), intent(in) :: this
  integer,          intent(in) :: i

  real(real64) :: v(2), u(2)
  integer :: i1, i2, n


  n = this%nnodes()
  if (this%is_closed()) then
     i1 = modulo(i-1, n-1)
     i2 = modulo(i+1, n-1)
  else
     i1 = max(0,i-1)
     i2 = min(i+1,n-1)
  endif

  v    = this%node(i2) - this%node(i1)
  u(1) = -v(2)
  u(2) =  v(1)
  u    = u / sqrt(sum(u**2))

  end function secant_normal
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_distance(this, p, k, t, d, vn)
  !
  ! compute distance between *p* and polygon
  !
  ! returns:
  !    k     >= 0 index of closest node on polygon (if point is near convex corner)
  !          < -index-1 of closest segment on polygon
  !
  !    t     relative position of footpoint on segment k
  !    d     distance to polygon
  !    vn    normal vector
  !
  class(polygon2d), intent(in   ) :: this
  real(real64),     intent(in   ) :: p(2)
  integer,          intent(  out) :: k
  real(real64),     intent(  out) :: t, d, vn(2)

  real(real64) :: x1(2), x2(2), ex(2), en(2), dx, di
  logical :: closed
  integer :: i, i1, i2, ip1, im1
  associate (x => this%implementation%values)


  ! 0. initialize
  d = huge(1.d0)
  closed = this%is_closed()


  ! 1. minimum distance to nodes
  i1 = lbound(x,2)
  i2 = ubound(x,2)
  do i=i1,i2
     x1 = x(:,i)
     di = norm2(x1-p)

     ! tangential and normal vectors
     im1 = i-1
     ip1 = i+1
     if (closed) then
        if (im1 < i1) im1 = i2-1
        if (ip1 > i2) ip1 = i1+1
     else
        im1 = max(im1,i1);   ip1 = min(ip1,i2)
     endif
     ex    = x(:,ip1) - x(:,im1);   ex = ex / norm2(ex)
     en(1) = ex(2);   en(2) = -ex(1)
     if (sum(en * (p-x1)) < 0.d0) di = -di

     if (abs(di) < abs(d)) then
        k = i
        d = di
        vn = en
     endif
  enddo


  ! 2. minimum distance to segments
  do i=lbound(x,2),ubound(x,2)-1
     x1 = x(:,i  )
     x2 = x(:,i+1)
     ex = x2 - x1;   dx = norm2(ex);   ex = ex / dx
     t  = sum(ex * (x2-p)) / dx

     ! footpoint is outside of segment i
     if (t < 0.d0  .or.  t > 1.d0) cycle

     en(1) = ex(2);   en(2) = -ex(1)
     di = sum(en * (p-x1))

     if (abs(di) < abs(d)) then
        k = -i-1
        d = di
        vn = en
     endif
  enddo


  end associate
  end subroutine aux_distance
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_distance(this, p) result(d)
  !
  ! calculate distance of point P to polygon
  !
  class(polygon2d), intent(in) :: this
  real(real64),     intent(in) :: p(2)
  real(real64)                 :: d

  real(real64) :: t, vn(2)
  integer :: k


  call this%aux_distance(p, k, t, d, vn)

  end function get_distance
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function orientation(this)
  !
  ! orientation of polygon: -1 (clockwise) or 1 (counter-clockwise)
  !
  class(polygon2d), intent(in) :: this
  integer                      :: orientation

  real(real64) :: dBA(2), dCA(2), x(2), xA(2), xB(2), xC(2)
  integer :: i, iB, n


  n = this%nnodes()
  if (this%is_closed()) then
     xB = -huge(1.d0)
     do i=0,n-1
        x = this%node(i)
        if (x(1) > xB(1)  .or.  (x(1) == xB(1)  .and.  x(2) > xB(2))) then
           iB = i
           xB = x
        endif
     enddo
     xA = this%node(modulo(iB-1, n-2))
     xC = this%node(modulo(iB+1, n-2))

  else
     xA = this%node(0)
     xC = this%node(n-1)
     xB = 0.d0
     do i=0,n-2
        xB = xB + (this%node(i) + this%node(i+1)) / 2
     enddo
     xB = xB / (n-1)
  endif


  dBA = xB - xA
  dCA = xC - xA
  orientation = 1;   if (dBA(1)*dCA(2) - dCA(1)*dBA(2) < 0.d0) orientation = -1

  end function orientation
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! Point inside polygon test (in 2D)
  ! Implementation of algorithm 6 in
  ! K. Hormann et al. Computational Geometry 20 (2001) 131-144
  !-----------------------------------------------------------------------------
  pure function polygon2d_winding_number(this, p) result(wn)
  class(polygon2d), intent(in) :: this
  real(real64),     intent(in) :: p(2)
  integer                      :: wn


  wn = winding_number(this%implementation%values, p, .false.)

  end function polygon2d_winding_number
  !-----------------------------------------------------------------------------
  pure function winding_number(x, p, wrap) result(wn)
  real(real64),     intent(in) :: x(:,:), p(2)
  logical,          intent(in) :: wrap
  integer                      :: wn

  real(real64) :: det
  integer :: i, i1, i2, ip1, w


  i1 = lbound(x,2)
  i2 = ubound(x,2);   if (.not.wrap) i2 = i2 - 1


  wn = 0
  do i=i1,i2
     ip1 = i + 1;   if (wrap .and. i == i2) ip1 = i1

     ! if crossing
     if (x(2,i) < p(2)  .neqv.  x(2,ip1) < p(2)) then
        if (x(1,i) >= p(1)) then
           if (x(1,ip1) > p(1)) then
              ! modify wn
              wn = wn - 1;   if (x(2,ip1) > x(2,i)) wn = wn + 2
           else
              det = (x(1,i)-p(1))*(x(2,ip1)-p(2)) - (x(1,ip1)-p(1))*(x(2,i)-p(2))
              ! if right crossing
              if (det > 0.d0  .eqv.  x(2,ip1) > x(2,i)) then
                 ! modify wn
                 wn = wn - 1;   if (x(2,ip1) > x(2,i)) wn = wn + 2
              elseif (det == 0.d0) then
                 return
              endif
           endif

        elseif (x(1,ip1) > p(1)) then
           det = (x(1,i)-p(1))*(x(2,ip1)-p(2)) - (x(1,ip1)-p(1))*(x(2,i)-p(2))
           ! if right crossing
           if (det > 0.d0  .eqv.  x(2,ip1) > x(2,i)) then
              ! modify wn
              wn = wn - 1;   if (x(2,ip1) > x(2,i)) wn = wn + 2
           elseif (det == 0.d0) then
              return
           endif
        endif
     endif
  enddo

  end function winding_number
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine update(P)
  !
  ! update bounding boxes for polygon segments
  !
  ! @todo: optional input for selected node
  !
  class(polygon2d), intent(inout) :: P

  real(real64) :: x1(2), x2(2)
  integer :: i, n


  n  = P%segments()
  x1 = P%node(0)
  do i=1,n
     x2 = P%node(i)

     P%rbox(1,i) = min(x1(1), x2(1))
     P%rbox(2,i) = max(x1(1), x2(1))
     P%zbox(1,i) = min(x1(2), x2(2))
     P%zbox(2,i) = max(x1(2), x2(2))

     x1 = x2
  enddo

  end subroutine update
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine reverse(this)
  class(polygon2d), intent(inout) :: this

  integer :: n


  n  = this%segments()
  call this%polygon%reverse()
  this%rbox = this%rbox(:,n:1:-1)
  this%zbox = this%zbox(:,n:1:-1)

  end subroutine reverse
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine xsect(this, p, arg2, xmode, x, r, s, n, count)
  class(polygon2d), intent(in)  :: this
  real(real64),     intent(in)  :: p(2), arg2(2)
  integer,          intent(in)  :: xmode
  real(real64),     intent(out) :: x(2), r, s
  integer,          intent(out) :: n
  integer,          intent(out), optional :: count

  real(real64) :: rbox(2), zbox(2), r0, r2(2), u(2), v(2), d1, d2, s1(2), s2(2), ri, si, xi(2), sumu2
  integer :: i, count_
  associate (sarray => this%implementation%values)


  n = -1
  s = -1.d0
  x = 0.d0
  select case(xmode)
  case(XSECT_SEGMENT)
     r2 = arg2
     r0 = 1.d0
     u  = r2 - p

  case(XSECT_RAY)
     u  = arg2
     r0 = huge(1.d0)
     r2 = p + r0 * u

  case default
     return
  end select
  sumu2 = sum(u**2)
  if (sumu2 == 0.d0) return
  v  = [-u(2), u(1)]

  count_ = 0
  r = r0
  rbox(1) = min(p(1), r2(1))
  rbox(2) = max(p(1), r2(1))
  zbox(1) = min(p(2), r2(2))
  zbox(2) = max(p(2), r2(2))
  do i=1,this%segments()
     ! S: segment s1->s2
     if (rbox(1) > this%rbox(2,i)) cycle
     if (rbox(2) < this%rbox(1,i)) cycle
     if (zbox(1) > this%zbox(2,i)) cycle
     if (zbox(2) < this%zbox(1,i)) cycle

     s1 = sarray(:,i-1)
     s2 = sarray(:,i)
     d1 = sum((s1-p)*v)
     d2 = sum((s2-p)*v)

     ! x1 and x2 on the same side of L
     if (d1*d2 > 0.d0) cycle

     ! x1 on L
     if (d1 == 0.d0) then
        si = 0.d0

     ! x2 on L
     elseif (d2 == 0.d0) then
        si = 1.d0

     ! intersection between x1 and x2
     else
        si = d1 / (d1 - d2)
     endif

     ! intersection point on L
     xi = s1 + si * (s2-s1)

     ! intersection between p and r2
     ri = sum((xi-p)*u) / sumu2
     if (ri >= 0.d0  .and.  ri <= r0) then
        count_ = count_ + 1
        if (ri <= r) then
           n = i-1
           r = ri
           s = si
           x = xi
        endif
     endif
  enddo
  if (present(count)) count = count_

  end associate
  end subroutine xsect
  !-----------------------------------------------------------------------------
  pure subroutine intersect_polygon2d(this, P, Xp, r, s, n, m)
  class(polygon2d), intent(in)  :: this, P
  real(real64),     intent(out) :: Xp(2), r, s
  integer,          intent(out) :: n, m


  do m=0,P%segments()-1
     call this%intersect(P%node(m), P%node(m+1), XSECT_SEGMENT, Xp, r, s, n)
     if (n >= 0) return
  enddo
  m = -1

  end subroutine intersect_polygon2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function intersects_segment(this, r1, r2, x, r, s, n)
  class(polygon2d), intent(in   ) :: this
  real(real64),     intent(in   ) :: r1(2), r2(2)
  real(real64),     intent(  out), optional :: x(2), r, s
  integer,          intent(  out), optional :: n
  logical                         :: intersects_segment

  real(real64) :: x_(2), r_, s_
  integer      :: n_


  call this%intersect(r1, r2, XSECT_SEGMENT, x_, r_, s_, n_)
  intersects_segment = n_ >= 0
  if (present(x)) x = x_
  if (present(r)) r = r_
  if (present(s)) s = s_
  if (present(n)) n = n_

  end function intersects_segment
  !-----------------------------------------------------------------------------
  function intersects_polygon2d(this, P, x, r, s, n, m)
  class(polygon2d), intent(in   ) :: this, P
  real(real64),     intent(  out), optional :: x(2), r, s
  integer,          intent(  out), optional :: n, m
  logical                         :: intersects_polygon2d

  real(real64) :: x_(2), r_, s_
  integer      :: n_, m_


  call this%intersect_polygon2d(P, x_, r_, s_, n_, m_)
  intersects_polygon2d = m_ >= 0
  if (present(x)) x = x_
  if (present(r)) r = r_
  if (present(s)) s = s_
  if (present(n)) n = n_
  if (present(m)) m = m_

  end function intersects_polygon2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function area(this) result(A)
  !
  ! calculate area enclosed by polygon
  !
  class(polygon2d), intent(in) :: this
  real(real64)                 :: A

  real(real64) :: dA, dx(2), x0(2), x1(2), x2(2)
  integer i
  associate (x => this%implementation%values)


  A = 0.d0
  do i=1,this%segments()
     x1 = x(:,i-1)
     x2 = x(:,i  )
     dx = x2 - x1
     x0 = 0.5d0*(x2 + x1)
     dA = 0.5d0 * (dx(2)*x0(1) - dx(1)*x0(2))

     A = A + dA
  enddo

  end associate
  end function area
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  recursive subroutine shift(this, ds, ierr, parallel)
  use moose_mpi, only: mpi_rank => rank, mpi_nproc => nproc, moose_mpi_sum
  use moose_algorithms, only: intersect_lines_2D
  class(polygon2d), intent(inout) :: this
  real(real64),     intent(in   ) :: ds
  integer,          intent(  out) :: ierr
  logical,          intent(in   ), optional :: parallel

  real(real64), allocatable :: w(:,:), cosa(:), d0(:)
  real(real64) :: u(2), v(2), x(2)
  real(real64) :: x1(2), x2(2), s1, s2, d
  logical :: periodic, parallel_exec
  integer :: i, ix, n, minseg, rank, nproc


  rank = 0;   nproc = 1
  if (present(parallel)) then
     if (parallel) then
        rank = mpi_rank;   nproc = mpi_nproc
     endif
  endif


  ! quick bugfix for ds < 0
  if (ds < 0.d0) then
     call this%reverse()
     call this%shift(abs(ds), ierr, parallel)
     call this%reverse()
     call this%update()
     return
  endif


  ierr = 0
  n = this%segments()
  allocate (w(2, 0:n), cosa(0:n), source=0.d0)


  ! set up boundary nodes
  periodic = this%is_closed()
  minseg = 1;   if (periodic) minseg = 3

  if (rank == 0) then
  if (periodic) then
     u = this%normal(0)
     v = this%normal(n-1)
     w(:,0) = bisector(u, v)
     w(:,n) = w(:,0)
     cosa(0) = sum(u*w(:,0))
     cosa(n) = cosa(0)
  else
     w(:,0) = this%normal(0)
     w(:,n) = this%normal(n-1)
     cosa(0) = 1.d0
     cosa(n) = 1.d0
  endif
  endif


  ! set up internal nodes
  do i=1+rank,n-1,nproc
     u = this%normal(i-1)
     v = this%normal(i)
     w(:,i) = bisector(u, v)
     cosa(i) = sum(u*w(:,i))
  enddo
  if (nproc > 1) then
     call moose_mpi_sum(w)
     call moose_mpi_sum(cosa)
  endif


  allocate (d0(0:nproc-1), source=0.d0);   d0(rank) = huge(1.d0)
  do i=rank,n-1,nproc
     x1 = this%node(i)
     x2 = this%node(i+1)
     call intersect_lines_2D(x1, w(:,i), x2, w(:,i+1), s1, s2, x, ix)
     if (ix /= 0) cycle
     if (s1 < 0.d0  .and.  s2 < 0.d0) cycle

     u = this%normal(i)
     d = sum((x-x1)*u)
     if (d < d0(rank)) d0(rank) = d
  enddo
  if (nproc > 1) then
     call moose_mpi_sum(d0)
  endif


  s1 = min(ds,minval(d0))
  s2 = ds-s1
  do i=0,n
     x = this%node(i)  +  s1 * w(:,i) / cosa(i)
     call this%set_node(i, x)
  enddo
  deallocate (w, cosa, d0)

  call this%remove_duplicate_nodes(1.d-10)
  if (this%segments() < minseg) then
     ierr = 1
     return
  endif
  if (s2 > 0.d0) then
     call this%shift(s2, ierr, parallel)
  endif
  call this%update()

  contains
  !.....................................................................
  function bisector(u, v) result(w)
  real(real64), intent(in) :: u(2), v(2)
  real(real64)             :: w(2)

  real(real64), parameter :: eps = 1.d-8

  real(real64) :: wmod


  w = u + v
  wmod = sqrt(sum(w**2))
  if (wmod < eps) print *, "WARNING: sharp angle detected in polygon2d%shift"
  w = w / max(wmod, eps)

  end function bisector
  !.....................................................................
  end subroutine shift
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(polygon2d), intent(inout) :: this


  call this%polygon%broadcast()
  call proc(0)%broadcast_allocatable(this%rbox)
  call proc(0)%broadcast_allocatable(this%zbox)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(polygon2d), intent(inout) :: this


  deallocate (this%rbox, this%zbox)
  call this%polygon%free()

  end subroutine free
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine sample_mfunc2d_polygon(F, basename, P)
  use moose_grids, only: cgrid
  use moose_mfunc, only: scalar_mfunc, aux_sample_mfunc2d
  class(scalar_mfunc), intent(in) :: F
  character(len=*),    intent(in) :: basename
  type(polygon2d),     intent(in) :: P

  type(cgrid) :: cgrid2d


  cgrid2d = P%cgrid()
  call aux_sample_mfunc2d(F, basename, cgrid2d)

  end subroutine sample_mfunc2d_polygon
  !-----------------------------------------------------------------------------

end module moose_polygon2d
