module moose_uqwork
  use iso_fortran_env
  use moose_uqmesh
  implicit none
  private


  ! workspace for generating an unstructured quadrilateral mesh ................
  type, extends(uqmesh), public :: uqwork
     real(real64), allocatable :: rwork_nodes(:,:), rwork_cells(:,:), rwork_aux(:,:)
     integer, allocatable :: iwork_nodes(:,:), iwork_cells(:,:), iwork_aux(:,:)

     ! number of nodes, cells, auxiliary nodes (from bisection of edges) and auxiliary cell connections
     ! in current and last layer
     integer :: mnodes, mcells, maux, mmulti
     integer :: lnodes, lcells, laux, lmulti

     contains
     procedure :: resize, auto_resize, finalize
     procedure :: io_metadata

     procedure :: any_bad_cell

     procedure :: add_snlayer
  end type uqwork



  logical, public :: debug_mode = .false.


  public :: &
     new_uqwork, loadtxt_uqwork, uqwork_contour, uqwork_layer, &
     area, bad_cell_shape, non_linearity, cleanup_nodes2, smooth_nodes, encode_bsect, mesh_contour


  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new_uqwork(nnodes, ncells, naux, nmulti, &
     iwork_nodes, iwork_cells, iwork_aux, rwork_nodes, rwork_cells, rwork_aux, &
     x1label, x2label, title) result(this)
  !
  ! allocate new workspace for mesh construction
  !
  integer, intent(in) :: nnodes, ncells
  integer, intent(in), optional :: naux, nmulti
  integer, intent(in), optional :: iwork_nodes, iwork_cells, iwork_aux
  integer, intent(in), optional :: rwork_nodes, rwork_cells, rwork_aux
  character(len=*), intent(in), optional :: x1label, x2label, title
  type(uqwork)        :: this


  this%uqmesh = uqmesh(nnodes, ncells, naux, nmulti, x1label=x1label, x2label=x2label, title=title)
  this%mnodes = 0
  this%mcells = 0
  this%maux = 0
  this%mmulti = 0
  this%lnodes = 0
  this%lcells = 0
  this%laux = 0
  this%lmulti = 0

  if (present(iwork_nodes)) allocate (this%iwork_nodes(iwork_nodes, 0:nnodes-1), source = 0)
  if (present(iwork_cells)) allocate (this%iwork_cells(iwork_cells, 0:ncells-1), source = 0)
  if (present(iwork_aux))   allocate (this%iwork_aux(iwork_aux, this%naux()), source = 0)
  if (present(rwork_nodes)) allocate (this%rwork_nodes(rwork_nodes, 0:nnodes-1), source = 0.d0)
  if (present(rwork_cells)) allocate (this%rwork_cells(rwork_cells, 0:ncells-1), source = 0.d0)
  if (present(rwork_aux))   allocate (this%rwork_aux(rwork_aux, this%naux()), source = 0.d0)

  end function new_uqwork
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadtxt_uqwork(filename) result(this)
  character(len=*), intent(in) :: filename
  type(uqwork)                 :: this

  integer :: iwork(8)

  this%uqmesh = loadtxt_uqmesh(filename)
  iwork = this%metadata%getint_rank1("WORKSPACE", 8)
  call this%metadata%remove("WORKSPACE")

  this%mnodes = iwork(1)
  this%mcells = iwork(2)
  this%maux = iwork(3)
  this%mmulti =  iwork(4)
  this%lnodes = iwork(5)
  this%lcells = iwork(6)
  this%laux = iwork(7)
  this%lmulti = iwork(8)

  end function loadtxt_uqwork
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function uqwork_contour(C, m, periodic, iwork_nodes, iwork_cells) result(this)
  !
  ! initialize workspace with *m* nodes along contour *C*
  !
  ! optional input:
  !    iwork_nodes, iwork_cells   additional workspace for nodes and cells
  !
  use moose_error,    only: ERROR
  use moose_analysis, only: interp
  use moose_geometry, only: curve
  class(curve), intent(in) :: C
  integer,      intent(in) :: m
  logical,      intent(in) :: periodic
  integer,      intent(in), optional :: iwork_nodes, iwork_cells
  type(uqwork)             :: this


  this = new_uqwork(m, m, iwork_nodes=iwork_nodes, iwork_cells=iwork_cells)
  this%x(:,0:m-1) = C%discretization(m, tmap=C%arclength_map(), endpoint=.not.periodic)
  this%mnodes = m

  end function uqwork_contour
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function uqwork_layer(C1, C2, m, iwork_nodes, iwork_cells) result(this)
  !
  ! initialize workspace with layer of *m* cells between closed contours *C1* and *C2*
  !
  ! optional input:
  !    iwork_nodes, iwork_cells   additional workspace for nodes and cells
  !
  use moose_error,    only: ERROR
  use moose_analysis, only: interp
  use moose_geometry, only: curve
  class(curve), intent(in) :: C1, C2
  integer,      intent(in) :: m
  integer,      intent(in), optional :: iwork_nodes, iwork_cells
  type(uqwork)             :: this

  real(real64), allocatable :: t(:), e(:)
  integer :: i, ierr


  allocate (t(m), e(m))
  this = new_uqwork(2*m, m, iwork_nodes=iwork_nodes, iwork_cells=iwork_cells)
  this%x(:,m:2*m-1) = C2%arclength_discretization(m, endpoint=.false.)
  call C1%find_footpoints(m, this%x(:,m:2*m-1), t, this%x(:,0:m-1), e, ierr, tmap=C1%arclength_map())
  if (ierr /= 0) call ERROR("computation of footpoints failed", "generate_mmesh2", ierr)
  do i=0,m-1
     this%quads(1,i) = i
     this%quads(2,i) = modulo(i+1,m)
     this%quads(3,i) = m + modulo(i+1,m)
     this%quads(4,i) = m + i
     this%next_cell(:,2,i) = [1, modulo(i+1,m)]
     this%next_cell(:,4,i) = [1, modulo(i-1,m)]
  enddo
  this%mnodes = 2*m
  this%lnodes = m
  this%mcells = m

  end function uqwork_layer
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  function mesh_contour(x) result(C)
  !
  ! construct poloidally closed curve from mesh nodes
  !
  use moose_geometry, only: interp_curve
  real(real64), intent(in) :: x(:,:)
  type(interp_curve)       :: C

  real(real64), allocatable :: xtmp(:,:)
  integer :: m


  m = size(x,2)
  allocate (xtmp(0:m,2))
  xtmp(0:m-1,:) = transpose(x)
  xtmp(m,:) = x(:,1)
  C = interp_curve(xtmp)
  deallocate (xtmp)

  end function mesh_contour
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function area(x1, x2, x3, x4)
  !
  ! compute area of cell cross-section
  !
  use moose_math, only: wedge_product
  real(real64), intent(in) :: x1(2), x2(2), x3(2), x4(2)
  real(real64)             :: area


  area = abs(wedge_product(x3-x2, x2-x1)) + abs(wedge_product(x4-x1, x3-x4))

  end function area
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function bad_cell_shape(x1, x2, x3, x4)
  !
  ! check if cell is non-convex
  !
  use moose_math, only: wedge_product
  real(real64), intent(in) :: x1(2), x2(2), x3(2), x4(2)
  logical                  :: bad_cell_shape

  real(real64) :: a(4)
  integer :: isgn(4)


  bad_cell_shape = .false.
  a(1) = wedge_product(x3-x2, x2-x1)
  a(2) = wedge_product(x4-x1, x3-x4)
  a(3) = wedge_product(x4-x1, x2-x1)
  a(4) = wedge_product(x3-x2, x3-x4)
  isgn = 1;   where (a < 0) isgn = -1
  if (abs(sum(isgn)) < 4) bad_cell_shape = .true.

  end function bad_cell_shape
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function non_linearity(x1, x2, x3, x4)
  !
  ! compute non-linearity of cell shape
  !
  use moose_math, only: wedge_product
  real(real64), intent(in) :: x1(2), x2(2), x3(2), x4(2)
  real(real64)             :: non_linearity

  real(real64) :: a(2), b(2), c(2), d(0:2)


  a = (x3 + x4 - x1 - x2) / 4
  b = (x3 + x2 - x1 - x4) / 4
  c = (x1 + x3 - x2 - x4) / 4

  d(0) = wedge_product(a, b)
  d(1) = wedge_product(b, c)
  d(2) = wedge_product(c, a)

  non_linearity = (abs(d(1)) + abs(d(2))) / abs(d(0))

  end function non_linearity
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine adjust_spacing(t, hmin)
  use moose_error
  use moose_math, only: linspace
  real(real64), intent(inout) :: t(0:)
  real(real64), intent(in   ) :: hmin

  real(real64) :: dl(0:ubound(t,1)), dr(0:ubound(t,1)), davail, dreq, tl, tr
  integer :: i, il, ir, m


  m = size(t) - 1
  if (t(m) - t(0) - m * hmin < 0.d0) call ERROR("hmin is too large", "adjust_spacing")


  ! compute available space for adjusting nodes
  dl = available_space(t, hmin, 1)
  dr = available_space(t, hmin, -1)


  il = 0
  do
     if (il >= m) exit

     ! identify left and right bounds of cluster
     if (dt(il) >= hmin) then
        il = il + 1
        cycle
     endif
     do ir=il+1,m-1
        if (dt(ir) >= hmin) exit
     enddo

     ! check if enough space is available
     davail = dl(il) + dr(ir)
     dreq = (ir - il) * hmin - (t(ir) - t(il))
     if (dreq > davail) call ERROR("not anough space available for cluster", "adjust_spacing")

     ! update nodes in cluster
     tl = t(il) - dreq * dl(il) / davail
     tr = t(ir) + dreq * dr(ir) / davail
     t(il:ir) = linspace(tl, tr, ir - il + 1)

     ! update nodes to the left
     do i=il-1,1,-1
        if (dt(i) > hmin) exit
        t(i) = t(i+1) - max(dt(i), hmin)
     enddo

     ! update nodes to the right
     do i=ir+1,m-1
        if (dt(i-1) > hmin) exit
        ! TODO: check interference with next cluster
        t(i) = t(i-1) + max(dt(i-1), hmin)
     enddo

     ! continue
     il = ir + 1
  enddo

 8001 format("cluster: ",i4,2x,i4,3(2x,f8.3))
  contains
  !.............................................................................
  function available_space(t, hmin,idir) result(d)
  real(real64), intent(inout) :: t(0:)
  real(real64), intent(in   ) :: hmin
  integer,      intent(in   ) :: idir
  real(real64)                :: d(0:ubound(t,1))

  integer, parameter :: k(-1:1) = [0, 0, 1]

  integer :: i, i1, i2, m


  m = size(t) - 1
  d = 0.d0

  i2 = m * k(idir)
  i1 = m * k(-idir) + idir
  do i=i1,i2,idir
     d(i) = d(i-idir) + abs(t(i-idir) - t(i)) - hmin
  enddo

  end function available_space
  !.............................................................................
  function dt(i)
  integer, intent(in) :: i
  real(real64)        :: dt


  dt = t(i+1) - t(i)

  end function dt
  !.............................................................................
  end subroutine adjust_spacing
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine cleanup_nodes(x1, x2, icleanup)
  use moose_error
  use moose_algorithms, only: xsegments
  use moose_utils, only: str
  use moose_txtio, only: savetxt
  real(real64), intent(in   ) :: x1(:,0:)
  real(real64), intent(inout) :: x2(:,0:)
  integer,      intent(  out) :: icleanup

  real(real64) :: dx(2)
  integer :: i, i0, i1, i1m1, i2, i2mod, i2p1, icheck(0:ubound(x1,2)), n


  ! check radial edges
  n = size(x1, 2)
  icheck = 0
  i0 = 0
  i1 = 0
  scan_loop: do
     ! i1: first index of x-cluster
     ! find i2: last index of x-cluster
     xcluster_loop: do i2=i1+1,n-1
        if (xsegments(x1(:,i1), x2(:,i1), x1(:,i2), x2(:,i2))) then
           ! update flags for x-cluster from i1 to i2
           icheck(i1) = icheck(i1) + 1
           icheck(i1+1:i2) = icheck(i1)

           if (i2 == n-1) then
              ! continue search for i2 from beginning
              wrapped_xcluster_loop: do i2mod=0,i1-1
                 if (.not.xsegments(x1(:,i1), x2(:,i1), x1(:,i2mod), x2(:,i2mod))) exit
                 icheck(i1) = icheck(i1) + 1
                 icheck(i1+1:i2) = icheck(i1)
                 icheck(0:i2mod) = icheck(i1)

                 ! head bites the tail of (another) x-cluster ....
                 if (icheck(i2mod+1) /= 0) then
                    call savetxt("ERROR_X1", transpose(x1))
                    call savetxt("ERROR_X2", transpose(x2))
                    call savetxt("ERROR_ICHECK", reshape(1.d0 * icheck, [n, 1]))
                    call ERROR("wrapped x-cluster connects back to itself or another x-cluster")
                 endif
              enddo wrapped_xcluster_loop
              i0 = i2mod + 1
              exit scan_loop
           endif

        ! end of x-cluster detected, continue scan_loop from there
        else
           i1 = i2
           exit xcluster_loop
        endif
     enddo xcluster_loop
     if (i1 >= n-1) exit
  enddo scan_loop


  ! check poloidal edges
  do i=0,n-2
     !if (.not.bad_cell_shape(x1(:,i), x1(:,i+1), x2(:,i+1), x2(:,i))) cycle
     if (.not.xsegments(x1(:,i), x1(:,i+1), x2(:,i), x2(:,i+1))) cycle

     ! already marked?
     if (icheck(i) /= 0  .and.  icheck(i+1) /= 0) cycle

     ! mark cells for adjustment
     if (icheck(i) == 0) then
        ! new x-cluster, or at beginning of x-cluster
        i2 = i + 1 + icheck(i+1)
        if (i2 >= n) call ERROR("i2 >= n")
        icheck(i:i2) = icheck(i+1) + 1
     else
        ! at end of x-cluster (only first icheck in x-cluster is updated)
        i1 = i - icheck(i) + 1
        print *, "WARNING: icheck may be wrong"
        call savetxt("ICHECK_WARNING", reshape(1.d0 * icheck, [n, 1]))
        if (i1 < 0) call ERROR("i1 < 0")
        icheck(i1) = icheck(i) + 1
     endif
  enddo
  icleanup = sum(icheck)
  if (icleanup == 0) return


  ! interpolate nodes in x-cluster
  i1 = i0
  fix_loop: do
     i2 = i1 + icheck(i1)
     if (icheck(i1) /= 0) then
        i1m1 = modulo(i1 - 1, n)
        i2p1 = modulo(i2 + 1, n)
        dx = x2(:, i2p1) - x2(:, i1m1)
        do i=i1,i2
           x2(:, modulo(i, n)) = x2(:, i1m1) + dx * (i - i1 + 1) / (icheck(i1) + 2)
        enddo
     endif
     i1 = i2 + 1
     if (i1 >= n-1) exit
  enddo fix_loop

  end subroutine cleanup_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine cleanup_nodes2(x1, x2, periodic, ncluster)
  use moose_error
  use moose_algorithms, only: xsegments
  use moose_txtio
  use moose_utils
  real(real64), intent(in   ) :: x1(:,0:)
  real(real64), intent(inout) :: x2(:,0:)
  logical,      intent(in   ) :: periodic
  integer,      intent(  out) :: ncluster

  integer :: i, i0, i1, i2, i2mod, ibounds(0:1, 0:ubound(x1, 2)), iflag(0:ubound(x1, 2)), imod, ip1, iterm, &
             k, n, nmax


  n = size(x1, 2)
  nmax = n-1;   if (periodic) nmax = n


  ! 1. scan for x-cells
  iflag = 0
  do i=0,nmax-1
     ip1 = modulo(i+1, n)
     if (xsegments(x1(:,i), x2(:,i), x1(:,ip1), x2(:,ip1))) then
        if (debug_mode) then
           write (99, *) x2(:,i), i
           write (99, *) x2(:,ip1), ip1
        endif
        iflag(i) = iflag(i) + 1
        iflag(ip1) = iflag(ip1) + 1
     endif
  enddo
  if (all(iflag /= 0)) call ERROR("at least one good cell is required", "cleanup_nodes")


  ! 2. get bounds of x-cluster (adjacent x-cells)
  ibounds = 0
  ncluster = 0

  ! 2.1. initially scan backwards for x-cell at periodic boundary
  i0 = 0
  if (periodic .and. iflag(0) /= 0) then
     backward_scan: do imod=n-1,1,-1
        if (iflag(imod) == 0) exit
     enddo backward_scan
     i0 = imod
  endif

  ! 2.2. now scan forwards
  i = i0
  iterm = nmax;   if (periodic) iterm = i0
  scan_loop: do
     ip1 = mod(i + 1, n)

     ! beginning of x-cluster detected
     if (iflag(ip1) /= 0) then
        ibounds(0, ncluster) = i
        cluster_loop: do i2=ip1+1,i0+n-1
           i2mod = mod(i2, n)

           ! end of x-cluster detected
           if (.not.periodic .and. i2mod == nmax) then
              ibounds(1, ncluster) = nmax
              ncluster = ncluster + 1
              exit scan_loop
           endif
           if (iflag(i2mod) == 0) then
              ibounds(1, ncluster) = i2mod
              exit cluster_loop
           endif
        enddo cluster_loop
        ncluster = ncluster + 1

        ! continue from here
        i = i2mod
     endif
     i = mod(i + 1, n)

     ! terminate scan_loop?
     if (i == iterm) exit
  enddo scan_loop


  ! 3. remove x-cluster by interpolating between good nodes
  do k=0,ncluster-1
     if (debug_mode) then
        print *
        print *, "x-cluster detected:", ibounds(:, k)
     endif

     ! 3.1. check for x-edge(s) or z-edges in cluster and interpolate between adjacent nodes
     if (xedge_in_cluster(ibounds(:, k), i1, i2)) then
        if (debug_mode) then
           print *, "xedge in cluster: ", i1, i2
           call savetxt_cluster(i1, i2, "XCLUSTER_DETECTED1")
        endif
        call interpolate_vertices(i1, i2)
        if (debug_mode) call savetxt_cluster(i1, i2, "XCLUSTER_FIXED1")
        if (cluster_is_fixed(i1, i2)) cycle
     endif


     ! 3.2. find x-edge outside of cluster
     if (debug_mode) print *, "scannig for x-edge outside of cluster: ", ibounds(:, k)
     call find_xedge(ibounds(0, k), ibounds(1, k), i1, i2)
     if (debug_mode) call savetxt_cluster(i1, i2, "XCLUSTER_DETECTED2")
     call interpolate_vertices(i1, i2)
     if (debug_mode) call savetxt_cluster(i1, i2, "XCLUSTER_FIXED2")
     if (cluster_is_fixed(i1, i2)) cycle
  enddo

  contains
  !.............................................................................
  subroutine savetxt_cluster(i1, i2, filename)
  integer,          intent(in) :: i1, i2
  character(len=*), intent(in) :: filename

  real(real64) :: xtmp(2,i1:i2)


  if (i2 < n) then
     call savetxt(filename, transpose(x2(:,i1:i2)), append=.true.)
  else
     xtmp(:,i1:n-1) = x2(:,i1:n-1)
     xtmp(:,n:i2) = x2(:,0:i2-n)
     call savetxt(filename, transpose(xtmp), append=.true.)
  endif

  end subroutine savetxt_cluster
  !.............................................................................
  function xedge_in_cluster(ibounds, i1_out, i2_out)
  integer, intent(in   ) :: ibounds(0:1)
  integer, intent(  out) :: i1_out, i2_out
  logical                :: xedge_in_cluster

  real(real64) :: v1(2), v2(2)
  integer :: i1, i1mod, i1p1mod, i1p2mod, i2, i2mod, i2p1mod, imin, imax


  ! pass 1: scan for x-edges
  xedge_in_cluster = .false.
  i2_out = ibounds(0)
  i1_out = ibounds(1);   if (i1_out < i2_out) i1_out = i1_out + n
  imin = i2_out
  imax = i1_out - 1
  do i1=imin,imax
     i1mod = mod(i1, n)
     i1p1mod = mod(i1 + 1, n)
     do i2=i1+2,imax
        i2mod = mod(i2, n)
        i2p1mod = mod(i2 + 1, n)
        if (xsegments(x2(:,i1mod), x2(:,i1p1mod), x2(:,i2mod), x2(:,i2p1mod))) then
           xedge_in_cluster = .true.
           i1_out = min(i1_out, i1)
           i2_out = max(i2_out, i2+1)
        endif
     enddo
  enddo


  ! pass 2: scan for z-edges
  do i1=imin,imax-1
     i1mod = mod(i1, n)
     i1p1mod = mod(i1 + 1, n)
     i1p2mod = mod(i1 + 2, n)
     v1 = x2(:,i1p1mod) - x2(:,i1mod)
     v2 = x2(:,i1p2mod) - x2(:,i1p1mod)
     if (sum(v1*v2) < 0.d0) then
        xedge_in_cluster = .true.
        i1_out = min(i1_out, i1)
        i2_out = max(i2_out, i1+2)
     endif
  enddo

  end function xedge_in_cluster
  !.............................................................................
  subroutine interpolate_vertices(i1, i2)
  integer, intent(in) :: i1, i2

  real(real64) :: dx(2)
  integer :: i, imod


  dx = x2(:, mod(i2, n)) - x2(:, i1)
  do i=i1+1,i2-1
     imod = mod(i, n)
     x2(:, imod) = x2(:, i1) + dx * (i - i1) / (i2 - i1)
  enddo

  end subroutine interpolate_vertices
  !.............................................................................
  function cluster_is_fixed(i1, i2)
  integer, intent(in) :: i1, i2
  logical             :: cluster_is_fixed

  integer :: i, imod, ip1mod


  cluster_is_fixed = .true.
  do i=i1,i2-1
     imod = mod(i, n)
     ip1mod = mod(i + 1, n)
     if (xsegments(x1(:,imod), x2(:,imod), x1(:,ip1mod), x2(:,ip1mod))) then
        cluster_is_fixed = .false.
        exit
     endif
  enddo
  if (debug_mode) print *, "cluster is fixed: ", cluster_is_fixed

  end function cluster_is_fixed
  !.............................................................................
  subroutine find_xedge(iback, ifront, i1, i2)
  integer, intent(in   ) :: iback, ifront
  integer, intent(  out) :: i1, i2

  integer :: iback1, iback2, ifront1, ifront2, k, m


  mloop: do m=1,nmax-1
     ! TODO: truncate for periodic = .false.

     ! 1. extend in forward direction
     iback1 = modulo(iback - m, n)
     iback2 = modulo(iback - m + 1, n)
     do k=0,m
        ifront1 = modulo(ifront + k - 1, n)
        ifront2 = modulo(ifront + k, n)
        if (xsegments(x2(:,ifront1), x2(:,ifront2), x2(:,iback1), x2(:,iback2))) exit mloop
     enddo

     ! 2. extend in backward direction
     ifront1 = modulo(ifront + m - 1, n)
     ifront2 = modulo(ifront + m, n)
     do k=0,m-1
        iback1 = modulo(iback - k, n)
        iback2 = modulo(iback - k + 1, n)
        if (xsegments(x2(:,ifront1), x2(:,ifront2), x2(:,iback1), x2(:,iback2))) exit mloop
     enddo
  enddo mloop
  i1 = iback1
  i2 = ifront2;   if (i2 < i1) i2 = i2 + n

  end subroutine find_xedge
  !.............................................................................
  end subroutine cleanup_nodes2
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine smooth_nodes(x, periodic)
  use moose_math, only: wedge_product
  real(real64), intent(inout) :: x(:,0:)
  logical,      intent(in   ) :: periodic

  real(real64) :: dx(2), s01(2), s12(2), s23(2), w01, w12, w23, wsum
  integer :: i0, i1, i2, i3, imax, k1, k2, m


  m = size(x,2)
  imax = m - 4;   if (periodic) imax = m - 1
  i0 = 0
  do
     if (i0 > imax) exit

     i1 = mod(i0+1, m)
     i2 = mod(i0+2, m)
     i3 = mod(i0+3, m)
     s01 = x(:,i1) - x(:,i0)
     s12 = x(:,i2) - x(:,i1)
     s23 = x(:,i3) - x(:,i2)
     k1 = 1;   if (wedge_product(s01, s12) < 0.d0) k1 = -1
     k2 = 1;   if (wedge_product(s12, s23) < 0.d0) k2 = -1
     dx = x(:,i3) - x(:,i0)

     ! heal consecutive acute angles
     if (sum(s01 * s12) < 0.d0  .and.  sum(s12 * s23) < 0.d0) then
        x(:,i1) = x(:,i0) + 1.d0/3.d0 * dx
        x(:,i2) = x(:,i0) + 2.d0/3.d0 * dx
        cycle
     endif

     ! move to next segment if curvature is consistent
     if (k1 * k2 > 0) then
        i0 = i0 + 1
        cycle
     endif

     ! remove wiggle
     w01 = norm2(s01)
     w12 = norm2(s12)
     w23 = norm2(s23)
     wsum = w01 + w12 + w23
     x(:,i1) = x(:,i0) + w01 / wsum * dx
     x(:,i2) = x(:,i0) + (w01 + w12) / wsum * dx
     i0 = i0 + 2
  enddo

  end subroutine smooth_nodes
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine resize(this, nnodes, ncells, naux, nmulti)
  use moose_table, only: resize_array
  class(uqwork), intent(inout) :: this
  integer,       intent(in   ), optional :: nnodes, ncells, naux, nmulti


  ! nodes
  if (present(nnodes)) then
     this%n = nnodes
     call this%metadata%set("NODES", nnodes)
     call resize_array(this%x, 2, nnodes)
     if (allocated(this%iwork_nodes)) call resize_array(this%iwork_nodes, 2, nnodes)
     if (allocated(this%rwork_nodes)) call resize_array(this%rwork_nodes, 2, nnodes)
  endif


  ! cells
  if (present(ncells)) then
     call this%metadata%set("CELLS", ncells)
     call resize_array(this%quads, 2, ncells)
     call resize_array(this%next_cell, 3, ncells)
     if (allocated(this%iwork_cells)) call resize_array(this%iwork_cells, 2, ncells)
     if (allocated(this%rwork_cells)) call resize_array(this%rwork_cells, 2, ncells)
  endif


  ! auxiliary nodes (from bisecton of edges)
  if (present(naux)) then
  if (naux > 0) then
     call this%metadata%set("AUX_NODES", naux)
     if (this%naux() == 0) then
         allocate (this%aux_nodes(2, naux), source = 0)
     else
         call resize_array(this%aux_nodes, 2, naux)
     endif
     if (allocated(this%iwork_aux)) call resize_array(this%iwork_aux, 2, naux)
     if (allocated(this%rwork_aux)) call resize_array(this%rwork_aux, 2, naux)

  elseif (this%naux() > 0) then
     call this%metadata%remove("AUX_NODES")
     deallocate (this%aux_nodes)
     if (allocated(this%iwork_aux)) deallocate(this%iwork_aux)
     if (allocated(this%rwork_aux)) deallocate(this%rwork_aux)
  endif
  endif


  ! cell edges with multiple neighbors
  if (present(nmulti)) then
  if (nmulti > 0) then
     call this%metadata%set("MULTI_NEXT", nmulti)
     if (this%nmulti() == 0) then
         allocate (this%multi_next(2, nmulti), source = 0)
     else
         call resize_array(this%multi_next, 2, nmulti)
     endif

  elseif (this%nmulti() > 0) then
     call this%metadata%remove("MULTI_NEXT")
     deallocate (this%multi_next)
  endif
  endif

  end subroutine resize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine auto_resize(this, nnodes, ncells, naux, nmulti)
  class(uqwork), intent(inout) :: this
  integer,       intent(in   ), optional :: nnodes, ncells, naux, nmulti

  integer :: m


  if (present(nnodes)) then
     m = this%nnodes()
     if (nnodes > m) call this%resize(nnodes = max(2*m, nnodes))
  endif


  if (present(ncells)) then
     m = this%ncells()
     if (ncells > m) call this%resize(ncells = max(2*m, ncells))
  endif


  if (present(naux)) then
     m = this%naux()
     if (naux > m) call this%resize(naux = max(2*m, naux))
  endif


  if (present(nmulti)) then
     m = this%nmulti()
     if (nmulti > m) call this%resize(nmulti = max(2*m, nmulti))
  endif

  end subroutine auto_resize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine finalize(this)
  !
  ! finalize size of workspace
  !
  class(uqwork), intent(inout) :: this


  call this%resize(this%mnodes, this%mcells, this%maux, this%mmulti)

  ! remove tag for merging
  this%next_cell(2, 3, this%lcells:this%mcells-1) = 0

  end subroutine finalize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function io_metadata(this)
  use moose_dict
  class(uqwork), intent(in) :: this
  type(dict)                :: io_metadata


  io_metadata = this%uqmesh%io_metadata()
  call io_metadata%set("WORKSPACE", [&
     this%mnodes, this%mcells, this%maux, this%mmulti, &
     this%lnodes, this%lcells, this%laux, this%lmulti])

  end function io_metadata
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function any_bad_cell(this, cell_range)
  !
  ! check if any cell has bad shape
  !
  class(uqwork), intent(in) :: this
  integer,       intent(in), optional :: cell_range(2)
  logical                   :: any_bad_cell

  real(real64) :: x1(2), x2(2), x3(2), x4(2)
  integer :: i, i1, i2


  i1 = 0
  i2 = this%mcells - 1
  if (present(cell_range)) then
     i1 = cell_range(1)
     i2 = cell_range(2)
  endif


  any_bad_cell = .false.
  do i=i1,i2
     x1 = this%node(this%quads(1,i))
     x2 = this%node(this%quads(2,i))
     x3 = this%node(this%quads(3,i))
     x4 = this%node(this%quads(4,i))
     if (bad_cell_shape(x1, x2, x3, x4)) then
        any_bad_cell = .true.
        return
     endif
  enddo

  end function any_bad_cell
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine add_snlayer(this, dr, dp, rmin, merging)
! TODO: mask for "inactive cells"
  !
  ! add layer of cells to mesh by expanding contour along secant normal
  !
  !    dr    radial step size
  !    dp    reference width for merging/splitting of cells
  !    rmin  minimal node spacing relative to dp
  !
  ! optional:
  !
  !    merging_     merge two adjacent edges if the sum of the two < dp
  !
  use moose_error
  use moose_utils, only: str
  use moose_txtio, only: savetxt
  use moose_geometry, only: interp_curve
  class(uqwork), intent(inout) :: this
  real(real64),  intent(in   ) :: dr, dp
  real(real64),  intent(in   ), optional :: rmin
  logical,       intent(in   ), optional :: merging

  type(interp_curve) :: C
  real(real64), allocatable :: t(:)
  real(real64) :: n(2), s(2), s1, cluster, x0(2), x1(2), rmin_
  logical :: merging_
  integer :: i, i0, i1, i2, ii, ii1, ii2, inext, icell0, inode0, isub, &
     m, maxmerge, mnext, mmerge, minsert, mk, mremove, nnodes, ncells, naux, nmulti, nmerge, icleanup


  ! user defined parameters
  rmin_ = 0.8d0;   if (present(rmin)) rmin_ = rmin
  merging_ = rmin_ < 0.5d0;   if (present(merging)) merging_ = merging
  if (debug_mode) then
     print *, "DEBUG add_snlayer:"
     print *, "arguments: ", m, dr, dp, rmin, merging_
     print *, "lnodes, mnodes: ", this%lnodes, this%mnodes
     print *, "lcells, mcells: ", this%lcells, this%mcells
     print *, "laux, maux:     ", this%laux, this%maux
     print *, "lmulti, multi:  ", this%lmulti, this%mmulti
     call this%savetxt("UQWORK")
     call savetxt("XOLD", transpose(this%x(:,this%lnodes:this%mnodes+1)))
  endif
  ! check for merged cells in parent layer
  mremove = 0
  i = this%lcells
  do
     if (i == this%mcells) exit
     if (this%next_cell(2, 3, i) > 0) then
        mmerge = 2**this%next_cell(2, 3, i) - 1
        i = i + mmerge
        mremove = mremove + mmerge
     endif
     i = i + 1
  enddo
  m = this%mcells - this%lcells - mremove


  ! increase workspace, if necessary
  nnodes = this%mnodes
  ncells = this%mcells
  naux = this%maux
  nmulti = this%mmulti
  call this%auto_resize(nnodes = nnodes + m, ncells = ncells + m)


  ! construct new layer by moving along normal to secant
  icell0 = this%lcells
  inode0 = this%lnodes
  do i=0,m-1
     i1 = inode0 + modulo(i - 1, m)
     i2 = inode0 + modulo(i + 1, m)

     s = this%x(:,i2) - this%x(:,i1)
     n = [s(2), -s(1)] / norm2(s)
     this%x(:,nnodes + i) = this%x(:,inode0 + i) + n * dr
  enddo
  if (debug_mode) call savetxt("XNEW", transpose(this%x(:,nnodes:nnodes+m-1)))


  ! check for intersecting segments
  i = 0
  do
     !call cleanup_nodes(this%x(:,inode0:inode0+m-1), this%x(:,nnodes:nnodes+m-1), icleanup)
     call cleanup_nodes2(this%x(:,inode0:inode0+m-1), this%x(:,nnodes:nnodes+m-1), .true., icleanup)
     if (debug_mode) call savetxt("XNEW_cleanup"//str(i), transpose(this%x(:,nnodes:nnodes+m-1)))
     if (icleanup == 0) exit
     i = i + 1
     if (i == 8) then
        print *, "WARNING: intersecting segments not fixed after 8 iterations"
        exit
     endif
  enddo


  ! enforce minimal spacing & smooth nodes
  if (rmin_ > 0.d0) then
     C = mesh_contour(this%x(:,nnodes:nnodes+m-1))
     allocate (t(0:m), source=C%t)
     call adjust_spacing(t, rmin_ * dp)
     this%x(:,nnodes:nnodes+m-1) = C%eval(t(0:m-1))
     if (debug_mode) call savetxt("XNEW_spacing"//str(i), transpose(this%x(:,nnodes:nnodes+m-1)))
     deallocate (t)
  endif
  call smooth_nodes(this%x(:,nnodes:nnodes+m-1), .true.)
  if (debug_mode) call savetxt("XNEW_smooth"//str(i), transpose(this%x(:,nnodes:nnodes+m-1)))


  ! define cells
  i = 0
  inext = 0
  mnext = m
  mremove = 0
  maxmerge = int(floor(log(1.d0*m) / log(2.d0)))
  do
     if (i == m) exit
     this%quads(1, ncells) = inode0 + i
     this%quads(4, ncells) = nnodes + inext

     i0 = nnodes + inext
     i1 = nnodes + mod(inext + 1, mnext)
     s1 = norm2(this%x(:,i1) - this%x(:,i0))
     mk = max(0, int(floor(log(s1/dp) / log(2.d0))))
     minsert = 2**mk - 1
     if (this%next_cell(2, 3, icell0) > 0  .and.  minsert > 0) then
        call ERROR("merging from previous layer and splitting at the same edge!")
     endif
     call connect_cells(1)


     ! scan for cluster of short edges for merging
     cluster = s1
     nmerge_loop: do nmerge=0,maxmerge
        ! cluster would exceed remaining number of cells in layer
        if (inext + 2**(nmerge+1) > mnext) exit

        ! extend cluster to next order
        ii1 = i1
        do ii=2**nmerge+1,2**(nmerge+1)
           ii2 = nnodes + mod(inext + ii, mnext)
           cluster = cluster + norm2(this%x(:,ii2) - this%x(:,ii1))
        enddo
        if (cluster > dp) exit
     enddo nmerge_loop


     ! merge 2 (or more) edges
     if (merging_  .and.  nmerge > 0) then
        if (rmin_ >= 0.5d0) print *, "WARNING: merging required despite rmin >= 0.5"
        mmerge = 2**nmerge - 1

        ! delete unnecessary node(s)
        this%x(:,nnodes+inext+1:) = eoshift(this%x(:,nnodes+inext+1:), shift=mmerge, dim=2)
        mnext = mnext - mmerge
        mremove = mremove + mmerge

        ! bisect edge between inext and inext + 1
        call this%auto_resize(naux = naux + mmerge)
        this%quads(3,ncells:ncells+mmerge-1) = bisect(i0, nnodes + mod(inext + 1, mnext), nmerge)
        this%quads(4,ncells+1:ncells+mmerge) = this%quads(3,ncells:ncells+mmerge-1)
        this%next_cell(2, 3, ncells:ncells+mmerge) = nmerge   ! temporary tag

        ! loop over remaining cells in cluster
        do isub=1,mmerge
           this%quads(2, ncells) = inode0 + mod(i + 1, m)
           this%quads(1, ncells + 1) = this%quads(2, ncells)
           i = i + 1
           icell0 = icell0 + 1
           ncells = ncells + 1
           call connect_cells(2)
        enddo

     ! bisect edge
     elseif (minsert > 0) then
        call this%auto_resize(nnodes + mnext + minsert, this%mcells + mnext + minsert + mremove, naux + minsert, nmulti + minsert + 1)
        ! bisect edge between i and i + 1
        this%quads(2,ncells:ncells+minsert-1) = bisect(inode0 + i, inode0 + mod(i + 1, m), mk)
        this%quads(1,ncells+1:ncells+minsert) = this%quads(2,ncells:ncells+minsert-1)

        ! shift remaining nodes in memory
        x0 = this%x(:,i0)
        x1 = this%x(:,i1)
        this%x(:,nnodes+inext+1:) = eoshift(this%x(:,nnodes+inext+1:), shift=-minsert, dim=2)

        ! set up multi_next for parent cell
        this%next_cell(:, 3, icell0) = [1+minsert, nmulti+1]
        this%multi_next(:, nmulti+1) = [encode_bsect(mk, 0), ncells]

        ! insert nodes between inext and inext + 1
        do isub=1,minsert
           this%x(:,nnodes+inext+isub) = x0 + (x1 - x0) * isub / (minsert + 1)
           this%quads(3,ncells) = nnodes + inext + isub
           this%quads(4,ncells+1) = this%quads(3,ncells)
           ncells = ncells + 1
           this%next_cell(:, 1, ncells) = [1, icell0]
           this%multi_next(:, nmulti+1+isub) = [encode_bsect(mk, isub), ncells]
        enddo
        inext = inext + minsert
        mnext = mnext + minsert
        nmulti = nmulti + minsert + 1
     endif

     this%quads(2, ncells) = inode0 + mod(i + 1, m)
     this%quads(3, ncells) = nnodes + mod(inext + 1, mnext)
     i = i + 1
     icell0 = icell0 + 1
     inext = inext + 1
     ncells = ncells + 1
  enddo


  ! update number of nodes
  do i=this%mcells,ncells - 1
     this%next_cell(:, 2, i) = [1, this%mcells + modulo(i - this%mcells + 1, ncells - this%mcells)]
     this%next_cell(:, 4, i) = [1, this%mcells + modulo(i - this%mcells - 1, ncells - this%mcells)]
  enddo
  this%lnodes = this%mnodes
  this%lcells = this%mcells
  this%laux = this%maux
  this%lmulti = this%mmulti
  this%mnodes = nnodes + mnext
  this%mcells = ncells
  this%maux = naux
  this%mmulti = nmulti

  contains
  !.............................................................................
  recursive function bisect(i1, i2, m) result(k)
  !
  ! recursive bisection of edge from vertex i1 to i2 into 2**m segments
  !
  integer, intent(in) :: i1, i2, m
  integer             :: k(2**m-1)

  integer :: i


  i = 2**(m-1)
  naux = naux + 1
  this%aux_nodes(1, naux) = i1
  this%aux_nodes(2, naux) = i2
  k(i) = -naux
  if (m == 1) return

  k(:i-1) = bisect(i1, k(i), m-1)
  k(i+1:) = bisect(k(i), i2, m-1)

  end function bisect
  !.............................................................................
  subroutine connect_cells(k)
  integer :: k

  integer :: i, n, m


  n = this%next_cell(2, 3, icell0)
  if (.not.all(this%next_cell(2, 3, icell0:icell0+n-1) == n)) then
     print *, "inconsistent cluster for merging, icell0 = ", icell0
     print *, this%next_cell(2, 3, icell0:icell0+n-1)
     stop
  endif


  ! connect to cell in previous layer
  if (n == 0) then
     this%next_cell(:, 1, ncells) = [1, icell0]
     this%next_cell(:, 3, icell0) = [1, ncells]

  ! connect to cluster of cells in previous layer
  elseif (n > 0) then
     m = 2**n
     call this%auto_resize(nmulti = nmulti + m)
     this%next_cell(:, 1, ncells) = [m, nmulti + 1]
     do i=0,m-1
        this%next_cell(:, 3, icell0 + i) = [1, ncells]
        this%multi_next(:, nmulti + 1 + i) = [encode_bsect(n, i), icell0 + i]
     enddo
     icell0 = icell0 + m - 1
     nmulti = nmulti + m

  else
     print *, "this should not happen"
     stop
  endif

  end subroutine connect_cells
  !.............................................................................
  end subroutine add_snlayer
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine add_snlayer2(this, dr, dp, rmin, merging, mask)
  !
  ! add layer of cells to mesh by expanding contour along secant normal
  !
  !    dr    radial step size
  !    dp    reference width for merging/splitting of cells
  !    rmin  minimal node spacing relative to dp
  !
  ! optional:
  !
  !    merging_     merge two adjacent edges if the sum of the two < dp
  !
  use moose_error
  use moose_utils, only: str
  use moose_txtio, only: savetxt
  use moose_geometry, only: interp_curve
  class(uqwork), intent(inout) :: this
  real(real64),  intent(in   ) :: dr, dp
  real(real64),  intent(in   ), optional :: rmin
  logical,       intent(in   ), optional :: merging
  integer,       intent(in   ), optional :: mask

  logical :: masked_cell(-1:1, this%lcells:this%mcells-1), periodic
  real(real64) :: x(2, 0:this%mnodes-this%lnodes-1)
  integer, dimension(0:this%mcells-this%lcells-1) :: icell1, icell2, inode1, dcells, dnodes
  integer :: i0, i1, i2, icell, inext, inode, icleanup, j, k, ncluster, nnodes, ncells, naux, nmulti, m


  ! 1. set up masked_cell from this%iwork_cells(mask)
  masked_cell = .false.
  if (present(mask)) then
     ! 1.1. set mask for each cell
     do icell=this%lcells,this%mcells-1
        if (this%iwork_cells(mask, icell) /= 0) masked_cell(0, icell) = .true.
     enddo

     ! 1.2. set mask for left and right neighbor cells
     do icell=this%lcells,this%mcells-1
        ! check if left neighbor does not exist
        if (this%next_cell(1, 4, icell) == 0) then
           masked_cell(-1, icell) = .true.
        else
           inext = this%next_cell(2, 4, icell)
           ! check if left neighbor is masked
           if (masked_cell(0, inext)) masked_cell(-1, icell) = .true.
        endif

        ! check if right neighbor does not exist
        if (this%next_cell(1, 2, icell) == 0) then
           masked_cell(1, icell) = .true.
        else
           inext = this%next_cell(2, 2, icell)
           ! check if right neighbor is masked
           if (masked_cell(0, inext)) masked_cell(1, icell) = .true.
        endif
     enddo
  endif


  ! 2. identify clusters of masked / unmasked cells
  ncluster = 0
  dcells = 0
  dnodes = 2
  do icell=this%lcells,this%mcells-1
     if (masked_cell(0, icell)) cycle
     dcells(ncluster) = dcells(ncluster) + 1

     ! check for start of new cluster
     if (icell == this%lcells  .or.  masked_cell(-1, icell)) then
        icell1(ncluster) = icell
        ! skip over auxiliary nodes
        do
           if (this%quads(4, icell1(ncluster)) >= 0) exit
           if (this%next_cell(1, 4, icell1(ncluster)) == 0) call ERROR("auxiliary node at boundary 4")
           icell1(ncluster) = this%next_cell(2, 4, icell1(ncluster))
        enddo
        inode1(ncluster) = this%quads(4, icell1(ncluster))

     ! count internal nodes (skip auxiliary nodes)
     else
        if (this%quads(4, icell) >= 0) dnodes(ncluster) = dnodes(ncluster) + 1
     endif

     ! check for end of cluster
     if (icell == this%mcells - 1  .or.  masked_cell(1, icell)) then
        icell2(ncluster) = icell
        ! skip over auxiliary nodes
        do
           if (this%quads(3, icell2(ncluster)) >= 0) exit
           if (this%next_cell(1, 2, icell2(ncluster)) == 0) call ERROR("auxiliary node at boundary 2")
           icell2(ncluster) = this%next_cell(2, 2, icell2(ncluster))
        enddo
        ncluster = ncluster + 1
     endif
  enddo

  ! post processing: beginning of first cluster connected to end of last cluster?
  periodic = .false.
  if (icell1(0) == this%lcells  .and.  .not. masked_cell(-1, this%lcells)  .and. &
          all(this%next_cell(:, 4, this%lcells) == [1, this%mcells - 1])) then
     ! merge first and last cluster
     if (ncluster /= 1) then
        icell1(0) = icell1(ncluster-1)
        inode1(0) = inode1(ncluster-1)
        dcells(0) = dcells(0) + dcells(ncluster-1)
        dnodes(0) = dnodes(0) + dnodes(ncluster-1) - 1
        ncluster = ncluster - 1

     ! periodic boundary conditions
     else
        periodic = .true.
        masked_cell(1, this%mcells-1) = .true.
        dnodes(0) = dnodes(0) - 1
     endif
  endif


  ! 3. construct new layer by moving along normal to secant
  m = this%mnodes - this%lnodes
  nnodes = this%mnodes;   this%lnodes = nnodes
  ncells = this%mcells;   this%lcells = ncells
  naux   = this%maux;     this%laux   = naux
  nmulti = this%mmulti;   this%lmulti = nmulti
  do k=0,ncluster-1
     call this%auto_resize(nnodes + dnodes(k))
     ! first node in cluster
     i0 = inode1(k)
     i1 = i0
     i2 = i0 + 1;   if (i2 == this%lnodes) i2 = i2 - m
     icell = icell1(k)
     if (this%next_cell(1, 4, icell) /= 0) i1 = this%quads(4, this%next_cell(2, 4, icell))
     call new_node(this%x(:,i0), this%node(i1), this%x(:, i2)) ! i1 may be virtual node


     ! intermediate nodes
     do inode=1,dnodes(k)-2
        i1 = i0
        i0 = i2
        i2 = i0 + 1;   if (i2 == this%lnodes) i2 = i2 - m
        call new_node(this%x(:, i0), this%x(:, i1), this%x(:, i2))
     enddo


     ! last node in cluster
     i1 = i0
     i0 = i2
     icell = icell2(k)
     if (periodic) then
        i2 = this%quads(3, icell)
     elseif (this%next_cell(1, 2, icell) /= 0) then
        i2 = this%quads(3, this%next_cell(2, 2, icell))
     endif
     call new_node(this%x(:, i0), this%x(:, i1), this%node(i2)) ! i2 may be virtual node


     ! post-processing of nodes
     x = cshift(this%x(:, this%lnodes-m:this%lnodes-1), inode1(k) + m - this%lnodes, dim=2)
     j = 0
     do
        call cleanup_nodes2(x(:, 0:dnodes(k)-1), this%x(:, this%mnodes:nnodes-1), periodic, icleanup)
        if (icleanup == 0) exit
        j = j + 1
        if (j == 8) then
           print *, "WARNING: intersecting segments not fixed after 8 iterations"
           exit
        endif
     enddo
     call smooth_nodes(this%x(:,this%mnodes:nnodes-1), periodic)


     call cell_def(icell1(k), icell, dnodes(k), dcells(k))
     this%mnodes = nnodes
     this%mcells = ncells
     this%maux   = naux
     this%mmulti = nmulti
  enddo

  contains
  !.............................................................................
  subroutine new_node(x0, x1, x2)
  real(real64), intent(in) :: x0(2), x1(2), x2(2)

  real(real64) :: n(2), s(2)


  s = x2 - x1
  n = [s(2), -s(1)] / norm2(s)
  this%x(:, nnodes) = x0 + n * dr
  nnodes = nnodes + 1

  end subroutine new_node
  !.............................................................................
  subroutine cell_def(icell1, icell2, new_nodes, new_cells)
  integer, intent(in   ) :: icell1, icell2
  integer, intent(inout) :: new_nodes, new_cells

  real(real64) :: cluster, ds, x3(2), x4(2)
  integer :: i1, i2, i3, ii3, i4, ii4, ia, ib, iparent, inew, inext, k, kk, kmax, mk, minsert, mmerge, merging_level


  call this%auto_resize(ncells = ncells + new_cells)

  ! 1. connect parent cell(s) to child cell(s)
  iparent = icell1
  k = 0
  do
     ! define nodes in left half of cell
     i1 = this%quads(4, iparent)
     i3 = this%mnodes + k + 1;   if (periodic) i3 = this%mnodes + mod(k + 1, new_nodes)
     i4 = this%mnodes + k
     this%quads(1, ncells) = i1
     this%quads(4, ncells) = i4


     ! evaluate length of upper cell edge -> splitting level
     ds = norm2(this%x(:, i3) - this%x(:, i4))
     mk = max(0, int(floor(log(ds / dp) / log(2.d0))))
     minsert = 2**mk - 1
     if (this%next_cell(2, 3, iparent) > 0  .and.  minsert > 0) then
        print *, "iparent, next_cell, minsert = ", iparent, this%next_cell(2, 3, iparent), minsert
        call ERROR("merging from previous layer and splitting at the same edge!")
     endif
     merging_level = 0
     ! connect to lower neighbor
     call connect_lower(iparent)


     ! A. split edge
     if (minsert > 0) then
        new_cells = new_cells + minsert
        new_nodes = new_nodes + minsert
        call this%auto_resize(this%mnodes + new_nodes, this%mcells + new_cells, naux + minsert, nmulti + minsert + 1)
        ! bisect edge between i1 and i2
        i2 = this%quads(3, iparent)
        this%quads(2, ncells:ncells+minsert-1) = bisect(i1, i2, mk)
        this%quads(1, ncells+1:ncells+minsert) = this%quads(2, ncells:ncells+minsert-1)

        ! shift remaining nodes in memory
        x4 = this%x(:, i4)
        x3 = this%x(:, i3)
        this%x(:, i4+1:) = eoshift(this%x(:, i4+1:), shift=-minsert, dim=2)

        ! set up multi_next for parent cell
        this%next_cell(:, 3, iparent) = [1 + minsert, nmulti + 1]
        this%multi_next(:, nmulti + 1) = [encode_bsect(mk, 0), ncells]

        ! insert nodes between i4 and i3
        do inew=1,minsert
           this%x(:, i4+inew) = x4 + (x3 - x4) * inew / (minsert + 1)
           this%quads(3, ncells) = i4 + inew
           this%quads(4, ncells+1) = this%quads(3, ncells)
           ncells = ncells + 1
           this%next_cell(:, 1, ncells) = [1, iparent]
           this%multi_next(:, nmulti+1+inew) = [encode_bsect(mk, inew), ncells]
        enddo
        nnodes = nnodes + minsert
        nmulti = nmulti + minsert + 1
        k = k + minsert
        i3 = this%mnodes + k + 1;   if (periodic) i3 = this%mnodes + mod(k + 1, new_nodes)


     ! B. check if merging is required
     else
        cluster = ds
        kmax = new_nodes;   if (periodic) kmax = new_nodes + 1
        do
           ! cluster would exceed remaining number of cells in layer
           if (k + 2**(merging_level+1) >= kmax) exit

           ! extend cluster to next order
           ii4 = i3
           do kk=2**merging_level+1,2**(merging_level+1)
              ii3 = this%mnodes + k + kk;   if (periodic) ii3 = this%mnodes + mod(k + kk, new_nodes)
              cluster = cluster + norm2(this%x(:, ii3) - this%x(:, ii4))
              ii4 = ii3
           enddo
           if (cluster > dp) exit
           merging_level = merging_level + 1
        enddo
     endif
     if (merging_level > 0) then
        mmerge = 2**merging_level - 1

        ! delete unnecessary node(s)
        ia = this%mnodes + k + 1
        ib = this%mnodes + new_nodes - 1
        this%x(:, ia:ib) = eoshift(this%x(:, ia:ib), shift=mmerge, dim=2)
        new_nodes = new_nodes - mmerge
        nnodes = nnodes - mmerge

        ! bisect edge between i3 and i4
        i3 = this%mnodes + k + 1;   if (periodic) i3 = this%mnodes + mod(k + 1, new_nodes)
        i4 = this%mnodes + k
        call this%auto_resize(naux = naux + mmerge)
        this%quads(3, ncells:ncells+mmerge-1) = bisect(i4, i3, merging_level)
        this%quads(4, ncells+1:ncells+mmerge) = this%quads(3,ncells:ncells+mmerge-1)
        this%next_cell(2, 3, ncells:ncells+mmerge) = merging_level   ! set temporary tag

        ! loop over remaining cells in cluster
        do kk=1,mmerge
           this%quads(2, ncells) = this%quads(3, iparent)
           this%quads(1, ncells + 1) = this%quads(2, ncells)
           iparent = this%next_cell(2, 2, iparent)
           ncells = ncells + 1
           call connect_lower(iparent)
        enddo
     endif


     ! define nodes in right half of cell (after merging / splitting)
     this%quads(2, ncells) = this%quads(3, iparent)
     this%quads(3, ncells) = i3


     ! move on to next cell
     ncells = ncells + 1
     k = k + 1
     if (iparent == icell2) exit
     iparent = this%next_cell(2, 2, iparent)
  enddo


  ! 2. connect to left and right neighbors (once splitting / merging is done)
  do k=this%mcells,ncells-2
     this%next_cell(:, 2, k) = [1, k+1]
     this%next_cell(:, 4, k+1) = [1, k]
  enddo
  if (periodic) then
     this%next_cell(:, 2, ncells-1) = [1, this%mcells]
     this%next_cell(:, 4, this%mcells) = [1, ncells-1]
  endif

  end subroutine cell_def
  !.............................................................................
  recursive function bisect(i1, i2, m) result(k)
  !
  ! recursive bisection of edge from vertex i1 to i2 into 2**m segments
  !
  integer, intent(in) :: i1, i2, m
  integer             :: k(2**m-1)

  integer :: i


  i = 2**(m-1)
  naux = naux + 1
  this%aux_nodes(1, naux) = i1
  this%aux_nodes(2, naux) = i2
  k(i) = -naux
  if (m == 1) return

  k(:i-1) = bisect(i1, k(i), m-1)
  k(i+1:) = bisect(k(i), i2, m-1)

  end function bisect
  !.............................................................................
  subroutine connect_lower(iparent)
  integer, intent(inout) :: iparent

  integer :: i, n, m


  ! retrieve tag for merging (if applicable)
  n = this%next_cell(2, 3, iparent)
  if (.not.all(this%next_cell(2, 3, iparent:iparent+n-1) == n)) then
     print *, "iparent, ncells, n = ", iparent, ncells, n
     print *, this%next_cell(2, 3, iparent:iparent+n-1)
     call ERROR("inconsistent cell cluster for merging", "add_snlayer2")
  endif


  ! connect to cell in previous layer
  if (n == 0) then
     this%next_cell(:, 1, ncells) = [1, iparent]
     this%next_cell(:, 3, iparent) = [1, ncells]

  ! connect to cluster of cells in previous layer
  elseif (n > 0) then
     m = 2**n
     call this%auto_resize(nmulti = nmulti + m)
     this%next_cell(:, 1, ncells) = [m, nmulti + 1]
     do i=0,m-1
        this%next_cell(:, 3, iparent + i) = [1, ncells]
        this%multi_next(:, nmulti + 1 + i) = [encode_bsect(n, i), iparent + i]
     enddo
     iparent = iparent + m - 1
     nmulti = nmulti + m

  else
     call ERROR("this should not happen", "add_snlayer2")
  endif

  end subroutine connect_lower
  !.............................................................................
  end subroutine add_snlayer2
  !-----------------------------------------------------------------------------

end module moose_uqwork
