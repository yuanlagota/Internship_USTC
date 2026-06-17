module flare_mmesh_vacuum
  use iso_fortran_env
  use moose_mpi
  use moose_geometry
  use flare_mmesh
  implicit none


  type(polygon2d), allocatable :: vacuum_boundary(:), plasma_boundary(:)
  integer, private :: nr, np, nt


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_workspace(iz, iside)
  !
  ! initialize workspace for mesh generator for vacuum domain on boundary
  ! *iside* in zone *iz*
  !
  integer, intent(in) :: iz, iside

  integer :: it


  nr = workspace(iz)%n(1)
  np = workspace(iz)%n(2)
  nt = workspace(iz)%n(3)

  ! 1. allocate vacuum boundary
  allocate (vacuum_boundary(0:nt-1))

  ! 2. boundary of plasma domain
  allocate (plasma_boundary(0:nt-1))
  do it=0,nt-1
     plasma_boundary(it) = mesh_surface(iz, it, workspace(iz)%r_surf_pl_trans_range(iside))
  enddo

  end subroutine init_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_workspace()
  !
  ! clean up workspace for mesh generator for vacuum domain
  !
  integer :: it

  do it=0,nt-1
     call vacuum_boundary(it)%free()
     call plasma_boundary(it)%free()
  enddo
  deallocate (vacuum_boundary, plasma_boundary)

  end subroutine free_workspace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mesh_surface(iz, it, ir) result(P)
  !
  ! polygonal representation of mesh surface
  !
  integer, intent(in) :: iz, it, ir
  type(polygon2d)     :: P


  P = polygon2d(workspace(iz)%n(2)-1)
  associate (x => P%implementation%values)
  x(1,:) = workspace(iz)%r(ir,:,it)
  x(2,:) = workspace(iz)%z(ir,:,it)
  end associate
  call P%update()

  end function mesh_surface
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function mesh_rpath(iz, it, ip, ir1, ir2) result(P)
  !
  ! polygonal representation of radial path
  !
  integer, intent(in) :: iz, it, ip, ir1, ir2
  type(polygon2d)     :: P


  P = polygon2d(ir2-ir1)
  associate (x => P%implementation%values)
  x(1,:) = workspace(iz)%r(ir1:ir2,ip,it)
  x(2,:) = workspace(iz)%z(ir1:ir2,ip,it)
  end associate
  call P%update()

  end function mesh_rpath
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_model_boundary(iz, i)
  !
  ! initialize vacuum boundary for zone *iz* from *i*-th model boundary
  !
  use moose_error
  use moose_utils, only: ordinal, str
  use moose_math,  only: pi
  use flare_model, only: boundary, bfield
  integer, intent(in) :: iz, i

  real(real64) :: phi, x0(2)
  integer :: it, nt


  if (rank == 0) then
     print 1000, ordinal(i)
     print *
     call check_index(i, [1, boundary%nsurfaces], "i")
  endif
 1000 format(3x,"- Initializing from ",a," model boundary")


  nt = workspace(iz)%n(3)
  do it=0,nt-1
     phi = workspace(iz)%phi(it)
     select type (S => boundary%surfaces(i)%geometry)
     type is (axisurf)
        vacuum_boundary(it) = S%P
     type is (torosurf)
        vacuum_boundary(it) = S%polygon2d(phi)
     end select

     if (vacuum_boundary(it)%segments() <= 0) then
        call ERROR("failed to initialize vacuum boundary at toroidal index "//str(it))
     endif

     if (vacuum_boundary(it)%is_closed()) then
        x0 = bfield%equi%magnetic_axis(phi)
        if (vacuum_boundary(it)%winding_number(x0) == -1) call vacuum_boundary(it)%reverse()
     else
        print *, "it, phi [deg] = ", it, phi / pi * 180
        call ERROR("boundary is not closed")
     endif
  enddo

  end subroutine init_model_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_boundary2d(source)
  !
  ! load 2-D boundary contour from *source*
  !
  use moose_utils,  only: split
  use moose_r3grid, only: length_scale
  character(len=*), intent(in) :: source

  character(len=len(source)) :: filename, units
  integer :: it, nt


  if (rank == 0) then
     call split(source, filename, units, set=':', default='m')
     print 1000, trim(filename)
     print *

     vacuum_boundary(0) = polygon2d(filename, scale=length_scale(units))
  endif
  call vacuum_boundary(0)%broadcast()
 1000 format(3x,"- Initializing from user defined 2-D outline ",a)

  nt = size(vacuum_boundary)
  do it=1,nt-1
     vacuum_boundary(it) = vacuum_boundary(0)
  enddo

  end subroutine load_boundary2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_boundary3d(iz, source)
  !
  ! load 3-D boundary for zone *iz* from *source*
  !
  use moose_utils,  only: split
  integer,          intent(in) :: iz
  character(len=*), intent(in) :: source

  character(len=len(source)) :: filename, units
  type(torosurf) :: T
  real(real64)   :: phi
  integer :: it, nt


  if (rank == 0) then
     call split(source, filename, units, set=':', default='m')
     print 1000, trim(filename)
     print *

     T = torosurf(filename, units, convert_units="m")
  endif
  call T%broadcast()
 1000 format(3x,"- Initializing from user defined 3-D outline ",a)


  nt = size(vacuum_boundary)
  do it=0,nt-1
     phi = workspace(iz)%phi(it)
     vacuum_boundary(it) = T%polygon2d(phi)
  enddo
  call T%free()

  end subroutine load_boundary3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_plasma_boundary(iz, iside)
  !
  ! initialize vacuum boundary for zone *iz* from lower/upper radial boundary of
  ! plasma domain
  !
  integer, intent(in) :: iz, iside

  integer :: it, nt


  if (rank == 0) then
     print 1000
     print *
  endif
 1000 format(3x,"- Initializing from boundary of plasma domain")


  nt = workspace(iz)%n(3)
  do it=0,nt-1
     vacuum_boundary(it) = mesh_surface(iz, it, workspace(iz)%r_surf_pl_trans_range(iside))
  enddo

  end subroutine init_plasma_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine expand_boundary(d)
  !
  ! expand boundary contours (right hand shift) by *d* cm
  !
  use moose_error
  real(real64), intent(in) :: d

  character(len=256) :: err
  integer :: ierr, it, nt


  if (rank == 0) then
     print 1000, d
     print *
  endif
 1000 format(3x"- Expanding boundary by ",f0.2," cm (right hand shift)")


  nt = size(vacuum_boundary)
  do it=0,nt-1
     call vacuum_boundary(it)%shift(-d/1.d2, ierr)
     if (ierr /= 0) then
        write (err, 9000) it, ierr
        if (rank == 0) call ERROR(err)
     endif
  enddo
 9000 format("polygon2d%shift failed at toroidal index ",i0," with ierr = ",i0)

  end subroutine expand_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine rzbuffer(width)
  !
  ! adjust nodes on boundary for minimum distance of *width* cm to plasma boundary.
  !
  real(real64), intent(in) :: width

  integer :: it, nt


  if (rank == 0) then
     print 1000, width
     print *
  endif
 1000 format(3x,"- Adjusting nodes on boundary for buffer zone of ",f0.2," cm")


  nt = size(vacuum_boundary)
  do it=0,nt-1
     call make_buffer_zone(vacuum_boundary(it), plasma_boundary(it), width/1.d2)
  enddo

  end subroutine rzbuffer
  !-----------------------------------------------------------------------------
  subroutine make_buffer_zone(P_adjust, P_fixed, width)
  use moose_error
  use moose_rlist
  type(polygon2d), intent(inout) :: P_adjust
  type(polygon2d), intent(in   ) :: P_fixed
  real(real64),    intent(in   ) :: width

  type(polygon2d) :: pwork(2)
  type(rlist2) :: zip
  real(real64), allocatable :: work(:,:)
  integer, allocatable :: iwork(:,:)
  real(real64) :: r, s, x(2)
  logical :: cycle_iwork
  integer :: count, i, ierr, j, k(2), n(2), ncross


  pwork(1) = P_adjust
  pwork(2) = P_fixed
  call pwork(2)%shift(-width, ierr)
  if (ierr /= 0) call ERROR("failed to create margin for plasm boundary")


  ! find intersections
  n(1) = pwork(1)%segments()
  n(2) = pwork(2)%segments()
  allocate (work(2,0:n(1)-1), iwork(2,0:n(1)-1))
  ncross = 0
  iwork = -1
  do i=0,n(1)-1
     call pwork(2)%intersect(pwork(1)%node(i), pwork(1)%node(i+1), XSECT_SEGMENT, x, r, s, j, count=count)
     if (j < 0  .or.  mod(count, 2) == 0) cycle

     work(:, ncross) = x
     iwork(:, ncross) = [i, j]
     ncross = ncross + 1
  enddo
  if (mod(ncross, 2) == 1) then
     print *, "ncross = ", ncross
     call pwork(1)%savetxt("ERROR_PWORK1")
     call pwork(2)%savetxt("ERROR_PWORK2")
     call P_fixed%savetxt("ERROR_PWORK0")
     call ERROR("odd number of margin points detected")
  endif
  cycle_iwork = any(iwork(2,1:ncross-1) - iwork(2,0:ncross-2) < 0)


  ! zip segments
  zip = rlist2()
  k = 0  ! keep track of current node
  j = 0  ! keep track of intersection points
  i = 1  ! switch between workspaces 1 <-> 2
  if (P_fixed%get_distance(pwork(1)%node(0)) < P_fixed%get_distance(pwork(2)%node(0))) i = 2
  do
     !if (k(i) > n(i)) exit   ! always exit for open contours
     if (k(i) > n(i)) then
        if (i == 2  .and.  cycle_iwork) then
           k(i) = 1
        else
           exit
        endif
     endif
     call zip%append(pwork(i)%node(k(i)))

     ! at intersection point?
     if (k(i) == iwork(i, j)) then
        call zip%append(work(:, j))
        i = 3 - i
        k(i) = iwork(i, j) + 1
        j = j + 1

     ! otherwise, continue with next point
     else
        k(i) = k(i) + 1
     endif
  enddo
  P_adjust = polygon2d(zip)
  call zip%free()

  end subroutine make_buffer_zone
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine interpolated_mesh(iz, iside)
  !
  ! construct interpolated mesh in zone *iz* between vacuum and plasma boundary
  ! based on spacing along plasma boundary
  !
  integer, intent(in) :: iz, iside

  type(interp_curve) :: C
  real(real64), allocatable :: s(:), xtmp(:,:,:)
  real(real64) :: r, x1(2), x2(2)
  integer :: ir, ir1, ir2, ir11, ir22, ip, it, np, nt


  if (rank == 0) then
     print 1000
     print *
  endif
 1000 format(3x,"- Constructing interpolated mesh based on spacing along plasma boundary")


  np = workspace(iz)%n(2)
  nt = workspace(iz)%n(3)

  ir1 = workspace(iz)%r_surf_pl_trans_range(iside)
  ir2 = 0;   if (iside == 2) ir2 = workspace(iz)%n(1) - 1
  ir11 = min(ir1,ir2)
  ir22 = max(ir1,ir2)

  allocate (s(0:np-1), xtmp(ir11:ir22, 0:np-1, 2), source=0.d0)
  do it=0,nt-1
     C = interp_polygon(vacuum_boundary(it))
     s = plasma_boundary(it)%accumulated_lengths(normalized=.true.) * vacuum_boundary(it)%length()

     xtmp = 0.d0
     do ip=rank,np-1,nproc
        x1 = plasma_boundary(it)%node(ip)
        x2 = C%eval(s(ip))

        do ir=ir1,ir2,(-1)**iside
           r = 1.d0 * (ir - ir1) / (ir2 - ir1)
           xtmp(ir,ip,:) = x1 + r * (x2 - x1)
        enddo
     enddo
     call moose_mpi_sum(xtmp)
     workspace(iz)%r(ir11:ir22,:,it) = xtmp(:,:,1)
     workspace(iz)%z(ir11:ir22,:,it) = xtmp(:,:,2)

     call C%free()
  enddo
  deallocate (s, xtmp)

  end subroutine interpolated_mesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine normal_trace(iz, iside, m, xlike_block, xlike_cell)
  !
  ! divide domain into m blocks by following normal direction from plasma boundary
  !
  use moose_error
  use moose_algorithms, only: xsegments
  use moose_math,     only: linspace
  use moose_analysis, only: interp, pchip
  use flare_control,  only: verbose, diagnostic_mode
  use moose_qmesh
  use moose_utils, only: str
  integer, intent(in)    :: iz, iside, m
  logical, intent(  out) :: xlike_block, xlike_cell

  type(interp_curve) :: C
  type(interp) :: wfunc
  type(qmesh) :: blocks
  real(real64), allocatable :: lv(:), lp(:)
  character(len=128) :: tmp
  real(real64) :: r, s, w(0:m), x(2,0:1), x1(2), x2(2), x3(2), x4(2), xi(2)
  integer :: i, ir, ir1, ir2, ir11, ir22, ip, it, iu, n, nvac


  if (rank == 0) then
     print 1000
     print 1001, m
     print *
  endif
 1000 format(3x,"- Constructing interpolated mesh by tracing normal direction from plasma boundary")
 1001 format(8x,"number of blocks: ", i0)


  np = workspace(iz)%n(2)
  nt = workspace(iz)%n(3)
  ir1 = workspace(iz)%r_surf_pl_trans_range(iside)
  ir2 = 0;   if (iside == 2) ir2 = workspace(iz)%n(1) - 1
  ir11 = min(ir1,ir2)
  ir22 = max(ir1,ir2)

  xlike_block = .false.
  xlike_cell = .false.
  do it=0,nt-1
     if (verbose) print *, it

     ! find intersection points & coordinates on vacuum boundary
     nvac = vacuum_boundary(it)%segments()
     allocate (lv(0:nvac), source = vacuum_boundary(it)%accumulated_lengths())
     allocate (lp(0:np-1), source = plasma_boundary(it)%accumulated_lengths())
     C = interp_curve(lp, plasma_boundary(it)%nodes())
     blocks = qmesh(2, m+1)
     do i=0,m
        x  = C%deriv(min(C%b*i/m, C%b), 1)
        x1 = x(:,0)
        x2 = x1 + [x(2,1), -x(1,1)]

        call vacuum_boundary(it)%intersect(x1, [x(2,1), -x(1,1)], XSECT_RAY, xi, r, s, n)
        if (n < 0) then
           print *, "it = ", it
           call plasma_boundary(it)%savetxt("ERROR_PLASMA_BOUNDARY")
           call vacuum_boundary(it)%savetxt("ERROR_VACUUM_BOUNDARY")
           open  (newunit=iu, file="ERROR_NORMAL")
           write (iu, *) x1
           write (iu, *) x2
           close (iu)
           call ERROR("plasma boundary normal does not intersect vacuum boundary")
        endif
        blocks%x(0,i,:) = x1
        blocks%x(1,i,:) = xi
        w(i) = lv(n) + s * (lv(n+1) - lv(n))
     enddo
     if (diagnostic_mode) call blocks%savetxt("BLOCKS_"//str(it))
     call C%free()


     ! test for x-like blocks
     do i=0,m-1
        if (xsegments(blocks%x(0,i,:), blocks%x(1,i,:), blocks%x(0,i+1,:), blocks%x(1,i+1,:))) then
           if (.not.xlike_block) print *, "x-like block detected: iblock, it = ", i, it
           xlike_block = .true.
        endif
     enddo


     ! construct increasing order & continuous mapping
     do i=1,m
        if (w(i) - w(i-1) < 0.d0) w(i:) = w(i:) + lv(nvac)
     enddo
     wfunc = pchip(linspace(0.d0, lp(np-1), m+1), w)
     if (diagnostic_mode) call wfunc%savetxt("WFUNC_"//str(it))


     ! interpolate mesh nodes
     C = linear_interp(lv, vacuum_boundary(it)%nodes())
     do ip=0,np-1
        x1 = plasma_boundary(it)%node(ip)
        x2 = C%eval(mod(wfunc%eval(lp(ip)), lv(nvac)))
        do ir=ir1,ir2,(-1)**iside
           r  = 1.d0 * (ir - ir1) / (ir2 - ir1)
           xi = x1 + r * (x2 - x1)
           workspace(iz)%r(ir,ip,it) = xi(1)
           workspace(iz)%z(ir,ip,it) = xi(2)
        enddo
     enddo
     deallocate (lp, lv)
     call C%free()


     ! check for x-like cells
     do ip=0,np-2
        x1(1) = workspace(iz)%r(ir1,ip,it)
        x1(2) = workspace(iz)%z(ir1,ip,it)
        x2(1) = workspace(iz)%r(ir1,ip+1,it)
        x2(2) = workspace(iz)%z(ir1,ip+1,it)
        x3(1) = workspace(iz)%r(ir2,ip+1,it)
        x3(2) = workspace(iz)%z(ir2,ip+1,it)
        x4(1) = workspace(iz)%r(ir2,ip,it)
        x4(2) = workspace(iz)%z(ir2,ip,it)
        if (bad_shape(x1, x2, x3, x4)) then
           if (.not.xlike_cell) print *, "x-like cell detected: ip, it = ", ip, it
           xlike_cell = .true.
        endif
     enddo
  enddo

  end subroutine normal_trace
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine auto_blocks(iz, iside)
  use moose_error
  integer, intent(in) :: iz, iside

  logical :: xlike_block, xlike_cell
  integer :: m, mmin, mmax


  m = 16
  mmin = 1
  mmax = workspace(iz)%n(2) - 1
  do
     call normal_trace(iz, iside, m, xlike_block, xlike_cell)
     ! success?
     if (.not.xlike_block  .and.  .not.xlike_cell) exit

     ! unrecoverable error?
     !if ((xlike_block .and. xlike_cell) .or. mmin == mmax) then
     if (mmax - mmin <= 1) then
        call ERROR("automatic construction of interpolated blocks failed")
     endif

     ! decrease number of blocks if x-like blocks are detected
     if (xlike_block) then
        mmax = m
        m = (mmin + m) / 2

     ! increase number of blocks if x-like cells are detected
     elseif (xlike_cell) then
        mmin = m
        m = (m + mmax) / 2
     endif
  enddo

  end subroutine auto_blocks
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine quasi_orthogonal_mesh(iz, iside)
  !
  ! construct quasi orthogonal mesh in vacuum domain on side *iside* in zone *iz*
  !
  use moose_error
  use moose_utils,   only: str
  use moose_math,    only: pi
  use moose_qmesh_generator
  use flare_control, only: verbose
  integer, intent(in) :: iz, iside

  type(qmesh) :: tmp
  type(interp_curve) :: C
  real(real64), allocatable :: xtmp(:,:,:,:)
  real(real64) :: phi
  integer :: ir1, ir11, ir2, ir22, it, np, nt


  if (rank == 0) then
     print 1000
     print *
  endif
  call mpi_barrier_world()
 1000 format(3x,"- Constructing quasi-orthogonal mesh ...")


  np = workspace(iz)%n(2)
  nt = workspace(iz)%n(3)
  ir1 = workspace(iz)%r_surf_pl_trans_range(iside)
  ir2 = 0;   if (iside == 2) ir2 = workspace(iz)%n(1) - 1
  ir11 = min(ir1,ir2)
  ir22 = max(ir1,ir2)


  if (rank == 0) print 2000
  allocate (xtmp(ir11:ir22, 0:np-1, 0:nt-1, 2), source=0.d0)
  tmp = qmesh(ir22-ir11+1,np)
  do it=rank,nt-1,nproc
     phi = workspace(iz)%phi(it) / pi * 180.d0
     print *, it, phi
     if (.not.vacuum_boundary(it)%is_closed()) then
        call vacuum_boundary(it)%savetxt("ERROR_VACUUM_BOUNDARY")
        call ERROR("contour of vacuum boundary is not closed (check ERROR_VACUUM_BOUNDARY)")
     endif

     ! initialize plasma boundary in tmp
     tmp%x(ir1-ir11,:,1) = workspace(iz)%r(ir1,:,it)
     tmp%x(ir1-ir11,:,2) = workspace(iz)%z(ir1,:,it)

     ! construct mesh
     if (verbose) debug_dummyU = "DUMMY_POTENTIAL"//str(it)
     C = interp_polygon(vacuum_boundary(it))
     call construct_submesh(tmp, ir1-ir11, ir2-ir11, C, container=16)
     call C%free()

     ! save mesh nodes to xtmp
     xtmp(ir11:ir22,:,it,:) = tmp%x
  enddo
  call moose_mpi_sum(xtmp)
 2000 format(8x,"index, toroidal angle [deg]")


  ! copy mesh nodes to workspace
  workspace(iz)%r(ir11:ir22,:,:) = xtmp(:,:,:,1)
  workspace(iz)%z(ir11:ir22,:,:) = xtmp(:,:,:,2)

  end subroutine quasi_orthogonal_mesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine indices_for_plate_aligned_mesh(iz, iside, j, jsign, k)
  use flare_model, only: bfield
  integer, intent(in   ) :: iz, iside
  integer, intent(  out) :: j, jsign, k

  integer :: kselect(-1:1), ksign


  kselect(-1) = 0
  kselect( 1) = workspace(iz)%n(3)-1
  ksign = int(bfield%equi%Bp_sign * bfield%equi%Bt_sign)


  select case(iside)
  ! lower boundary
  case (1)
     j = 0;   jsign = -1
     k = kselect(ksign)

  ! upper boundary
  case (2)
     j = workspace(iz)%n(2)-1;   jsign = 1
     k = kselect(-ksign)

  end select
  !print *, "indices_for_plate_aligned_mesh(", iz, ", ", iside, ")"
  !print *, j, jsign, k

  end subroutine indices_for_plate_aligned_mesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine interpolate_nodes(iz, i1, i2, j, k)
  !
  ! interpolate nodes in radial direction between *i1* and *i2* in zone *iz* at
  ! poloidal index *j* and toroidal index *k*
  !
  integer, intent(in) :: iz, i1, i2, j, k

  type(polygon2d) :: rpath
  real(real64), allocatable :: r(:)
  real(real64) :: x(2), x1(2), x2(2)
  integer :: i


  rpath = mesh_rpath(iz, k, j, i1, i2)
  allocate (r(i1:i2), source=rpath%accumulated_lengths(normalized=.true.))


  x1 = [workspace(iz)%r(i1,j,k), workspace(iz)%z(i1,j,k)]
  x2 = [workspace(iz)%r(i2,j,k), workspace(iz)%z(i2,j,k)]
  do i=i1+1,i2-1
     x = x1 + r(i) * (x2-x1)
     workspace(iz)%r(i,j,k) = x(1)
     workspace(iz)%z(i,j,k) = x(2)
  enddo


  deallocate (r)
  call rpath%free()

  end subroutine interpolate_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine update_nodes(iz, j0, jsign, k0)
  !
  ! update nodes in zone *iz* at poloidal index *j0* in toroidal direction based
  ! nodes at toroidal index *k0*
  !
  use flare_mmesh_parameters
  integer, intent(in) :: iz, j0, jsign, k0

  real(real64) :: x(2), v(2)
  integer :: i, k


  do i=0,workspace(iz)%n(1)-1
     x = [workspace(iz)%r(i,j0,k0), workspace(iz)%z(i,j0,k0)]
     v = x - [workspace(iz)%r(i,j0+jsign,k0), workspace(iz)%z(i,j0+jsign,k0)]
     v = v / sqrt(sum(v**2)) * polo_ext / 100

     do k=0,workspace(iz)%n(3)-1
        workspace(iz)%r(i,j0,k) = x(1) + v(1)
        workspace(iz)%z(i,j0,k) = x(2) + v(2)
     enddo
  enddo

  end subroutine update_nodes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine close_mesh_domain(iz, iside)
  !
  ! extend boundary cells on poloidal side *iside* in zone *iz* order to close
  ! mesh domain
  !
  use moose_math
  use flare_control, only: verbose
  integer, intent(in) :: iz, iside

  type(polygon2d) :: rpath
  real(real64) :: t(2), v(2), a, amax
  integer :: i, ii, imax, j, jsign, jselect(2), k, nr


  ! set reference path
  call indices_for_plate_aligned_mesh(iz, iside, j, jsign, k)
  nr = workspace(iz)%n(1)-1
  rpath = mesh_rpath(iz, k, j, 0, nr)


  i = 0
  do
     if (i == nr) exit

     ! tangent vector from node i -> i+1
     t = rpath%node(i+1) - rpath%node(i)
     ! find max. deviation from line
     amax = 0.d0
     do ii=i+2,nr
        v = rpath%node(ii) - rpath%node(i)
        a = -jsign * (atan2(t(2), t(1)) - atan2(v(2), v(1)));   if (abs(a) > pi) a = a - sign(pi2,a)
        if (a > amax) then
           amax = a
           imax = ii
        endif
     enddo

     ! adjustment necessary?
     if (amax > 0.d0) then
        call interpolate_nodes(iz, i, imax, j, k)
        if (verbose) print 1001, i, imax
        i = imax
        cycle
     endif

     ! continue with next node
     i = i + 1
  enddo
  call update_nodes(iz, j, -jsign, k)
 1001 format(8x,'interpolating nodes in concave domain ir = ',i0,' -> ',i0)


  ! update p_surf_pl_trans_range
  jselect = [1, workspace(iz)%n(2)-2]
  workspace(iz)%p_surf_pl_trans_range(iside) = jselect(iside)

  ! cleanup
  call rpath%free()

  end subroutine close_mesh_domain
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine close_mesh_domain_usr(iz, iside, source)
  use moose_utils,  only: split
  use moose_r3grid, only: length_scale
  integer,          intent(in) :: iz, iside
  character(len=*), intent(in) :: source

  character(len=len(source)) :: filename, units
  type(polygon2d) :: P
  type(interp_curve) :: C
  real(real64), allocatable :: t(:)
  real(real64) :: dr, dz, x(2)
  integer :: i, j0, j1, jselect(2), k, nr


  ! load contour from file
  if (rank == 0) then
     call split(source, filename, units, set=':', default='m')
     print 1000, trim(filename)
     print *

     P = polygon2d(filename, scale=length_scale(units))
  endif
  call P%broadcast()
  C = interp_polygon(P)
 1000 format(3x,"- Setting divertor leg closure from user defined 2-D contour ",a)


  ! update p_surf_pl_trans_range
  if (iside == 1) then
     j0 = 0
     j1 = 1
  else
     j0 = workspace(iz)%n(2) - 1
     j1 = j0 - 1
  endif
  workspace(iz)%p_surf_pl_trans_range(iside) = j1


  nr = workspace(iz)%n(1)
  allocate (t(0:nr-1))
  do k=0,workspace(iz)%n(3)-1
     ! pass 1: compute t
     do i=1,nr-1
        dr = workspace(iz)%r(i,j1,k) - workspace(iz)%r(i-1,j1,k)
        dz = workspace(iz)%z(i,j1,k) - workspace(iz)%z(i-1,j1,k)
        t(i) = t(i-1) + sqrt(dr**2 + dz**2)
     enddo
     t = t / t(nr-1)

     ! pass 2: interpolate
     do i=0,nr-1
        x = C%eval(C%a + t(i) * (C%b - C%a))
        workspace(iz)%r(i,j0,k) = x(1)
        workspace(iz)%z(i,j0,k) = x(2)
     enddo
  enddo

  end subroutine close_mesh_domain_usr
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine sample_Bmod_in_vacuum_domain(iz, iside)
  !
  ! sample magnetic field strength in vacuum domain of radial side *iside* in
  ! zone *iz*
  !
  use moose_math
  use flare_control
  use flare_model, only: bfield
  integer, intent(in) :: iz, iside

  real(real64), allocatable :: b(:,:,:)
  real(real64) :: x(3)
  integer :: i, ir, ir1, ir11, ir2, ir22, ip, it, n, np, nt


  if (rank == 0) print 1000
 1000 format(3x,"- Sampling magnetic field strength:")


  ir1 = workspace(iz)%r_surf_pl_trans_range(iside)
  ir2 = 0;   if (iside == 2) ir2 = workspace(iz)%n(1) - 1
  ir11 = min(ir1,ir2)
  ir22 = max(ir1,ir2)

  np = workspace(iz)%n(2)
  nt = workspace(iz)%n(3)
  allocate (b(ir11:ir22, 0:np-1, 0:nt-1), source=0.d0)

  i = 0
  n = (ir22-ir11+1) * np * nt
  call progress_bar(0, n)
  do it=0,nt-1
  do ip=0,np-1
  do ir=ir11,ir22
     i = i + 1;   if (mod(i,nproc) /= rank) cycle

     x(1) = workspace(iz)%r(ir,ip,it)
     x(2) = workspace(iz)%z(ir,ip,it)
     x(3) = workspace(iz)%phi(it)

     b(ir,ip,it) = bfield%bmod(x)
     call progress_bar(i, n)
  enddo
  enddo
  enddo
  call moose_mpi_sum(b)
  call finalize_progress_bar()


  workspace(iz)%b(ir11:ir22,:,:) = b
  deallocate (b)

  end subroutine sample_Bmod_in_vacuum_domain
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine mesh_generator(iz, iside, recipe, sample_Bmod)
  use moose_mpi
  use moose_error
  use moose_utils, only: str, substring, nsubstrings
  use flare_model
  use flare_tasks
  integer,          intent(in) :: iz, iside
  character(len=*), intent(in) :: recipe
  logical,          intent(in) :: sample_Bmod

  character(len=*), parameter :: side(2) = ['lower', 'upper']

  character(len=256) :: err
  character(len=len(recipe)) :: cmd, dummy, instruction, filename, tag
  real(real64) :: phi, d, w
  logical :: xlike_block, xlike_cell
  integer :: i, it, ios, n, m


  call begin_task()
  if (rank == 0) then
     print 1000, iz, side(iside)
     print *
  endif
  n = nsubstrings(recipe);   if (n == 0) call ERROR("undefined instructions")
  call init_workspace(iz, iside)
 1000 format(1x,"Constructing vacuum domain in zone ",i0," on ",a," radial boundary")


  ! 1. initialize vacuum boundary geometry
  instruction = substring(recipe, 1)
  read (instruction, *) cmd
  select case (cmd)
  ! A. initialize from model boundary
  case ('MODEL_BOUNDARY')
     read  (instruction, *, iostat=ios) dummy, i;   call assert_ios()
     call init_model_boundary(iz, i)


  ! B. initialize from user defined 2-D outline
  case ('LOAD2D')
     read  (instruction, *, iostat=ios) dummy, filename;   call assert_ios()
     call load_boundary2d(filename)


  ! C. initialize from user defined 3-D outline
  case ('LOAD3D')
     read  (instruction, *, iostat=ios) dummy, filename;   call assert_ios()
     call load_boundary3d(iz, filename)


  ! D. initialize from boundary of plasma domain
  case ('PLASMA_BOUNDARY')
     call init_plasma_boundary(iz, iside)


  ! invalid command
  case default
     write (err, 9001) trim(cmd)
     if (rank == 0) call ERROR(err)

  end select
 9001 format("invalid command ",a," for boundary initialization")


  ! 2. process boundary geometry
  do i=2,n-1
     instruction = substring(recipe,i)
     if (instruction == "") cycle

     read (instruction, *) cmd
     select case (cmd)
     ! A. expand boundary contours (right hand shift)
     case ('EXPAND')
        read (instruction, *, iostat=ios) dummy, d;   call assert_ios()
        call expand_boundary(d)


     ! B. adjust nodes to include buffer zone
     case ('RZBUFFER')
        read (instruction, *, iostat=ios) dummy, w;   call assert_ios()
        call rzbuffer(w)


     ! C. flip orientation
     case ('FLIP')
        do it=0,size(vacuum_boundary)-1
           call vacuum_boundary(it)%reverse()
        enddo


     ! Z. export geometry
     case ('EXPORT')
        read (instruction, *, iostat=ios) dummy, tag;   call assert_ios()
        do it=0,size(vacuum_boundary)-1
           call vacuum_boundary(it)%savetxt("VACUUM_BOUNDARY_"//trim(tag)//"_"//str(iz)//"_"//str(it)//".txt")
        enddo

     ! invalid command
     case default
        write (err, 9002) trim(cmd)
        call ERROR(err)

     end select
  enddo
  if (diagnostic_mode .and. rank == 0) then
     print 3001
     do i=0,nt-1
        call plasma_boundary(i)%savetxt("PLASMA_BOUNDARY"//str(i))
        call vacuum_boundary(i)%savetxt("VACUUM_BOUNDARY"//str(i))
     enddo
  endif
 3001 format(3x,"- Saving contour to VACUUM_BOUNDARY*")
 9002 format("invalid command ",a)


  ! 3.0. verify contours
  do it=0,nt-1
     if (vacuum_boundary(it)%intersects(plasma_boundary(it))) then
        call vacuum_boundary(it)%savetxt("ERROR_VACUUM_BOUNDARY")
        call plasma_boundary(it)%savetxt("ERROR_PLASMA_BOUNDARY")
        call ERROR("vacuum boundary at it = "//str(it)//" intersects plasma boundary")
     endif
  enddo

  ! 3.1. construct mesh
  instruction = substring(recipe,n)
  read (instruction, *) cmd
  select case (trim(cmd))
  ! A. export relative spacing to boundary & interpolate mesh nodes
  case ('INTERPOLATE')
     call interpolated_mesh(iz, iside)


  ! B. divide domain into m blocks by following normal direction from plasma boundary
  case ('NORMAL_TRACE')
     read  (instruction, *, iostat=ios) dummy, m;   call assert_ios()
     call normal_trace(iz, iside, m, xlike_block, xlike_cell)


  ! B.2.
  case ('AUTO_BLOCKS')
     call auto_blocks(iz, iside)


  ! C. construct quasi-orthogonal mesh
  case ('QUASI_ORTHOGONAL')
     call quasi_orthogonal_mesh(iz, iside)


  ! invalid command
  case default
     write (err, 9003) trim(cmd)
     call ERROR(err)

  end select
 9003 format("invalid mesh construction command ",a)


  ! 4. sample Bmod
  if (sample_Bmod) call sample_Bmod_in_vacuum_domain(iz, iside)

  call free_workspace()
  call finalize_task()

  contains
  !.............................................................................
  subroutine assert_ios()


  if (ios == 0) return
  call ERROR("missing or invalid parameter for "//trim(cmd))

  end subroutine
  !.............................................................................
  end subroutine mesh_generator
  !-----------------------------------------------------------------------------

end module flare_mmesh_vacuum
