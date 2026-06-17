! Magnetic Flux Surface mesh
module mfs_mesh
  use iso_fortran_env
  use moose_geometry, only: curve, hypersurf2d
  use mesh_interface
  implicit none
  private

  integer, parameter, public :: &
     RIGHT_HANDED =  1, &
     LEFT_HANDED  = -1, &
     LOWER = -1, &
     UPPER =  1, &
     ASCENT  = 1, &
     ASCENT_LEFT  = 1, &
     ASCENT_RIGHT = 2, &
     DESCENT      = 3, &
     DESCENT_CORE = 3, &
     DESCENT_PFR  = 4, &
     LOWER_TO_UPPER = LOWER, &
     UPPER_TO_LOWER = UPPER

  integer, parameter, public :: &
     RADIAL         = 1, &
     POLOIDAL       = 2, &
     AUTOMATIC      = -1024

  real(real64), parameter :: &
     COMPRESSION_FACTOR = 2.0d0


  ! orthogonal, flux surface aligned mesh with strike point adjustment
  type, public :: t_mfs_mesh
     real(real64), pointer :: mesh(:,:,:) => null()
     integer :: nr, np
     integer :: &
        ir0 = -1, & ! reference radial surface for mesh generation
        ip0 = -1    ! reference poloidal surface for mesh generation

     ! "upstream" poloidal neighbor
     class(t_mfs_mesh), pointer  :: upn => null()
     integer                     :: upn_side

     contains
     procedure :: initialize
     procedure :: connect_to
     procedure :: setup_boundary_nodes
     procedure :: plot_boundary
     procedure :: make_orthogonal_grid
     procedure :: make_interpolated_mesh
     procedure :: make_interpolated_submesh
     procedure :: make_divertor_grid
     procedure :: copy
     procedure :: arclength ! calculate arclength on flux surface ir between nodes ip1 and ip2
     procedure :: push_poloidal
     procedure :: upstream_adjust
     procedure :: upstream_adjust_divertor_leg
     procedure :: save
  end type t_mfs_mesh


  type(hypersurf2d), target, public :: guiding_contour
  type(hypersurf2d), pointer, public :: guiding_contour_L, guiding_contour_R


  contains
!=======================================================================



!=======================================================================
  subroutine initialize(this, nr, np, phi)
  class(t_mfs_mesh)        :: this
  integer, intent(in)      :: nr, np
  real(real64), intent(in) :: phi


  this%nr = nr
  this%np = np
  allocate (this%mesh(0:nr, 0:np, 2), source=0.d0)

  end subroutine initialize
!=======================================================================



!=======================================================================
! Connect this bock to next block (M)
! direction:	boundary direction (radial or poloidal)
! side:		lower-to-upper or upper-to-lower boundary
!=======================================================================
  subroutine connect_to(this, M, direction, side)
  class(t_mfs_mesh), target   :: this
  class(t_mfs_mesh)   :: M
  integer, intent(in) :: direction, side

  real(real64), dimension(:,:,:), pointer :: M1, M2
  integer :: ir, ir1, ir2, ip, ip1, ip2


  M1 => this%mesh
  M2 => M%mesh



  select case(direction)
  ! CONNECT IN RADIAL DIRECTION
  case(RADIAL)
     ! 1. check poloidal resolution of interface
     if (this%np .ne. M%np) then
        write (6, 9000)
        write (6, 9010) this%np, M%np
        stop
     endif

     ! 2. setup radial indices
     select case(side)
     case(LOWER_TO_UPPER)
        ir1 = 0
        ir2 = M%nr
     case(UPPER_TO_LOWER)
        ir1 = this%nr
        ir2 = 0
     case default
        write (6, 9000)
        write (6, 9002) side
        stop
     end select

     ! 3. copy boundary nodes
     do ip=0,this%np
        M2(ir2, ip, :) = M1(ir1, ip, :)
     enddo
     M%ir0 = ir2

  ! CONNECT IN POLOIDAL DIRECTION
  case(POLOIDAL)
     ! 1. check radial resolution of interface
     if (this%nr .ne. M%nr) then
        write (6, 9000)
        write (6, 9011) this%nr, M%nr
        stop
     endif

     ! 2. setup poloidal indices
     select case(side)
     case(LOWER_TO_UPPER)
        ip1 = 0
        ip2 = M%np
     case(UPPER_TO_LOWER)
        ip1 = this%np
        ip2 = 0
     case default
        write (6, 9000)
        write (6, 9002) side
        stop
     end select

     ! 3. copy boundary nodes
     do ir=0,this%nr
        M2(ir, ip2, :) = M1(ir, ip1, :)
     enddo
     M%ip0 = ip2
     M%upn => this
     M%upn_side = side

  case default
     write (6, 9000)
     write (6, 9001) direction
     stop
  end select

 9000 format('error in t_mfs_mesh%connect_to')
 9010 format('invalid poloidal resolution in neighboring blocks: ', i0, 2x, i0)
 9011 format('invalid radial resolution in neighboring blocks: ', i0, 2x, i0)
 9001 format('invalid boundary direction: ', i0)
 9002 format('invalid boundary side id: ', i0)
  end subroutine connect_to
!=======================================================================



!=======================================================================
! Setup boundary nodes in mesh
! boundary_side:	lower or upper boundary
! boundary_type:	radial or poloidal block boundary

! C_boundary:	boundary definition
! spacings:	spacing function for nodes on boundary
!=======================================================================
  subroutine setup_boundary_nodes(this, boundary_side, boundary_type, C_boundary, spacings, i1, i2, debug)
  use moose_quantiles
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: boundary_side, boundary_type
  class(curve),    intent(in) :: C_boundary
  class(qfunc),    intent(in), optional :: spacings
  integer,         intent(in), optional :: i1, i2
  logical,         intent(in), optional :: debug

  real(real64), dimension(:,:,:), pointer :: M
  real(real64) :: tau, x(2)
  integer      :: ir, ip, i11, i22, idir
  logical      :: screen_output


  screen_output = .false.
  if (present(debug)) then
     screen_output = debug
  endif


  M => this%mesh


  select case(boundary_side)
  case(LOWER)
     ir = 0
     ip = 0
  case(UPPER)
     ir = this%nr
     ip = this%np
  case default
     write (6, *) 'error in t_base_mesh%setup_boundary_nodes:'
     write (6, *) 'invalid boundary side ', boundary_side
     stop
  end select


  idir = 1
  select case(boundary_type)
  case(RADIAL)
     this%ir0 = ir
     i11      = 0
     i22      = this%np
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     if (i22 < i11) idir = -1
     do ip=i11,i22,idir
        x = C_boundary%eval(ip-i11, i22-i11, spacings)
        M(ir, ip, :) = x
        if (screen_output) write (6, *) x
     enddo
  case(POLOIDAL)
     this%ip0 = ip
     i11      = 0
     i22      = this%nr
     if (present(i1)) i11 = i1
     if (present(i2)) i22 = i2
     if (i22 < i11) idir = -1
     do ir=i11,i22,idir
        if (present(spacings)) then
           tau = spacings%qquantile(ir-i11, i22-i11)
           if (boundary_side == UPPER) tau = 1.d0 - spacings%qquantile(i22-ir, i22-i11)
        else
           tau = 1.d0 * (ir-i11) / (i22-i11)
           if (boundary_side == UPPER) tau = 1.d0 - 1.d0 * (i22-ir) / (i22-i11)
        endif
        x = C_boundary%eval(C_boundary%a + tau * (C_boundary%b - C_boundary%a))
        M(ir, ip, :) = x
        if (screen_output) write (6, *) x
     enddo
  end select

  end subroutine setup_boundary_nodes
!=======================================================================



!=======================================================================
  subroutine plot_boundary(this, direction, side, prefix, iz)
  class(t_mfs_mesh)            :: this
  integer, intent(in)          :: direction, side, iz
  character(len=*), intent(in) :: prefix

  integer, parameter :: iu = 99

  real(real64), dimension(:,:,:), pointer :: M
  character(len=72) :: filename
  integer :: ir, ip


  M => this%mesh

  write (filename, 1000) prefix, iz, side
 1000 format(a,'_Z',i0,'_side',i0)
  open  (iu, file=filename)
  select case(direction)
  case(RADIAL)
     select case(side)
     case(LOWER)
        ir = 0
     case(UPPER)
        ir = this%nr
     end select
     do ip=0,this%np
        write (iu, *) M(ir, ip, :)
     enddo

  case(POLOIDAL)
     select case(side)
     case(LOWER)
        ip = 0
     case(UPPER)
        ip = this%np
     end select
     do ir=0,this%nr
        write (iu, *) M(ir, ip, :)
     enddo
  end select
  close (iu)

  end subroutine plot_boundary
!=======================================================================



!=======================================================================
! Make (quasi) orthogonal grid
!=======================================================================
  subroutine make_orthogonal_grid(this, rrange, prange, periodic, addX, side, debug)
  use flare_model, only: equi2d
  use flare_rpath2d
  class(t_mfs_mesh)             :: this
  integer, intent(in), optional :: rrange(2), prange(2), addX(2), side
  logical, intent(in), optional :: periodic, debug

  real(real64), dimension(:,:,:), pointer :: M
  type(rpath2d_curve) :: R
  real(real64)  :: psiN(0:this%nr), x(2), psiN_final
  integer       :: ir, ir0, ir1, ir2, ir_final, ip, ip0, ip1, ip2, ipp, direction, ix, ipx
  integer       :: inverse, xdir
  logical       :: screen_output


  screen_output = .false.
  if (present(debug)) then
     screen_output = debug
  endif


  M => this%mesh


  ! set poloidal range
  ip0 = this%ip0
  if (ip0 == 0) then
     ip1 = 1
     ip2 = this%np
     ipp = this%np
  elseif (ip0 == this%np) then
     ip1 = 0
     ip2 = this%np - 1
     ipp = 0
  else
     write (6, 9000)
     write (6, 9001) ip0
     stop
  endif
  if (present(periodic)) then
  if (periodic) then
     ip1 = 1
     ip2 = this%np - 1
  endif
  endif
  if (present(prange)) then
     ip1 = prange(1)
     ip2 = prange(2)
  endif


  ! set radial range
  ir0 = this%ir0
  if (ir0 == 0) then
     ir1        = 1
     ir2        = this%nr
     direction  = ASCENT
     inverse    = DESCENT
  elseif (ir0 == this%nr) then
     ir1        = 0
     ir2        = this%nr - 1
     direction  = DESCENT
     inverse    = ASCENT
  else
     write (6, 9000)
     write (6, 9002) ir0
     stop
  endif
  if (present(rrange)) then
     ir1 = rrange(1)
     ir2 = rrange(2)
  endif


  ! setup reference psiN values
  x         = M(ir0,ip0,:)
  psiN(ir0) = equi2d%psiN(x)
  do ir=ir1,ir2
     x        = M(ir,ip0,:)
     psiN(ir) = equi2d%psiN(x)
     if (screen_output) write (6, *) ir, x, psiN(ir)
  enddo

  select case(direction)
  case(ASCENT)
     ir_final = ir2
  case(DESCENT)
     ir_final = ir1
  end select
  xdir = 1
  if (present(side)) then
     direction = direction + side
     inverse   = inverse   + side
     if (side == 1) xdir = -1
  endif
  psiN_final = psiN(ir_final)


  ! additional X-point to take into account
  ipx = -1
  if (present(addX)) then
     ix  = addX(1)
     ipx = addX(2)

     ! set poloidal index to opposite boundary from reference boundary ip0
     if (ipx == AUTOMATIC) ipx = (1 - ip0 / this%np) * this%np
  endif


  ! set up nodes in poloidal range ip1->ip2 and radial range ir1->ir2
  write (6, 1000) ir0, ir1, ir2, ip0, ip1, ip2
  do ip=ip1,ip2
     if (ip == ipx .and. ix > 0) then
        if (screen_output) write (6, *) ip, ix
        R = rpath2d_curveX(ix, -xdir, RPATH2D_PSIN, psiN_final)
     elseif (ip == ipx .and. ix < 0) then
        if (screen_output) write (6, *) ip, ix
        R = rpath2d_curveX(abs(ix), -xdir, RPATH2D_PSIN, psiN(ir0))
        psiN(ir_final) = R%b ! update psiN(ir_final) due to finite accuracy
     else
        if (screen_output) write (6, *) ip, M(ir0,ip,:)
        R = rpath2d_curve(M(ir0,ip,:), RPATH2D_PSIN, psiN_final)
     endif

     do ir=ir1,ir2
        x = R%eval(psiN(ir))
        M(ir,ip,:) = x
     enddo
  enddo


  ! set up periodic boundaries
  if (present(periodic)) then
  if (periodic) then
     M(:,ipp,:) = M(:,ip0,:)
  endif
  endif


 1000 format(8x,'generate orthogonal mesh: (',i0,': ',i0,' -> ',i0,') x (',i0,': ',i0,' -> ',i0,')')
 9000 format('error in t_mfs_mesh%make_orthogonal_grid')
 9001 format('invalid poloidal reference index ', i0)
 9002 format('invalid radial reference index ', i0)
  end subroutine make_orthogonal_grid
!=======================================================================



!=======================================================================
! Generate interpolated mesh to inner simulation boundary (2 -> ir2)
!=======================================================================
  subroutine make_interpolated_mesh(this, ir2, Sr, C0, C1, psiN1_max, prange)
  use moose_quantiles
  use moose_geometry, only: hypersurf2d
  use flare_model, only: equi2d
  use flare_rpath2d
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: ir2
  class(qfunc),    intent(in) :: Sr
  type(hypersurf2d), intent(in), target :: C0, C1
  real(real64),    intent(in) :: psiN1_max
  integer,         intent(in), optional :: prange(2)

  logical, parameter :: Debug = .false.
  !logical, parameter :: Debug = .true.

  real(real64), dimension(:,:,:), pointer :: M
  type(rpath2d_curve) :: Rtmp
  real(real64)  :: psiN2, eta, x(2)
  integer       :: i, ir1, ip1, ip2, j, nr


  ! set defaults for optional input
  ip1 = 0
  ip2 = this%np
  if (present(prange)) then
     ip1 = prange(1)
     ip2 = prange(2)
  endif


  M  => this%mesh
  nr  = this%nr
  ir1 = 1
  write (6, 1001) ir1+1, ir2-1

  ! sanity check
  psiN2 = equi2d%psiN(M(ir2,0,:)) ! radial location of innermost unperturbed flux surface
  if (psiN2 < psiN1_max) then
     write (6, 9000) ir2, psiN2
     write (6, 9001) psiN1_max
     stop
  endif

  if (Debug) print *, "for poloidal indices ", ip1, " -> ", ip2
  do j=ip1,ip2
     x = M(ir2,j,:)
     if (Debug) print *, j, x
     Rtmp = rpath2d_curve(x, RPATH2D_ARCLENGTH, -huge(1.d0), boundary=C1)

     ! interpolated surfaces
     do i=ir1,ir2-1
        !eta = Sr%qquantile(i-1, nr-1) / Sr%qquantile(ir2-1, nr-1)
        eta = 1.d0 * (i-ir1) / (ir2-ir1)
        x = Rtmp%eval(Rtmp%a + eta * (Rtmp%b-Rtmp%a))

        M(i,j,:) = x
     enddo
     call Rtmp%free()

     ! innermost surface
     x = M(ir1,j,:)
     Rtmp = rpath2d_curve(x, RPATH2D_ARCLENGTH, -huge(1.d0), boundary=C0)
     M(0,j,:) = Rtmp%eval(Rtmp%a)
     !write (93, *) M(0,j,:)
  enddo
  if (Debug) print *, "... done"

 1001 format(8x,'interpolating from inner boundary to 1st unperturbed flux surface: ', &
             i0, ' -> ', i0)
 9000 format('error: last unperturbed flux surface at radial index ', i0, ' is at psiN = ', &
             f0.3, ' but it must be completely outside of inner simulation boundary!'// &
             'try using a larger n_interpolate!')
 9001 format('outer most point on inner simulation boundary is at psiN = ', f0.3)
  end subroutine make_interpolated_mesh
!=======================================================================



!=======================================================================
  subroutine make_interpolated_submesh(this, rrange, prange)
  class(t_mfs_mesh)           :: this
  integer,         intent(in) :: rrange(2), prange(2)

  real(real64), dimension(:,:,:), pointer :: M
  real(real64) :: x(2), s, x1(2), x2(2)
  integer      :: ir, ip


  M => this%mesh

  do ir=rrange(1),rrange(2)
     x1 = M(ir,prange(1),:)
     x2 = M(ir,prange(2),:)
     do ip=prange(1)+1,prange(2)-1
        s = 1.d0 * (ip-prange(1)) / (prange(2)-prange(1))
        x = x1 + s * (x2-x1)
        M(ir,ip,:) = x
     enddo
  enddo

  end subroutine make_interpolated_submesh
!=======================================================================



!=======================================================================
! Generate nodes in the base plane so that field lines connect to the
! strike point x0
!
! Z:     toroidal discretization used for strike point adjustment
!=======================================================================
  subroutine align_strike_points(x0, Z, M)
  use moose_error
  use moose_math, only: pi, pi2
  use flare_mmesh_parameters, only: toroidal_discretization
  use flare_fieldline
  real(real64), intent(in)  :: x0(2)
  type(toroidal_discretization), intent(in)  :: Z
  real(real64), intent(out) :: M(0:Z%nt, 2)

  type(fdriver) :: F
  real(real64)      :: ts, y0(3), y(3), Dphi
  integer           :: idir, it, it0, it_sub, it_end, ierr


  ! set parameters for field line tracing
  ! (taken from subroutine trace_nodes)
  ts = pi2 / 3600.d0



  ! set initial point
  y0(1:2)    = x0
  y0(3)      = Z%phi(Z%it_base) / 180.d0 * pi
  it0        = Z%it_base
  M(it0,1:2) = y0(1:2)

  F = fdriver(stop_at_boundary=.false.)
  do idir=-1,1,2
     call F%reset()

     y = y0
     ! trace from base location to zone boundaries
     it_end = 0
     if (idir > 0) it_end = Z%nt
     do it=Z%it_base+idir,it_end, idir
        ierr = F%evolve3(y, Z%phi(it) / 180.d0 * pi)
        if (ierr .ne. 0) then
           write (6, 9000) ierr
           write (6, *) 'reference point: ', y0
           write (6, *) 'it, idir = ', it, idir
           write (6, *) 'it_end   = ', it_end
           write (6, *) 'it_base  = ', Z%it_base
           write (6, *) 'Dphi[deg]= ', Dphi / pi * 180.d0
           write (6, *) 'present location: ', y
           stop
        endif

        M(it,1:2) = y(1:2)
     enddo
  enddo
 9000 format('error in subroutine align_strike_point: ',//, &
             't_fieldline%trace_Dphi returned ierr = ', i0)
  end subroutine align_strike_points
!=======================================================================



!=======================================================================
! Generate grid in divertor legs
!
! TODO: upstream = 0, downstream = np and flip orientation at the end?
!
! Rside: location of seed mesh
! npA_range:   cell range for quasi-orthogonal grid (upstream)
!=======================================================================
!  subroutine make_divertor_grid(this, R, Rside, Sr, P, Pside, Sp, Z, ierr)
  subroutine make_divertor_grid(this, Rside, ip0, npA_range, Sp, U, Sr, Z, ir_skip, ierr)
  use moose_error
  use moose_math,     only: pi
  use moose_analysis, only: ufunc
  use moose_quantiles
  use moose_geometry, only: contour, interp_curve, hypersurf2d
  use flare_model,    only: equi2d
  use flare_fluxsurf2d
  use flare_mmesh_parameters, only: toroidal_discretization, npSP_subres, npSP_extend, auto_adjust_mesh
  class(t_mfs_mesh)            :: this
!  type(t_curve),   intent(in)  :: R, P
!  integer,         intent(in)  :: Rside, Pside
!  type(t_spacing), intent(in)  :: Sr, Sp
  integer,         intent(in)  :: Rside, ip0, npA_range(2), ir_skip
  class(qfunc),    intent(in)  :: Sp, Sr
  class(ufunc),    intent(in)  :: U
  type(toroidal_discretization),    intent(in)  :: Z
  integer,         intent(out) :: ierr

  integer, parameter :: nskip = 0, nguard = 1, iu_err = 66
  logical, parameter :: Debug = .false.

  type(contour) :: Ftmp
  type(interp_curve) :: C, F
  type(toroidal_discretization) :: TSP
  type(poly_qfunc) :: S

  real(real64), dimension(:,:,:), pointer :: M
  real(real64), dimension(:,:), allocatable :: MSP
  real(real64)  :: psiN(0:this%nr), x(2), psiN_final, tau, xi, xi0, L0, L, dphi, a
  real(real64)  :: Ladjust, Lu1, Lu2, Lu3
  integer       :: ir, ir0, ir1, ir2, ip, ips, ipt, dir, ioffset, np_SP
  integer       :: it, its, it_start, it_end, dirT, downstream, nsub, np_skip
  integer       :: ipu, nextend

  type(hypersurf2d), pointer :: guiding_contour


  ! check intersection of R with guiding surface
!  if (R%intersect_curve(C_guide, x, tau)) then
!     L    = R%length() * tau
!     write (6, 9000) L
!     ierr = 1
!     return
!  endif


  ! set up effective resolution for strike point area
  nsub  = npSP_subres
  nextend = npSP_extend * nsub;   if (nextend < 0) nextend = Z%nt / 2 * nsub
  np_SP = Z%nt * nsub / (nskip+1) + nguard + nextend


  ! check resolution
  np_skip = ip0;  if (Rside == UPPER) np_skip = this%np - ip0
  if (np_SP > this%np-np_skip) call ERROR("poloidal grid resolution too small")


  ! setup poloidal boundary with reference nodes for flux surfaces
  !call this%setup_boundary_nodes(POLOIDAL, Rside, R, Sr)
  select case(Rside)
  case(LOWER)
     dir        = RIGHT_HANDED
     downstream = 1
     ips        = this%np - np_SP
     ipt        = this%np
     ipu        = 0
     guiding_contour => guiding_contour_L
  case(UPPER)
     dir        = LEFT_HANDED
     downstream = -1
     ips        = np_SP
     ipt        = 0
     ipu        = this%np
     guiding_contour => guiding_contour_R
     ! ip_ds    = 0
  end select


  ! find downstream direction
  dirT       = dir * equi2d%Bp_sign * equi2d%Bt_sign
  it_start   = np_SP
  it_end     = 0
  if (dirT > 0) then
     it_start = 0
     it_end   = np_SP
  endif


  ! initialize toroidal discretization for strike point adjustment
  TSP%nt = np_SP
  allocate (TSP%phi(0:np_SP))

  ! extend target aligned mesh
  ioffset = nextend; if (dirT > 0) ioffset = nguard
  if (dirT < 0) then
     ! add guard cell(s) beyond target
     dphi = Z%phi(Z%nt) - Z%phi(Z%nt-1)
     do its=1,nguard
        TSP%phi(np_SP-nguard + its) = Z%phi(Z%nt) + its*dphi
     enddo

     ! extend alignment for mapping to next zone
     ! NOTE: this should be determined by the toroidal resolution in the next zone
     dphi = (Z%phi(1) - Z%phi(0)) / nsub
     do its=1,nextend
        TSP%phi(nextend - its) = Z%phi(0) - its*dphi
     enddo
  else
     ! add guard cell(s) beyond target
     dphi = Z%phi(1) - Z%phi(0)
     do its=1,nguard
        TSP%phi(nguard - its) = Z%phi(0) - its*dphi
     enddo

     ! extend alignment for mapping to next zone
     dphi = (Z%phi(Z%nt) - Z%phi(Z%nt-1)) / nsub
     do its=1,nextend
        TSP%phi(np_SP-nextend + its) = Z%phi(Z%nt) + its*dphi
     enddo
  endif

  ! set up toroidal discretization for strike point adjustment
  TSP%it_base = Z%it_base * nsub / (nskip+1)  +  ioffset
  do it=0,Z%nt
     TSP%phi(it*nsub + ioffset) = Z%phi(it)
  enddo
  ! add sub-resolution
  do it=0,Z%nt-1
     dphi = Z%phi(it+1) - Z%phi(it)

     do its=1,nsub-1
        TSP%phi(it*nsub + its + ioffset) = Z%phi(it) + 1.d0 * its/nsub * dphi
     enddo
  enddo


  M => this%mesh
  allocate (MSP(0:TSP%nt, 2))
  ! separatrix leg
  ir = this%ir0
  if (ir >= 0  .and.  ir /= ir_skip) then
     ! separatrix strike point is already known from setup_boundary_nodes
     C = interp_curve(M(ir, min(ipu,ipt):max(ipu,ipt), :), reverse= ipu > ipt)
     x = C%eval(C%b)
     L0 = C%b - C%a

     ! generate nodes from which field lines connect to strike point x
     call align_strike_points(x, TSP, MSP)


     ! setup downstream strike point nodes
     do it=it_start,it_end,dirT
        ip = ips + abs(it - it_end)*dir
        M(ir,ip,:) = MSP(it,:)
        if (Debug) write (88, *) M(ir,ip,:)
     enddo


     ! length on flux surface for strike point mesh
     L  = 0.d0
     do it=TSP%it_base+dirT,it_end,dirT
        L = L + sqrt(sum(  (MSP(it-dirT,:)-MSP(it,:))**2))
     enddo

     ! aligned strike point mesh extends beyond upstream reference location?
     ! -> adjust upstream location by L-L0
     Ladjust = L - L0
     if (Ladjust > 0.d0) then
        write (6, *) "error: separatrix strike point mesh extends beyond upstream reference!"
        stop
     endif


     ! make suggestion for poloidal spacing
     a = (abs(ipu-ips)*sqrt(sum( (M(ir,ips,:)-M(ir,ips+dir,:))**2 ))/(L0-L) - 1.d0) / (1.d0 - 1.d0/abs(ipu-ips))
     S = quadratic_qfunc(a)


     ! resample separatrix leg from X-point to first node of strike point mesh
     do ip=ipu+dir,ips-dir,dir
        tau = 1.d0 * (ip - ipu) / (ips - ipu)
        xi  = S%eval(tau) * (L0-L) / L0

        x = C%eval(C%a + xi * (C%b-C%a))
        M(ir,ip,:) = x
     enddo
     call S%free()
  endif


     ! generate quasi-orthogonal mesh on upstream end
     if (npA_range(2) >= npA_range(1)) then
        call this%make_orthogonal_grid(prange=npA_range)

        ! account for upstream adjustment
        if (Debug) print *, "UPSTREAM_ADJUST_DIVERTOR_LEG"
        call this%upstream_adjust_divertor_leg(U, Sr, npA_range, ips)
     endif



  ! all other flux surfaces
  if (Debug) print *, "DOWNSTREAM MESH"
  do ir=0,this%nr
     if (ir == ir_skip) cycle

     ! separatrix leg (already taken care of)
     if (ir == this%ir0) cycle

     ! divertor leg of flux surface
     x = M(ir,ip0,:)

     ! generate flux surface from "upstream" location x to target
     Ftmp = fluxsurf2d_contour(x, -dir*equi2d%Bp_sign, ierr, guiding_contour)
     ! A. successfull trace of flux surface to target
     if (ierr == -3) then
        F = Ftmp%interp()
!        select case(dir)
!        case(RIGHT_HANDED)
!           x = F%x(F%n_seg, :)
!        case(LEFT_HANDED)
!           x = F%x(0,:)
!        end select
        x  = F%eval(F%b)
        !write (99, *) x
        L0 = F%b - F%a
        !write (71, *) F%eval(F%a), ir
        !write (72, *) F%eval(F%b), ir

     ! B. flux surface out of bounds, trace backwards to target and adjust upstream mesh
     elseif (ierr == -1) then
        write (6, *) 'x0 = ', x
        call Ftmp%free()
        ! adjust "upstream" location
        ! 1. trace back to boundary
        Ftmp = fluxsurf2d_contour(x, -dir*equi2d%Bp_sign, ierr, guiding_contour)
        F = Ftmp%interp()

        ! point on boundary
        x  = F%eval(F%b)
        L0 = -F%b + F%a

     ! C. UNKOWN situation
     else
        call Ftmp%savetxt("FLUXSURF2D_CONTOUR")
        print *, "ierr = ", ierr
        print *, "flux surface contour is stored in FLUXSURF2D_CONTOUR"
        call ERROR("unexpected situation on divertor target")
     endif
     !if (Debug) call F%plot(filename='F.plt', append=.true.)
     call Ftmp%free()


     ! generate nodes from which field lines connect to strike point x
     call align_strike_points(x, TSP, MSP)


     ! setup downstream strike point nodes
     do it=it_start,it_end,dirT
        ip = ips + abs(it - it_end)*dir
        M(ir,ip,:) = MSP(it,:)
        if (Debug) write (88, *) M(ir,ip,:)
     enddo


     ! length on flux surface for strike point mesh
     L  = 0.d0
     do it=TSP%it_base+dirT,it_end,dirT
        L = L + sqrt(sum(  (MSP(it-dirT,:)-MSP(it,:))**2))
        !write (92, *) MSP(it, :)
     enddo

     ! aligned strike point mesh extends beyond upstream reference location?
     ! -> adjust upstream location by L-L0
     Ladjust = L - L0
     if ((Ladjust > 0.d0)  .and.  .not. auto_adjust_mesh) then
        write (6, 9020) ir
        open  (iu_err, file='error_strike_point_mesh1.plt')
        write (iu_err, *) M(ir,ip0,:)
        close (iu_err)
        open  (iu_err, file='error_strike_point_mesh2.plt')
        do it=0,TSP%nt
           write (iu_err, *) MSP(it,:)
        enddo
        close (iu_err)
        open  (iu_err, file='error_upstream_nodes.plt')
        do ir1=0,this%nr
           write (iu_err, *) M(ir1,ip0,:)
        enddo
        close (iu_err)
        stop
     endif
     Lu1 = this%arclength(ir,ip0,ipu)
     Lu2 = 2*L - (L0+Lu1)
     if (Ladjust-Lu1 > -L) then
        if (.not.associated(this%upn)) then
           write (6, *) 'error: pointer to upstream mesh not associated!'
           stop
        endif
        call adjust_upstream(L)

     else
     ! interpolate nodes on flux surface between upstream orthogonal grid nodes
     ! and downstream strike point nodes
        a = (abs(ip0-ips)*sqrt(sum( (M(ir,ips,:)-M(ir,ips+dir,:))**2 ))/(L0-L) - 1.d0) / (1.d0 - 1.d0/abs(ip0-ips))
        S = quadratic_qfunc(a)
        do ip=ip0+dir,ips-dir,dir
           tau = 1.d0 * (ip - ip0) / (ips - ip0)
           xi  = S%eval(tau) * (L0-L) / L0

           x = F%eval(F%a + xi * (F%b-F%a))
           M(ir,ip,:) = x
        enddo
     endif
  enddo

  ierr = 0
  deallocate (MSP)
 9000 format('ERROR: reference path for radial discretization crosses guiding surface!', &
             'intersection at L = ', f0.3)
 9020 format('ERROR: aligned strike point mesh extends beyond upstream reference point ', &
             'at radial index ', i0//,&
             'see error_strike_point_mesh.plt')
  contains
  !---------------------------------------------------------------------
  subroutine adjust_upstream(L)
  real(real64), intent(in) :: L

  real(real64) :: s
  type(fluxsurf2d) :: F


  ! adjust divertor grid
  x = M(ir,ips,:)
  F = fluxsurf2d(x, -dir*equi2d%Bp_sign, param=FLUXSURF2D_PARAM_ARCLENGTH, boundary=guiding_contour)

  do ip=ips-dir,ipu,-dir
     s = 1.d0 * (ip-ips) / (ipu-ips) * L
     call F%out_of_bounds_check(s, "adjust_upstream(L)")
     x = F%eval(s)
     M(ir,ip,:) = x
  enddo


  ! adjust upstream grid
  call this%upn%push_poloidal(ir, dir, Lu2, Lu2*2)

  end subroutine adjust_upstream
  !---------------------------------------------------------------------
  end subroutine make_divertor_grid
!=======================================================================



!=======================================================================
! copy mesh M onto this mesh at node index ir0, ip0
!=======================================================================
  subroutine copy(this, ir0, ip0, M)
  class(t_mfs_mesh)            :: this
  integer,          intent(in) :: ir0, ip0
  type(t_mfs_mesh), intent(in) :: M

  integer :: ir, ip


  ! check boundaries
  if (ir0 < 0  .or.  ir0 > this%nr  .or.  ip0 < 0  .or.  ip0 > this%np) then
     write (6, 9000);  write (6, 9001) ir0, ip0, this%nr, this%np;  stop
  endif

  ! check size
  if (ir0+M%nr > this%nr  .or.  ip0+M%np > this%np) then
     write (6, 9000);  write (6, 9002) ir0+M%nr, ip0+M%np, this%nr, this%np;  stop
  endif

  ! copy mesh
  do ir=0,M%nr
  do ip=0,M%np
     this%mesh(ir0+ir, ip0+ip, :) = M%mesh(ir, ip, :)
  enddo
  enddo

 9000 format('error in t_mfs_mesh%copy:')
 9001 format('initial node is outside of mesh!'//'ir0, ip0 = ',i0,', ',i0// &
             'nr,  np  = ',i0,', ',i0)
 9002 format('upper node is outside of mesh!'//'ir,  ip  = ',i0,', ',i0// &
             'nr,  np  = ',i0,', ',i0)
  end subroutine copy
!=======================================================================



!=======================================================================
! calculate arclength on grid/flux surface ir between nodes ip1 and ip2
!=======================================================================
  function arclength(this, ir, ip1, ip2) result(l)
  class(t_mfs_mesh)   :: this
  integer, intent(in) :: ir, ip1, ip2
  real(real64)        :: l

  real(real64), dimension(:,:,:), pointer :: M
  integer :: ip


  M => this%mesh

  ! check inputs
  call irange_check(ir,  0, this%nr, 'ir', 't_mfs_mesh%arclength')
  call irange_check(ip1, 0, this%np, 'ip', 't_mfs_mesh%arclength')
  call irange_check(ip2, 0, this%np, 'ip', 't_mfs_mesh%arclength')


  l = 0.d0
  do ip=min(ip1,ip2)+1,max(ip1,ip2)
     l = l + sqrt(sum( (M(ir,ip,:)-M(ir,ip-1,:))**2 ))
  enddo

  end function arclength
!=======================================================================



!=======================================================================
  subroutine push_poloidal(this, ir, side, L0, Lpush)
  use moose_error
  use moose_geometry, only: interp_curve
  use flare_fluxsurf2d
  class(t_mfs_mesh)        :: this
  integer,      intent(in) :: ir, side
  real(real64), intent(in) :: L0, Lpush

  type(interp_curve) :: C
  real(real64), dimension(:,:,:), pointer :: M
  real(real64), allocatable :: w(:)
  real(real64) :: L, dl, x(2), s
  integer :: ip, ip1, ip2, ipdir, ierr


  select case(side)
  case(1)
     ip1 = this%np
     ip2 = 0
  case(-1)
     ip1 = 0
     ip2 = this%np
  case default
     write (6, *) 'error: invalid argument side = ', side, '!'
     stop
  end select
  ipdir = -side


  ! find index ip2 for Lpush
  allocate (w(0:this%np), source=0.d0)
  L = 0.d0
  M => this%mesh
  do ip=ip1+ipdir,ip2,ipdir
     dl    = sqrt(sum( (M(ir,ip,:)-M(ir,ip-ipdir,:))**2 ))
     w(ip) = w(ip-ipdir) + dl
     L     = L + dl

     if (L > Lpush) exit
  enddo
  ip2 = ip
  w   = w / Lpush


  ! approximate flux surface from interpolation along mesh nodes
  C = interp_curve(M(ir,:,:), reverse=ip1>ip2)


  do ip=ip1,ip2-ipdir,ipdir
     s = (L0 + w(ip)*(Lpush-L0))
     call C%out_of_bounds_check(s, "push_poloidal")
     x = C%eval(s)
     M(ir,ip,:) = x
  enddo

  deallocate (w)
  end subroutine push_poloidal
!=======================================================================


!=======================================================================
  subroutine upstream_adjust(this, side, U, V, Sr)
  use moose_analysis
  class(t_mfs_mesh)         :: this
  integer,       intent(in) :: side
  type(interp), intent(in) :: U, V
  class(qfunc), intent(in) :: Sr

  real(real64) :: r, x(2), L0, Lpush, w
  integer      :: ir


  do ir=0,this%nr
     w = Sr%qquantile(ir, this%nr)
     L0 = U%eval(U%a + w * (U%b - U%a))
     if (L0 <= 0) cycle

     Lpush = V%eval(V%a + w * (V%b - V%a))
     call this%push_poloidal(ir, side, L0, Lpush)
  enddo

  end subroutine upstream_adjust
!=======================================================================


!=======================================================================
  subroutine upstream_adjust_divertor_leg(this, U, Sr, np_range, ips)
  use moose_analysis, only: ufunc, qfunc
  use flare_model, only: equi2d
  use flare_fluxsurf2d
  class(t_mfs_mesh)         :: this
  class(ufunc), intent(in) :: U
  class(qfunc), intent(in) :: Sr
  integer,       intent(in) :: np_range(2), ips

  type(fluxsurf2d)   :: F
  real(real64), dimension(:,:,:), pointer :: M
  real(real64) :: s, x(2), dl, w(0:this%np), L
  integer      :: idir, ir, ir0, ir1, ir2, ip, ip1, ip2, dp, ierr


  M => this%mesh


  if (this%ir0 == 0) then
     idir = 1
     ir1 = 1
     ir2 = this%nr
  else
     idir = -1
     ir1 = this%nr - 1
     ir2 = 0
  endif


  ! find index ir0 with U > 0
  do ir0=ir1,ir2,idir
     s = U%eval(U%a + Sr%qquantile(ir0, this%nr) * (U%b - U%a))
     if (s > 0.d0) exit
     if (ir0 == ir2) return
  enddo


  ! set weights from last flux surface
  w = 0.d0
  ip1 = this%ip0
  if (ip1 == 0) then
     ip2 = maxval(np_range)
     dp  = 1
  else
     ip2 = minval(np_range)
     dp  = -1
  endif
  do ip=ip1+dp,ip2,dp
     dl    = sqrt(sum( (M(ir0-idir,ip,:)-M(ir0-idir,ip-dp,:))**2 ))
     w(ip) = w(ip-dp) + dl
     write (6, *) ip, M(ir0-idir,ip,:)-M(ir0-idir,ip-dp,:)
  enddo
  !write (6, *) 'DEBUG: ', ir0, ip1, ip2, L, w(ip1:ip2)

  ! reference length on separatrix leg
  x = M(ir0-idir,ip1,:)
  F = fluxsurf2d(x, -dp*equi2d%Bp_sign, param=FLUXSURF2D_PARAM_ARCLENGTH, boundary=guiding_contour)
  L = F%b

  ! move nodes
  do ir=ir0,ir2,idir
     x = M(ir,ip1,:)
     F = fluxsurf2d(x, -dp*equi2d%Bp_sign, param=FLUXSURF2D_PARAM_ARCLENGTH, boundary=guiding_contour)

     do ip=ip1+dp,ip2,dp
        s = w(ip) * F%b / L
        x = F%eval(s)
        M(ir,ip,:) = x
     enddo
  enddo


  end subroutine upstream_adjust_divertor_leg
!=======================================================================


!=======================================================================
subroutine irange_check(i, i1, i2, label, procedure_name)
  integer, intent(in) :: i, i1, i2
  character(len=*), intent(in) :: label, procedure_name


  if (i < i1  .or.  i > i2) then
     write (6, 9000) trim(label), i, i1, trim(label), i2, procedure_name
     stop
  endif

 9000 format('error: ',a,' = ',i0,' when ',i0,' <= ',a,' <= ',i0,' required in ',a,'!')

end subroutine irange_check
!=======================================================================


!=======================================================================
  subroutine save(this, basename, iz)
  use moose_grids, only: qmesh
  class(t_mfs_mesh), intent(in) :: this
  character(len=*),  intent(in) :: basename
  integer,           intent(in) :: iz

  type(qmesh) :: Q
  character(len=128) :: filename


  write (filename, 1000) basename, iz
 1000 format(a,i0,".dat")


  Q = qmesh(this%mesh(:,:,1), this%mesh(:,:,2))
  call Q%savetxt(filename)
  call Q%free()

  end subroutine save
!=======================================================================

end module mfs_mesh
