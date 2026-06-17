module flare_jorek
  use iso_fortran_env
  use moose_analysis
  use flare_bfield
  use flare_equi2d
  implicit none
  private


  ! search hint for coordinate transformation ..................................
  type :: search_hint
     real(real64) :: r = 0.d0, z = 0.d0, s = 0.d0, t = 0.d0
     integer :: id = 0
  end type
  ! search_hint ................................................................



  ! mesh for coordinate transformation (r,z) <-> (id,s,t) ......................
  type, public :: jorek_mesh
     real(real64), allocatable :: r(:,:,:), z(:,:,:), s(:,:,:)
     integer, allocatable :: n(:,:)
     integer :: nelem

     ! auxiliary mesh for backward transform
     real(real64) :: rmin, rmax, zmin, zmax
     real(real64), allocatable :: smap(:,:), tmap(:,:)
     integer, allocatable :: imap(:,:)
     integer :: naux
     type(search_hint), pointer :: p1

     contains
     procedure :: broadcast

     procedure :: eval => forward_transform
     procedure :: backward_transform, J1d
  end type jorek_mesh
  ! jorek_mesh .................................................................



  ! scalar data field on mesh ..................................................
  type, extends(scalar_mfunc2d), public :: jorek_field
     type(jorek_mesh), pointer :: mesh
     real(real64), pointer :: values(:,:,:)

     ! near axis expansion
     real(real64) :: psi0, psir, psiz, r0, z0, dr, dz

     contains
     procedure :: broadcast => field_broadcast
     procedure :: near_axis, out_of_bounds

     ! evaluation at (id,s,t)
     procedure :: local_eval, local_deriv, local_hessian
     ! evaluation at (r,z)
     procedure :: eval, deriv, hessian
  end type jorek_field
  ! jorek_field ................................................................



  ! implementation of axisymmetric equilibrium .................................
  type, extends(equi2d), public :: jorek_equi2d
     real(real64) :: F0

     type(jorek_mesh), pointer :: mesh
     real(real64), pointer :: values(:,:,:,:)
     integer, allocatable :: n(:)
     integer :: nmode

     contains
     procedure :: broadcast => broadcast_jorek_equi2d

     procedure :: FpsiN, Fx, FdF
  end type jorek_equi2d
  ! jorek_equi2d ...............................................................



  ! implementation of perturbation field (toroidal mode) .......................
  type, extends(magnetic_field), public :: jorek_bfield
     type(jorek_field) :: cos_psi, sin_psi
     integer :: k, n
     real(real64) :: amplitude

     contains
     procedure :: aux_broadcast
     procedure :: broadcast => broadcast_jorek_bfield
     procedure :: eval => eval_bfield
     procedure :: jac
  end type jorek_bfield
  ! jorek_bfield ...............................................................



  public :: &
     basis_functions_2d, &
     load_jorek_equi2d, init_jorek_bfield

  contains
  !-----------------------------------------------------------------------------


! jorek_mesh ===================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function new_jorek_mesh(n, backward_transform) result(this)
  integer,          intent(in) :: n
  logical,          intent(in) :: backward_transform
  type(jorek_mesh)             :: this


  allocate (this%r(4,4,n), source=0.d0)
  allocate (this%z(4,4,n), source=0.d0)
  allocate (this%s(4,4,n), source=0.d0)
  allocate (this%n(4,n), source=0)
  this%nelem = n


  if (backward_transform) then
     allocate (this%p1)
  else
     this%naux = 0
     nullify (this%p1)
  endif

  end function new_jorek_mesh
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_backward_transform(this, dirname)
  class(jorek_mesh), intent(inout) :: this
  character(len=*),  intent(in   ) :: dirname

  character(len=len(dirname)+20) :: filename
  logical :: ex


  filename = dirname//"/.backward_transform"
  inquire (file=filename, exist=ex)
  if (ex) then
     call aux_load(this, filename)
  else
     call aux_init(this, 1024, filename)
  endif

  end subroutine init_backward_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init(this, naux, filename)
  use moose_netcdf
  use flare_control
  class(jorek_mesh), intent(inout) :: this
  integer,           intent(in   ) :: naux
  character(len=*),  intent(in   ) :: filename

  type(netcdf_dataset) :: N
  real(real64) :: d, dmin, dr, dz, r, raux, z, zaux, s1, t1
  integer :: i, iaux, id1, istat, jaux, dimid


  this%rmin = minval(this%r(:,1,:))
  this%rmax = maxval(this%r(:,1,:))
  this%zmin = minval(this%z(:,1,:))
  this%zmax = maxval(this%z(:,1,:))
  this%naux = naux
  allocate (this%imap(0:naux-1, 0:naux-1), source=0)
  allocate (this%smap(0:naux-1, 0:naux-1), source=0.5d0)
  allocate (this%tmap(0:naux-1, 0:naux-1), source=0.5d0)


  dr = (this%rmax - this%rmin) / naux
  dz = (this%zmax - this%zmin) / naux
  this%rmin = this%rmin - dr / 100
  this%rmax = this%rmax + dr / 100
  this%zmin = this%zmin - dz / 100
  this%zmax = this%zmax + dz / 100


  ! for each auxiliary cell, find mesh cell with shortest distance
  print 1000
  do iaux=0,naux-1
  do jaux=0,naux-1
     call progress_bar(iaux*naux+jaux, naux**2)
     raux = this%rmin + (iaux + 0.5d0) * dr
     zaux = this%zmin + (jaux + 0.5d0) * dz
     dmin = huge(1.d0)
     do i=1,this%nelem
        r = sum(this%r(:,1,i)) / 4
        z = sum(this%z(:,1,i)) / 4
        d = (r-raux)**2 + (z-zaux)**2
        if (d < dmin) then
           dmin = d
           this%imap(iaux, jaux) = i
        endif
     enddo
  enddo
  enddo
  call finalize_progress_bar()
 1000 format(3x,"- setting up search hints for backward transform")


  ! pass 2: compute coordinate transformation for each auxiliary cell
  do iaux=0,naux-1
  do jaux=0,naux-1
     raux = this%rmin + (iaux + 0.5d0) * dr
     zaux = this%zmin + (jaux + 0.5d0) * dz
     call this%backward_transform(raux, zaux, id1, s1, t1, istat)
     if (istat /= 0) cycle

     this%imap(iaux, jaux) = id1
     this%smap(iaux, jaux) = s1
     this%tmap(iaux, jaux) = t1
  enddo
  enddo


  N = netcdf_create(filename)
  call N%def_dim("naux", naux, dimid)
  call N%def_var("rmin", NF90_DOUBLE)
  call N%def_var("rmax", NF90_DOUBLE)
  call N%def_var("zmin", NF90_DOUBLE)
  call N%def_var("zmax", NF90_DOUBLE)
  call N%def_var("imap", NF90_INT, [dimid, dimid])
  call N%def_var("smap", NF90_DOUBLE, [dimid, dimid])
  call N%def_var("tmap", NF90_DOUBLE, [dimid, dimid])
  call N%enddef()

  call N%put_var("rmin", this%rmin)
  call N%put_var("rmax", this%rmax)
  call N%put_var("zmin", this%zmin)
  call N%put_var("zmax", this%zmax)
  call N%put_var("imap", this%imap)
  call N%put_var("smap", this%smap)
  call N%put_var("tmap", this%tmap)
  call N%close()

  end subroutine aux_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_load(this, filename)
  use moose_netcdf
  class(jorek_mesh), intent(inout) :: this
  character(len=*),  intent(in   ) :: filename

  type(netcdf_dataset) :: N


  N = netcdf_open(filename)
  this%naux = N%dim("naux")
  allocate (this%imap(0:this%naux-1, 0:this%naux-1))
  allocate (this%smap(0:this%naux-1, 0:this%naux-1))
  allocate (this%tmap(0:this%naux-1, 0:this%naux-1))
  call N%get_var("rmin", this%rmin)
  call N%get_var("rmax", this%rmax)
  call N%get_var("zmin", this%zmin)
  call N%get_var("zmax", this%zmax)
  call N%get_var("imap", this%imap)
  call N%get_var("smap", this%smap)
  call N%get_var("tmap", this%tmap)
  call N%close()

  end subroutine aux_load
  !-----------------------------------------------------------------------------


! type-bound procedures
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(jorek_mesh), intent(inout) :: this


  call proc(0)%broadcast(this%nelem)
  call proc(0)%broadcast_allocatable(this%r)
  call proc(0)%broadcast_allocatable(this%z)
  call proc(0)%broadcast_allocatable(this%s)
  call proc(0)%broadcast_allocatable(this%n)
  call proc(0)%broadcast(this%naux)
  if (this%naux > 0) then
     call proc(0)%broadcast(this%rmin)
     call proc(0)%broadcast(this%rmax)
     call proc(0)%broadcast(this%zmin)
     call proc(0)%broadcast(this%zmax)
     call proc(0)%broadcast_allocatable(this%imap)
     call proc(0)%broadcast_allocatable(this%smap)
     call proc(0)%broadcast_allocatable(this%tmap)

     if (rank > 0) allocate (this%p1)
  else
     if (rank > 0) nullify (this%p1)
  endif

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function forward_transform(this, id, s, t) result(x)
  !
  ! forward coordinate transformation: (R,Z)(id,s,t)
  !
  class(jorek_mesh), intent(in) :: this
  integer,           intent(in) :: id
  real(real64),      intent(in) :: s, t
  real(real64)                  :: x(2)

  real(real64) :: H(4,4,0:5)
  integer :: k, l


  H = basis_functions_2d(s, t, 0)
  x(1) = sum(this%r(:,:,id) * this%s(:,:,id) * H(:,:,0))
  x(2) = sum(this%z(:,:,id) * this%s(:,:,id) * H(:,:,0))

  end function forward_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine backward_transform(this, r, z, id, s, t, istat, iterations)
  !
  ! istat =    0: success
  !          1-4: out of bounds at upper/lower radial/poloidal boundary
  !            5: detJ = 0
  !            6: backward transform not initialized
  !           99: (r,z) out of bounds of auxiliary mesh
  class(jorek_mesh), intent(in   ) :: this
  real(real64),      intent(in   ) :: r, z
  integer,           intent(  out) :: id, istat
  real(real64),      intent(  out) :: s, t
  integer,           intent(  out), optional :: iterations

  real(real64), parameter :: eps = 1.d-10
  integer, parameter :: n = 256

  type(search_hint), pointer :: hint
  real(real64) :: ds, dt
  integer :: i, iaux, jaux


  if (.not.associated(this%p1)) then
     istat = 6
     return
  endif


  hint => this%p1
  ! re-use previous results, if applicable
  if (r == hint%r  .and.  z == hint%z) then
     id = hint%id
     s = hint%s
     t = hint%t
     istat = 0
     if (present(iterations)) iterations = 0
  endif


  iaux = int(this%naux * (r - this%rmin) / (this%rmax - this%rmin))
  jaux = int(this%naux * (z - this%zmin) / (this%zmax - this%zmin))
  if (iaux < 0 .or. jaux < 0 .or. iaux >= this%naux .or. jaux >= this%naux) then
     istat = 99
     return
  endif


  id = this%imap(iaux, jaux)
  s = this%smap(iaux, jaux)
  t = this%tmap(iaux, jaux)
  do i=1,n
     ! compute step size (ds,dt) approximation in order to reach (r,z) from current position
     call this%J1d(id, s, t, r, z, ds, dt, istat)
     if (istat /= 0) return
     ! truncate small values
     if (abs(ds) < 1.d-12) ds = 0.d0
     if (abs(dt) < 1.d-12) dt = 0.d0

     ! take step (ds,dt)
     call step(s, t, ds, dt, istat)
     ! across cell surface
     if (istat > 0) then
        id = this%n(istat, id)
        if (id == 0) return ! out of bounds

        select case(istat)
        case(1) ! across lower poloidal surface
           t = 1.d0
        case(2) ! across upper radial surface
           s = 0.d0
        case(3) ! across upper poloidal surface
           t = 0.d0
        case(4) ! across lower radial surface
           s = 1.d0
        end select
     endif

     ! check if last step is within required accuracy of target
     if (abs(ds)+abs(dt) < eps) then
        istat = 0
        hint%r = r
        hint%z = z
        hint%s = s
        hint%t = t
        hint%id = id
        if (present(iterations)) iterations = i
        return
     endif
  enddo
  if (present(iterations)) iterations = 0

  end subroutine backward_transform
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine backward_transform_failed_error(x, procname, istat)
  use moose_error
  real(real64),     intent(in) :: x(3)
  character(len=*), intent(in) :: procname
  integer,          intent(in) :: istat


  print *, "r [m]     = ", x(1)
  print *, "z [m]     = ", x(2)
  print *, "phi [deg] = ", x(3) / pi * 180.d0
  call ERROR("backward transform failed", procname, istat)

  end subroutine backward_transform_failed_error
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine J1d(this, id, s, t, rdest, zdest, ds, dt, istat)
  class(jorek_mesh), intent(in   ) :: this
  integer,           intent(in   ) :: id
  real(real64),      intent(in   ) :: s, t, rdest, zdest
  real(real64),      intent(  out) :: ds, dt
  integer,           intent(  out) :: istat

  real(real64) :: det, dr, dz, H(4,4,0:5), r, z, jac(2,2)


  H = basis_functions_2d(s, t, 1)
  r = sum(this%r(:,:,id) * this%s(:,:,id) * H(:,:,0))
  z = sum(this%z(:,:,id) * this%s(:,:,id) * H(:,:,0))
  jac(1,1) = sum(this%r(:,:,id) * this%s(:,:,id) * H(:,:,1))
  jac(1,2) = sum(this%r(:,:,id) * this%s(:,:,id) * H(:,:,2))
  jac(2,1) = sum(this%z(:,:,id) * this%s(:,:,id) * H(:,:,1))
  jac(2,2) = sum(this%z(:,:,id) * this%s(:,:,id) * H(:,:,2))

  istat = 0
  det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
  if (det == 0.d0) then
     istat = 5
     return
  endif
  dr = r - rdest
  dz = z - zdest
  ds = - (jac(2,2) * dr - jac(1,2) * dz) / det
  dt = - (jac(1,1) * dz - jac(2,1) * dr) / det

  end subroutine J1d
  !-----------------------------------------------------------------------------
! jorek_mesh ===================================================================



! jorek_field ==================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_jorek_field(mesh, values) result(this)
  type(jorek_mesh), target, intent(in) :: mesh
  real(real64),     target, intent(in) :: values(:,:,:)
  type(jorek_field)                    :: this


  this = aux_init_field_ptr(mesh, values)
  call init_scalar_mfunc2d(this, [mesh%rmin, mesh%zmin], [mesh%rmax, mesh%zmax])
  call aux_init_near_axis_expansion(this)

  end function init_jorek_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function aux_init_field_ptr(mesh, values) result(this)
  type(jorek_mesh), target, intent(in) :: mesh
  real(real64),     target, intent(in) :: values(:,:,:)
  type(jorek_field)                    :: this


  this%mesh => mesh
  this%values => values

  end function aux_init_field_ptr
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_near_axis_expansion(this)
  use moose_linalg, only: linregress
  class(jorek_field), intent(inout) :: this

  real(real64), allocatable :: grad_psi(:,:)
  real(real64) :: ar, az, br, bz
  integer :: i, ncore


  ! 1. poloidal resolution in core domain
  do i=1,this%mesh%nelem
     if (this%mesh%n(4,i) /= 0) then
        ncore = i-1
        exit
     endif
  enddo

  ! 2. compute grad psi
  allocate (grad_psi(ncore,2))
  do i=1,ncore
     grad_psi(i,:) = this%local_deriv(i, 1.d0, 0.d0)
  enddo

  ! 3. linear fit for grad psi
  call linregress(this%mesh%r(2,1,1:ncore), grad_psi(:,1), ar, br)
  call linregress(this%mesh%z(2,1,1:ncore), grad_psi(:,2), az, bz)
  this%r0 = -br / ar;   this%psir = ar / 2
  this%z0 = -bz / az;   this%psiz = az / 2
  this%psi0 = (sum(this%values(2,1,1:ncore)) &
            - this%psir * sum((this%mesh%r(2,1,1:ncore)-this%r0)**2) &
            - this%psiz * sum((this%mesh%z(2,1,1:ncore)-this%z0)**2)) / ncore
  this%dr = (maxval(this%mesh%r(2,1,1:ncore)) - minval(this%mesh%r(2,1,1:ncore))) / 2
  this%dz = (maxval(this%mesh%z(2,1,1:ncore)) - minval(this%mesh%z(2,1,1:ncore))) / 2
  deallocate (grad_psi)

  end subroutine aux_init_near_axis_expansion
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine field_broadcast(this)
  use moose_mpi
  class(jorek_field), intent(inout) :: this


  call this%mfunc_broadcast()
  call proc(0)%broadcast(this%psi0)
  call proc(0)%broadcast(this%psir)
  call proc(0)%broadcast(this%psiz)
  call proc(0)%broadcast(this%r0)
  call proc(0)%broadcast(this%z0)
  call proc(0)%broadcast(this%dr)
  call proc(0)%broadcast(this%dz)

  end subroutine field_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function near_axis(this, x)
  class(jorek_field), intent(in) :: this
  real(real64),       intent(in) :: x(2)
  logical                        :: near_axis


  near_axis = .false.
  if (abs(x(1)-this%r0) < this%dr  .and.  abs(x(2)-this%z0) < this%dz) near_axis = .true.
  end function near_axis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function out_of_bounds(this, x)
  class(jorek_field), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  logical                        :: out_of_bounds

  real(real64) :: s, t
  integer :: id, istat


  out_of_bounds = .false.
  if (this%near_axis(x)) return

  call this%mesh%backward_transform(x(1), x(2), id, s, t, istat)
  if (istat /= 0) out_of_bounds = .true.

  end function out_of_bounds
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function local_eval(this, id, s, t) result(val)
  use moose_error
  class(jorek_field), intent(in) :: this
  integer,            intent(in) :: id
  real(real64),       intent(in) :: s, t
  real(real64)                   :: val

  real(real64) :: H(4,4,0:5)


  H = basis_functions_2d(s, t, 0)
  val = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,0))

  end function local_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function local_deriv(this, id, s, t) result(deriv)
  use moose_error
  class(jorek_field), intent(in) :: this
  integer,            intent(in) :: id
  real(real64),       intent(in) :: s, t
  real(real64)                   :: deriv(2)

  real(real64) :: det, dpsi_ds, dpsi_dt, H(4,4,0:5), jac(2,2)


  H = basis_functions_2d(s, t, 1)
  dpsi_ds = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  dpsi_dt = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))

  ! Jacobian
  jac(1,1) = sum(this%mesh%r(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  jac(1,2) = sum(this%mesh%r(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))
  jac(2,1) = sum(this%mesh%z(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  jac(2,2) = sum(this%mesh%z(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))
  det = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)
  if (det == 0.d0) call ERROR("det = 0", "jorek_field%local_deriv")

  deriv(1) = (dpsi_ds * jac(2,2) - dpsi_dt * jac(2,1)) / det
  deriv(2) = (dpsi_dt * jac(1,1) - dpsi_ds * jac(1,2)) / det

  end function local_deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function local_hessian(this, id, s, t) result(hessian)
  use moose_error
  use moose_math, only: xproduct => wedge_product
  class(jorek_field), intent(in) :: this
  integer,            intent(in) :: id
  real(real64),       intent(in) :: s, t
  real(real64)                   :: hessian(2,2)

  real(real64) :: det, det1, d1det(2), grad_psi(2), d1psi(2)
  real(real64) :: d2r(3), d2z(3), d2psi(3), M(2,2)
  real(real64) :: H(4,4,0:5), jac(2,2)
  integer :: k


  ! 1st order (s,t)-derivatives
  H = basis_functions_2d(s, t, 2)
  d1psi(1) = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  d1psi(2) = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))
  jac(1,1) = sum(this%mesh%r(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  jac(2,1) = sum(this%mesh%r(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))
  jac(1,2) = sum(this%mesh%z(:,:,id) * this%mesh%s(:,:,id) * H(:,:,1))
  jac(2,2) = sum(this%mesh%z(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2))
  det = jac(1,1) * jac(2,2) - jac(1,2) * jac(2,1)
  if (det == 0.d0) call ERROR("det = 0", "jorek_field%local_hessian")
  det1 = 1.d0 / det

  ! 1st order (r,z)-derivatives
  grad_psi(1) =  det1 * xproduct(d1psi, jac(:,2))
  grad_psi(2) = -det1 * xproduct(d1psi, jac(:,1))

  ! 2nd order (s,t)-derivatives
  do k=1,3
     d2r(k) = sum(this%mesh%r(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2+k))
     d2z(k) = sum(this%mesh%z(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2+k))
     d2psi(k) = sum(this%values(:,:,id) * this%mesh%s(:,:,id) * H(:,:,2+k))
  enddo

  ! mixed derivatives
  d1det(1) = xproduct(d2r(1:2), jac(:,2)) - xproduct(d2z(1:2), jac(:,1))
  d1det(2) = xproduct(d2r(2:3), jac(:,2)) - xproduct(d2z(2:3), jac(:,1))
  M(1,1) = -d1det(1)*grad_psi(1) - xproduct(d2z(1:2), d1psi) + xproduct(d2psi(1:2), jac(:,2))
  M(2,1) = -d1det(2)*grad_psi(1) - xproduct(d2z(2:3), d1psi) + xproduct(d2psi(2:3), jac(:,2))
  M(1,2) = -d1det(1)*grad_psi(2) + xproduct(d2r(1:2), d1psi) - xproduct(d2psi(1:2), jac(:,1))
  M(2,2) = -d1det(2)*grad_psi(2) + xproduct(d2r(2:3), d1psi) - xproduct(d2psi(2:3), jac(:,1))
  M = M * det1

  ! 2nd order (r,z)-derivatives
  hessian(1,1) =  det1 * xproduct(M(:,1), jac(:,2))
  hessian(1,2) =  det1 * xproduct(M(:,2), jac(:,2))
  hessian(2,2) = -det1 * xproduct(M(:,2), jac(:,1))
  hessian(2,1) =  hessian(1,2)

  end function local_hessian
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(psi)
  class(jorek_field), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: psi

  real(real64) :: s, t
  integer :: id, istat


  if (this%near_axis(x)) then
     psi = this%psi0 + this%psir * (x(1) - this%r0)**2 + this%psiz * (x(2) - this%z0)**2
  else
     call this%mesh%backward_transform(x(1), x(2), id, s, t, istat)
     if (istat /= 0) call backward_transform_failed_error(x, "jorek_field%eval", istat)
     psi = this%local_eval(id, s, t)
  endif

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x) result(dpsi)
  class(jorek_field), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: dpsi(this%ndim)

  real(real64) :: s, t
  integer :: id, istat


  if (this%near_axis(x)) then
     dpsi(1) = 2 * this%psir * (x(1) - this%r0)
     dpsi(2) = 2 * this%psiz * (x(2) - this%z0)
  else
     call this%mesh%backward_transform(x(1), x(2), id, s, t, istat)
     if (istat /= 0) call backward_transform_failed_error(x, "jorek_field%deriv", istat)
     dpsi = this%local_deriv(id, s, t)
  endif

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function hessian(this, x)
  class(jorek_field), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: hessian(this%ndim, this%ndim)

  real(real64) :: s, t
  integer :: id, istat


  if (this%near_axis(x)) then
     hessian(1,1) = 2 * this%psir
     hessian(2,2) = 2 * this%psiz
     hessian(1,2) = 0.d0
     hessian(2,1) = 0.d0
  else
     call this%mesh%backward_transform(x(1), x(2), id, s, t, istat)
     if (istat /= 0) call backward_transform_failed_error(x, "jorek_field%hessian", istat)
     hessian = this%local_hessian(id, s, t)
  endif

  end function hessian
  !-----------------------------------------------------------------------------
! jorek_field ==================================================================



! jorek_equi2d =================================================================
! constructors
  !-----------------------------------------------------------------------------
  function load_jorek_equi2d(filename) result(this)
  use moose_utils, only: dirname
  use moose_netcdf
  character(len=*), intent(in) :: filename
  type(jorek_equi2d)           :: this

  type(jorek_mesh), pointer :: mesh
  type(netcdf_dataset) :: N
  integer :: k, nelem


  N = netcdf_open(filename)

  ! load mesh
  nelem = N%dim("nelem")
  allocate (mesh, source=new_jorek_mesh(nelem, .true.))
  call N%get_var("r", mesh%r)
  call N%get_var("z", mesh%z)
  call N%get_var("s", mesh%s)
  call N%get_var("next", mesh%n)
  this%mesh => mesh
  call init_backward_transform(mesh, dirname(filename))
  call init_magnetic_field(this, mesh%rmin, mesh%rmax, mesh%zmin, mesh%zmax)
  print 1000, mesh%rmin, mesh%rmax, mesh%zmin, mesh%zmax
 1000 format(8x,'Bounding box:               ', &
             'R      = ',f8.3, ' m  ->  ',f8.3,' m',/36x, &
             'Z      = ',f8.3, ' m  ->  ',f8.3,' m')

  ! load toroidal mode numbers
  this%nmode = N%dim("nmode")
  allocate (this%n(this%nmode))
  call N%get_var("n", this%n)

  ! load poloidal flux values
  k = N%dim("k")
  allocate (this%values(4,4,nelem,k))
  call N%get_var("psi", this%values);   this%values = -this%values
  call N%get_var("F0", this%F0);   this%F0 = -this%F0
  call N%close()
  allocate (this%psi, source=init_jorek_field(mesh, this%values(:,:,:,1)))

  ! initialize magnetic axis
  this%delta_Psi = 0.d0
  call this%setup(mesh%r(1,1,1), mesh%z(1,1,1), dirname(filename))
  print 2001, this%Bt_axis
  print 2002, this%delta_Psi
 2001 format(8x,'Toroidal magnetic field:  Bt(axis) = ',f8.3,' T')
 2002 format(8x,'Delta Psi                          = ',e12.5)

  end function load_jorek_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_jorek_equi2d(this)
  use moose_mpi
  class(jorek_equi2d), intent(inout) :: this


  if (rank > 0) allocate (this%mesh)
  call this%mesh%broadcast()

  call proc(0)%broadcast(this%F0)
  call proc(0)%broadcast(this%nmode)
  call proc(0)%broadcast_allocatable(this%n)
  if (rank > 0) allocate (this%values(4,4,this%mesh%nelem,2*this%nmode+1))
  call proc(0)%broadcast(this%values)

  if (rank > 0) allocate (this%psi, source=aux_init_field_ptr(this%mesh, this%values(:,:,:,1)))
  call this%equi2d_broadcast()

  end subroutine broadcast_jorek_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function FpsiN(this, psiN) result(F)
  class(jorek_equi2d), intent(in) :: this
  real(real64),        intent(in) :: psiN
  real(real64)                    :: F


  F = this%F0

  end function FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function Fx(this, r) result(F)
  class(jorek_equi2d), intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: F


  F = this%F0

  end function Fx
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function FdF(this, r)
  class(jorek_equi2d), intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: FdF(0:1)


  FdF(0) = this%F0
  FdF(1) = 0.d0

  end function FdF
  !-----------------------------------------------------------------------------
! jorek_equi2d =================================================================



! jorek_bfield =================================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_jorek_bfield(n, amplitude, E) result(this)
  integer,             intent(in) :: n
  real(real64),        intent(in) :: amplitude
  class(jorek_equi2d), intent(in) :: E
  type(jorek_bfield)              :: this


  this%k = find_index(n)
  this%n = n
  this%amplitude = amplitude
  print *
  print 1000
  if (amplitude == 1.d0) then
     print 1001, n
  else
     print 1002, n, amplitude
  endif
 1000 format(3x,"- Magnetic field from JOREK")
 1001 format(8x,"n = ",i0)
 1002 format(8x,"n = ",i0," amplitude = ",e8.3)

  associate (mesh => E%mesh)
  call init_magnetic_field(this, mesh%rmin, mesh%rmax, mesh%zmin, mesh%zmax, n)
  this%cos_psi = init_jorek_field(mesh, E%values(:,:,:,2*this%k))
  this%sin_psi = init_jorek_field(mesh, E%values(:,:,:,2*this%k+1))
  end associate

  contains
  !.............................................................................
  function find_index(n) result(k)
  use moose_error
  use moose_utils
  integer, intent(in) :: n
  integer             :: k


  do k=1,E%nmode
     if (n == E%n(k)) return
  enddo
  call ERROR("jorek data set does not support toroidal mode number n = "//str(n))

  end function find_index
  !.............................................................................
  end function init_jorek_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine aux_broadcast(this, E, k)
  use moose_mpi
  class(jorek_bfield), intent(inout) :: this
  class(jorek_equi2d), intent(in   ) :: E
  integer,             intent(in   ) :: k


  if (rank > 0) then
     this%cos_psi = aux_init_field_ptr(E%mesh, E%values(:,:,:,2*k))
     this%sin_psi = aux_init_field_ptr(E%mesh, E%values(:,:,:,2*k+1))
  endif

  end subroutine aux_broadcast
  !-----------------------------------------------------------------------------
  subroutine broadcast_jorek_bfield(this)
  use moose_mpi
  class(jorek_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%n)
  call proc(0)%broadcast(this%k)
  call proc(0)%broadcast(this%amplitude)

  call this%cos_psi%broadcast()
  call this%sin_psi%broadcast()

  end subroutine broadcast_jorek_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_bfield(this, x) result(B)
  class(jorek_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: B(this%mdim)

  real(real64) :: cosnphi, sinnphi, dpsi(2)


  cosnphi = cos(this%n * x(3))
  sinnphi = sin(this%n * x(3))
  dpsi = cosnphi * this%cos_psi%deriv(x(1:2)) - sinnphi * this%sin_psi%deriv(x(1:2))
  B(1) = -dpsi(2) / x(1) * this%amplitude
  B(2) =  dpsi(1) / x(1) * this%amplitude
  B(3) = 0.d0

  end function eval_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(jorek_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: jac(this%mdim, this%ndim)


  jac = 0.d0
  ! NOT IMPLEMENTED

  end function jac
  !-----------------------------------------------------------------------------
! jorek_bfield =================================================================



! module procedures:
  !---------------------------------------------------------------------
  pure function basis_functions_1d(x, ideriv) result(h)
  !
  ! evalute 1-D basis functions and derivatives
  !
  ! h(1,0) = B0 + B1     ! B0..B3: Bernstein polynomials
  ! h(2,0) = B1
  ! h(3,0) = B3 + B2
  ! h(4,0) = B2
  !
  ! h(:,1) = d/dx h(:,0)
  ! h(:,2) = d^2/dx^2 h(:,0)
  !
  real(real64), intent(in) :: x
  integer,      intent(in) :: ideriv
  real(real64)             :: h(4, 0:2)

  real(real64) :: q, q2, x2, d1, d2


  x2 = x**2
  q  = 1 - x
  q2 = q**2

  h(1,0) = q2 * (1 + 2*x)
  h(2,0) = 3 * x * q2
  h(3,0) = x2 * (3 - 2*x)
  h(4,0) = 3 * x2 * q
  if (ideriv == 0) return

  ! 1st derivative
  d1 = 3*x - 1
  d2 = d1 - 1
  h(1,1) = - 6 * x * q
  h(2,1) = - 3 * q * d1
  h(3,1) = - h(1,1)
  h(4,1) = - 3 * x * d2
  if (ideriv == 1) return

  ! 2nd derivative
  h(1,2) = 6 * (2*x -1)
  h(2,2) = 6 * d2
  h(3,2) = - h(1,2)
  h(4,2) = - 6 * d1

  end function basis_functions_1d
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  pure function basis_functions_2d(s, t, ideriv) result(H)
  real(real64), intent(in) :: s, t
  integer,      intent(in) :: ideriv
  real(real64)             :: H(4,4,0:5)

  integer, parameter :: &
     ivertex(4) = [0,2,2,0], &
     jvertex(4) = [0,0,2,2], &
     idof(4) = [1,2,1,2], &
     jdof(4) = [1,1,2,2]

  real(real64) :: hs(4,0:2), ht(4,0:2)
  integer :: i, j, vertex, dof


  hs = basis_functions_1d(s, ideriv)
  ht = basis_functions_1d(t, ideriv)
  do vertex=1,4
  do dof=1,4
     i = ivertex(vertex) + idof(dof)
     j = jvertex(vertex) + jdof(dof)
     H(vertex,dof,0) = hs(i,0) * ht(j,0)
  enddo
  enddo
  if (ideriv == 0) return

  ! 1st derivative
  do vertex=1,4
  do dof=1,4
     i = ivertex(vertex) + idof(dof)
     j = jvertex(vertex) + jdof(dof)
     H(vertex,dof,1) = hs(i,1) * ht(j,0)
     H(vertex,dof,2) = hs(i,0) * ht(j,1)
  enddo
  enddo
  if (ideriv == 1) return

  ! 2nd derivative
  do vertex=1,4
  do dof=1,4
     i = ivertex(vertex) + idof(dof)
     j = jvertex(vertex) + jdof(dof)
     H(vertex,dof,3) = hs(i,2) * ht(j,0)
     H(vertex,dof,4) = hs(i,1) * ht(j,1)
     H(vertex,dof,5) = hs(i,0) * ht(j,2)
  enddo
  enddo

  end function basis_functions_2d
  !---------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine step(s, t, ds, dt, istat)
  !
  ! from (s,t) step in directoin (ds,dt)
  !
  ! istat = 0: (s+ds, t+dt)
  ! istat > 0: stop at cell boundary
  !
  real(real64), intent(inout) :: s, t
  real(real64), intent(in   ) :: ds, dt
  integer,      intent(  out) :: istat

  real(real64) :: d, fs, ft, s1, s2, t1, t2


  s1 = s + ds
  t1 = t + dt
  istat = 0

  ! 1. boundary check
  ! lower poloidal boundary
  if (t1 < 0.d0) then
     ft = - t / dt
     s2 = s + ft*ds
     if (s2 >= 0.d0  .and.  s2 <= 1.d0) then
        istat = 1
     endif
  endif
  ! upper radial boundary
  if (s1 > 1.d0) then
     fs = (1-s) / ds
     t2 = t + fs*dt
     if (t2 >= 0.d0  .and.  t2 <= 1.d0) then
        istat = 2
     endif
  endif
  ! upper poloidal boundary
  if (t1 > 1.d0) then
     ft = (1-t) / dt
     s2 = s + ft*ds
     if (s2 >= 0.d0  .and.  s2 <= 1.d0) then
        istat = 3
     endif
  endif
  ! lower radial boundary
  if (s1 < 0.d0) then
     fs = - s / ds
     t2 = t + fs*dt
     if (t2 >= 0.d0  .and.  t2 <= 1.d0) then
        istat = 4
     endif
  endif


  ! 2. step
  select case(istat)
  ! finish inside cell
  case(0)
     s = s1
     t = t1

  ! lower poloidal boundary
  case(1)
     s = s2
     t = 0.d0

  ! upper radial boundary
  case(2)
     s = 1.d0
     t = t2

  ! upper poloidal boundary
  case(3)
     s = s2
     t = 1.d0

  ! lower radial boundary
  case(4)
     s = 0.d0
     t = t2

  end select

  end subroutine step
  !-----------------------------------------------------------------------------

end module flare_jorek
