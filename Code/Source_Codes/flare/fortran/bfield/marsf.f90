module flare_marsf
  use iso_fortran_env
  use moose_algorithms, only: binary_search_L
  use moose_math,       only: pi, pi2
  use flare_control,  only: report
  use flare_bfield
  implicit none
  private


  integer, parameter :: B_FIELD = 1, A_FIELD = 2, X_FIELD = 3

  integer, parameter :: dp = selected_real_kind(12)


  ! coordinate transformation (s, chi) -> (R, Z) (SCHIMESH file) ...............
  type, public :: marsf_schimesh
     real(real64), dimension(:,:), allocatable :: r, z, s, chi
     real(real64) :: raxis, zaxis
     integer :: nr, nz

     contains
     procedure :: broadcast => schimesh_broadcast
     procedure :: free      => schimesh_free

     procedure :: find_mesh_index
     procedure :: forward_transform
  end type marsf_schimesh


  interface marsf_schimesh
     procedure :: load_schimesh
  end interface marsf_schimesh
  ! marsf_schimesh .............................................................



  ! MARS-F data set (BPLASMA file) .............................................
  type, public :: marsf_bplasma
     real(real64), dimension(:), allocatable :: cs, csm, q
     complex(kind=dp), dimension(:,:), allocatable :: rmi, zmi, rmm, zmm
     complex(kind=dp), dimension(:,:), allocatable :: b1, b2, b3, x1

     real(real64) :: rtmp(2)
     integer :: nn, ns, nsp, m1, m2, mmaxe, mmaxp, itmp

     contains
     procedure :: broadcast => bplasma_broadcast
     procedure :: free      => bplasma_free

     procedure :: save

     procedure :: find_s_index
     procedure :: eval_geometry, eval_x1
  end type marsf_bplasma


  interface marsf_bplasma
     procedure :: load_bplasma
  end interface marsf_bplasma
  ! marsf_bplasma ..............................................................



  ! implementation of magnetic fields from MARS-F data set .....................
  type, extends(magnetic_field), public :: marsf_bfield
     type(marsf_schimesh), pointer :: schimesh
     type(marsf_bplasma),  pointer :: bplasma

     contains
     procedure :: broadcast => broadcast_bfield
     procedure :: free      => free_bfield

     procedure :: eval
     procedure :: jac
  end type marsf_bfield


  interface marsf_bfield
     procedure :: init_bfield
     procedure :: load_bfield
  end interface marsf_bfield
  ! marsf_bfield ...............................................................



  public :: &
     merge_bplasma

  contains
  !-----------------------------------------------------------------------------


! type marsf_schimesh ==========================================================
! contructors:
  !-----------------------------------------------------------------------------
  function load_schimesh(filename) result(this)
  character(len=*), intent(in) :: filename
  type(marsf_schimesh)         :: this

  integer :: ir, iz, iu, nr, nz


  open  (newunit=iu, file=filename, action='read', form='formatted')
  read  (iu, *) nr, nz, this%raxis, this%zaxis

  this%nr = nr
  this%nz = nz
  allocate (this%r(nr,nz), this%z(nr,nz), this%s(nr,nz), this%chi(nr,nz))
  do ir=1,nr
  do iz=1,nz
     read  (iu, *) this%r(ir,iz), this%z(ir,iz), this%s(ir,iz), this%chi(ir,iz)
  enddo
  enddo
  close (iu)

  end function load_schimesh
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine schimesh_broadcast(this)
  use moose_mpi
  class(marsf_schimesh), intent(inout) :: this


  call proc(0)%broadcast(this%nr)
  call proc(0)%broadcast(this%nz)
  call proc(0)%broadcast(this%raxis)
  call proc(0)%broadcast(this%zaxis)

  call proc(0)%broadcast_allocatable(this%r)
  call proc(0)%broadcast_allocatable(this%z)
  call proc(0)%broadcast_allocatable(this%s)
  call proc(0)%broadcast_allocatable(this%chi)

  end subroutine schimesh_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine schimesh_free(this)
  class(marsf_schimesh), intent(inout) :: this


  deallocate (this%r, this%z, this%s, this%chi)

  end subroutine schimesh_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine find_mesh_index(this, r, z, ir, iz, ierr)
  !
  ! For given (r, z) find corresponding mesh index (ir, iz).
  ! ierr =
  !        0:    successful
  !        1:    error finding radial index ir
  !        2:    error finding vertical index iz
  !
  class(marsf_schimesh), intent(in)  :: this
  real(real64),          intent(in)  :: r, z
  integer,               intent(out) :: ir, iz, ierr


  ierr = 0
  ! scan radial direction
  if (r == this%r(1,1)) then
     ir = 1
  else
     ir = binary_search_L(this%r(:,1), r)
  endif
  if (ir < 1  .or.  ir >= this%nr) then
     ierr = 1
     return
  endif


  ! scan vertical direction
  if (z == this%z(1,1)) then
     iz = 1
  else
     iz = binary_search_L(this%z(1,:), z)
  endif
  if (iz < 1  .or.  iz >= this%nz) then
     ierr = 2
     return
  endif

  end subroutine find_mesh_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine forward_transform(this, r, z, s, chi, ierr)
  !
  ! For given (r, z) find corresponding (s, chi)
  !
  class(marsf_schimesh), intent(in)  :: this
  real(real64),          intent(in)  :: r, z
  real(real64),          intent(out) :: s, chi
  integer,               intent(out) :: ierr

  real(real64)   :: r1, r2, z1, z2, chi00, chi01, chi10, chi11
  double complex :: ctmp
  integer :: ir, iz


  ! find mesh index for (r, z)
  call this%find_mesh_index(r, z, ir, iz, ierr)
  if (ierr /= 0) return


  ! set up bounding box of cell
  r1 = this%r(ir,   1)
  r2 = this%r(ir+1, 1)
  z1 = this%z(1, iz  )
  z2 = this%z(1, iz+1)
  chi00 = this%chi(ir,   iz)
  chi01 = this%chi(ir,   iz+1)
  chi10 = this%chi(ir+1, iz)
  chi11 = this%chi(ir+1, iz+1)


  ! interpolate s within cell
  s  = ( this%s(ir,iz)*    (r2-r )*(z2-z ) + &
         this%s(ir,iz+1)*  (r2-r )*(z -z1) + &
         this%s(ir+1,iz)*  (r -r1)*(z2-z ) + &
         this%s(ir+1,iz+1)*(r -r1)*(z -z1) )/(r2-r1)/(z2-z1)


  ! interpolate chi within cell (taking into account discontinuity at 2pi -> 0)
  if ( (max(chi00,chi01,chi10,chi11) - &
        min(chi00,chi01,chi10,chi11)) < pi ) then
     chi  = ( chi00*(r2-r )*(z2-z ) + &
              chi01*(r2-r )*(z -z1) + &
              chi10*(r -r1)*(z2-z ) + &
              chi11*(r -r1)*(z -z1) )/(r2-r1)/(z2-z1)
  else
     ctmp = ( exp((0.,1.)*chi00)*(r2-r )*(z2-z ) + &
              exp((0.,1.)*chi01)*(r2-r )*(z -z1) + &
              exp((0.,1.)*chi10)*(r -r1)*(z2-z ) + &
              exp((0.,1.)*chi11)*(r -r1)*(z -z1) )/(r2-r1)/(z2-z1)
     chi  = datan2(imag(ctmp),real(ctmp))
     if (chi < 0.d0) chi = chi + pi2
  endif

  end subroutine forward_transform
  !-----------------------------------------------------------------------------
! type marsf_schimesh ==========================================================



! type marsf_bplasma ===========================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_new_bplasma(this)
  type(marsf_bplasma), intent(inout) :: this

  integer :: ns, nsp, mmaxe, mmaxp


  ns    = this%ns
  nsp   = this%nsp
  mmaxe = this%mmaxe
  mmaxp = this%mmaxp
  allocate (this%cs(ns), this%csm(ns), this%q(ns))
  allocate (this%rmi(ns,mmaxe), this%zmi(ns,mmaxe), this%rmm(ns,mmaxe), this%zmm(ns,mmaxe))
  allocate (this%b1(ns,this%mmaxp), this%b2(ns,this%mmaxp), this%b3(ns,this%mmaxp))
  allocate (this%x1(nsp+1,this%mmaxp))

  end subroutine aux_new_bplasma
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_bplasma(filename, amplitude, phase) result(this)
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: amplitude, phase
  type(marsf_bplasma)          :: this

  complex(kind=dp) :: ephi
  real(real64) :: rtmp8(8)
  integer :: i, iu, j, k, nsp, nsv, mmaxe


  open  (newunit=iu, file=filename, action='read', form='formatted')
  read  (iu, *) this%nn, mmaxe, this%m1, this%m2, nsp, nsv, this%itmp, this%rtmp
  this%mmaxe = mmaxe
  this%mmaxp = this%m2 - this%m1 + 1
  this%nsp   = nsp
  this%ns    = nsp + nsv
  call aux_new_bplasma(this)


  ! radial mesh and q-profile data
  do i=1,this%ns
     read  (iu, *) this%cs(i), this%csm(i), this%q(i)
  enddo


  ! coordinates mapping data
  do j=1,mmaxe
  do i=1,this%ns
     read  (iu, *) (rtmp8(k), k=1,8)
     this%rmi(i,j) = cmplx(rtmp8(1),rtmp8(2), kind=dp)
     this%zmi(i,j) = cmplx(rtmp8(3),rtmp8(4), kind=dp)
     this%rmm(i,j) = cmplx(rtmp8(5),rtmp8(6), kind=dp)
     this%zmm(i,j) = cmplx(rtmp8(7),rtmp8(8), kind=dp)
  enddo
  enddo


  ! b-field data
  do j=1,this%mmaxp
  do i=1,this%ns
     read  (iu, *) (rtmp8(k), k=1,6)
     this%b1(i,j) = cmplx(rtmp8(1),rtmp8(2), kind=dp)
     this%b2(i,j) = cmplx(rtmp8(3),rtmp8(4), kind=dp)
     this%b3(i,j) = cmplx(rtmp8(5),rtmp8(6), kind=dp)
  enddo
  enddo


  ! x-field data
  do j=1,this%mmaxp
  do i=1,nsp+1
     read  (iu, *) (rtmp8(k), k=1,2)
     this%x1(i,j) = cmplx(rtmp8(1),rtmp8(2), kind=dp)
  enddo
  enddo
  close (iu)


  if (present(amplitude)) then
     this%b1 = amplitude * this%b1
     this%b2 = amplitude * this%b2
     this%b3 = amplitude * this%b3
     this%x1 = amplitude * this%x1
  endif
  if (present(phase)) then
     ephi    = exp((0.d0,1.d0) * phase / 180.d0 * pi)
     this%b1 = this%b1 * ephi
     this%b2 = this%b2 * ephi
     this%b3 = this%b3 * ephi
     this%x1 = this%x1 * ephi
  endif

  end function load_bplasma
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function merge_bplasma(bplasma) result(this)
  class(marsf_bplasma), intent(in) :: bplasma(:)
  type(marsf_bplasma)              :: this

  integer :: i


  this = bplasma(1)
  do i=2,size(bplasma)
     ! verify compatibility
     if (bplasma(i)%nn    /= this%nn)    call incompatible("nn")
     if (bplasma(i)%ns    /= this%ns)    call incompatible("ns")
     if (bplasma(i)%nsp   /= this%nsp)   call incompatible("nsp")
     if (bplasma(i)%m1    /= this%m1)    call incompatible("m1")
     if (bplasma(i)%mmaxe /= this%mmaxe) call incompatible("mmaxe")
     if (bplasma(i)%mmaxp /= this%mmaxp) call incompatible("mmaxp")
     ! TODO: verify compatible cs, csm, q, rmi, zmi, rmm, zmm

     ! add field components
     this%b1 = this%b1 + bplasma(i)%b1
     this%b2 = this%b2 + bplasma(i)%b2
     this%b3 = this%b3 + bplasma(i)%b3
     this%x1 = this%x1 + bplasma(i)%x1
  enddo

  contains
  !.............................................................................
  subroutine incompatible(param)
  character(len=*), intent(in) :: param


  print 9000, param
  stop
 9000 format("ERROR: incompatible parameter ",a,"!")

  end subroutine incompatible
  !.............................................................................
  end function merge_bplasma
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine bplasma_broadcast(this)
  use moose_mpi
  class(marsf_bplasma), intent(inout) :: this


  call proc(0)%broadcast(this%nn)
  call proc(0)%broadcast(this%ns)
  call proc(0)%broadcast(this%nsp)
  call proc(0)%broadcast(this%m1)
  call proc(0)%broadcast(this%m2)
  call proc(0)%broadcast(this%mmaxe)
  call proc(0)%broadcast(this%mmaxp)

  if (rank > 0) call aux_new_bplasma(this)
  call proc(0)%broadcast(this%cs)
  call proc(0)%broadcast(this%csm)
  call proc(0)%broadcast(this%q)

  call proc(0)%broadcast(this%rmi)
  call proc(0)%broadcast(this%zmi)
  call proc(0)%broadcast(this%rmm)
  call proc(0)%broadcast(this%zmm)

  call proc(0)%broadcast(this%b1)
  call proc(0)%broadcast(this%b2)
  call proc(0)%broadcast(this%b3)
  call proc(0)%broadcast(this%x1)

  end subroutine bplasma_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine bplasma_free(this)
  class(marsf_bplasma), intent(inout) :: this


  deallocate (this%cs, this%csm)
  deallocate (this%rmi, this%zmi, this%rmm, this%zmm)
  deallocate (this%b1, this%b2, this%b3)
  deallocate (this%x1)

  end subroutine bplasma_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine save(this, filename)
  class(marsf_bplasma), intent(in) :: this
  character(len=*),     intent(in) :: filename

  integer :: i, iu, j, nsv


  open  (newunit=iu, file=filename)
  nsv = this%ns - this%nsp
  write (iu, 1001) this%nn, this%mmaxe, this%m1, this%m2, this%nsp, nsv, this%itmp, this%rtmp

  ! write radial mesh and q-profile data
  do i=1,this%ns
     write (iu, 1002) this%cs(i), this%csm(i), this%q(i)
  enddo

  ! write coordinates mapping data
  do j=1,this%mmaxe
  do i=1,this%ns
     write (iu, 1002) real(this%rmi(i,j)), aimag(this%rmi(i,j)), &
                      real(this%zmi(i,j)), aimag(this%zmi(i,j)), &
                      real(this%rmm(i,j)), aimag(this%rmm(i,j)), &
                      real(this%zmm(i,j)), aimag(this%zmm(i,j))
  enddo
  enddo

  ! write b-field data
  do j=1,this%mmaxp
  do i=1,this%ns
     write (iu, 1002) real(this%b1(i,j)), aimag(this%b1(i,j)), &
                      real(this%b2(i,j)), aimag(this%b2(i,j)), &
                      real(this%b3(i,j)), aimag(this%b3(i,j))
  enddo
  enddo

  ! write x-field data
  do j=1,this%mmaxp
  do i=1,this%nsp+1
     write (iu, 1002) real(this%x1(i,j)), aimag(this%x1(i,j))
  enddo
  enddo
  close (iu)
 1001 format(7(i4,1x),2(e17.10,1x))
 1002 format(8(e17.10,1x))

  end subroutine save
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine find_s_index(this, s, is, ierr)
  !
  ! For given s, find corresponding index is.
  ! ierr =
  !        0:    successful
  !        3:    s out of range
  !
  use moose_algorithms, only: binary_search_L
  class(marsf_bplasma), intent(in)  :: this
  real(real64),         intent(in)  :: s
  integer,              intent(out) :: is, ierr


  ierr = 0
  is   = binary_search_L(this%cs(:), s)
  if (is < 1  .or.  is >= this%ns) then
     ierr = 3
     return
  endif

  end subroutine find_s_index
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_x1(this, s, chi, is) result(x1)
  class(marsf_bplasma), intent(in) :: this
  real(real64),         intent(in) :: s, chi
  integer,              intent(in), optional :: is
  complex(kind=dp)                 :: x1

  complex(kind=dp) :: y0, y1, y2, y3
  real(real64) :: d, h1, h2, h3
  integer :: i, ierr, j


  x1 = 0.d0
  if (present(is)) then
     i = is
  else
     call this%find_s_index(s, i, ierr)
     if (ierr /= 0) return
  endif


  h1 = this%cs(i)    - s
  h2 = this%cs(i+1)  - s
  h3 = this%cs(i+2)  - s
  d  = (h1-h2)*(h2-h3)*(h3-h1)
  do j=1,this%mmaxp
     y1 = this%x1(i,j)
     y2 = this%x1(i+1,j)
     y3 = this%x1(i+2,j)
     y0 = (h3*h2*(h3-h2)*y1 + h1*h3*(h1-h3)*y2 + h2*h1*(h2-h1)*y3) / d
     x1 = x1 + y0 * exp((0.d0,1.d0)*(j-1+this%m1)*chi)
  enddo

  end function eval_x1
  !-----------------------------------------------------------------------------
! type marsf_bplasma ===========================================================



! type marsf_bfield ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_bfield(schimesh, bplasma) result(this)
  use flare_bfield, only: init_magnetic_field
  type(marsf_schimesh), target, intent(in) :: schimesh
  type(marsf_bplasma),  target, intent(in) :: bplasma
  type(marsf_bfield)                       :: this

  real(real64) :: rmin, rmax, zmin, zmax


  this%schimesh => schimesh
  this%bplasma  => bplasma

  rmin = minval(schimesh%r)
  rmax = maxval(schimesh%r)
  zmin = minval(schimesh%z)
  zmax = maxval(schimesh%z)
  call init_magnetic_field(this, rmin, rmax, zmin, zmax, abs(bplasma%nn))

  end function init_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_bfield(schimesh, bplasma, amplitude, phase) result(this)
  use moose_utils,  only: basename
  character(len=*), intent(in) :: schimesh, bplasma
  real(real64),     intent(in), optional :: amplitude, phase
  type(marsf_bfield)           :: this

  type(marsf_schimesh), pointer :: G
  type(marsf_bplasma),  pointer :: B
  real(real64) :: a


  if (report) then
     print *
     print 1000, basename(bplasma)
  endif
 1000 format(3x,"- Magnetic field from MARS-F: ", a)


  allocate (G, source=marsf_schimesh(schimesh))
  allocate (B, source=marsf_bplasma(bplasma, amplitude, phase))
  this = init_bfield(G, B)
  if (report) then
     if (present(phase)) then
        a = 1.d0;   if (present(amplitude)) a = amplitude
        print 1001, abs(this%bplasma%nn), a, phase
     elseif (present(amplitude)) then
        print 1002, abs(this%bplasma%nn), amplitude
     else
        print 1003, abs(this%bplasma%nn)
     endif
  endif
 1001 format(8x,"n = ",i0,", amplitude factor = ",f0.3,", phase = ",f0.3," deg")
 1002 format(8x,"n = ",i0,", amplitude factor = ",f0.3)
 1003 format(8x,"n = ",i0)

  end function load_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_bfield(this)
  use moose_mpi
  class(marsf_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  if (rank > 0) then
     allocate (this%schimesh)
     allocate (this%bplasma)
  endif
  call this%schimesh%broadcast()
  call this%bplasma%broadcast()

  end subroutine broadcast_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_bfield(this)
  class(marsf_bfield), intent(inout) :: this


  call this%schimesh%free()
  call this%bplasma%free()
  call this%mfunc_free()
  deallocate (this%schimesh, this%bplasma)

  end subroutine free_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(marsf_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: B(this%mdim)

  integer :: ierr


  call eval_field(this%schimesh, this%bplasma, B_FIELD, x, B, ierr)
  if (ierr /= 0) B = 0.d0

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(marsf_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: jac(this%mdim, this%ndim)


  jac = 0.d0
  print *, "ERROR: Jacobian of MARS-F magnetic field not implemented"
  stop

  end function jac
  !-----------------------------------------------------------------------------
! type marsf_bfield ============================================================


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine eval_geometry(this, s, i, chi, r, z, drds, drdc, dzds, dzdc)
  class(marsf_bplasma), intent(in   ) :: this
  real(real64),         intent(in   ) :: s, chi
  integer,              intent(in   ) :: i
  real(real64),         intent(  out) :: r, z, drds, drdc, dzds, dzdc

  real(real64)     :: h1, h2, h3, d
  complex(kind=dp) :: y0, y1, y2, y3, yp, ctmp
  integer :: j


  ! on the s-mesh, compute dR/ds,dZ/ds,dR/dchi,dZ/dchi, at (s0,chi0)
  ! using 3 points: 2 integer  points and 1 middle point
  h1   = this%cs(i)   - s
  h2   = this%csm(i)  - s
  h3   = this%cs(i+1) - s
  d    = (h1-h2)*(h2-h3)*(h3-h1)
  do j=1,this%mmaxe
     y1 = this%rmi(i,j)
     y2 = this%rmm(i,j)
     y3 = this%rmi(i+1,j)
     y0 = (h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
     yp = ((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d

     if (j == 1) then
        r    = real(y0)
        drds = real(yp)
        drdc = 0.d0
     else
        ctmp = exp((0.d0,1.d0)*(j-1)*chi)
        r    = r    + 2.d0*real(y0*ctmp)
        drds = drds + 2.d0*real(yp*ctmp)
        drdc = drdc + 2.d0*real(y0*(0.d0,1.d0)*(j-1)*ctmp)
     endif

     y1 = this%zmi(i,j)
     y2 = this%zmm(i,j)
     y3 = this%zmi(i+1,j)
     y0 = (h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
     yp = ((h2-h3)*(h2+h3)*y1+(h3-h1)*(h3+h1)*y2+(h1-h2)*(h1+h2)*y3)/d

     if (j == 1) then
        z    = real(y0)
        dzds = real(yp)
        dzdc = 0.d0
     else
        z    = z    + 2.d0*real(y0*ctmp)
        dzds = dzds + 2.d0*real(yp*ctmp)
        dzdc = dzdc + 2.d0*real(y0*(0.d0,1.d0)*(j-1)*ctmp)
     endif
  enddo

  end subroutine eval_geometry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine eval_field(schimesh, bplasma, type_field, x, f, ierr)
  class(marsf_schimesh), intent(in   ) :: schimesh
  class(marsf_bplasma),  intent(in   ) :: bplasma
  integer,               intent(in   ) :: type_field
  real(real64),          intent(in   ) :: x(3)
  real(real64),          intent(  out) :: f(3)
  integer,               intent(  out) :: ierr

  real(real64)     :: s, chi, sg22, h1, h2, h3, h4, h5, d, drds, drdc, dzds, dzdc, rtmp, ztmp, jac
  complex(kind=dp) :: y0, y1, y2, y3, ctmp, b10, b20, b30, x10, br, bz, bp
  integer :: i, j


  ! compute (s, chi) coordinates at (R, Z) = x(1:2)
  call schimesh%forward_transform(x(1), x(2), s, chi, ierr)
  if (ierr /= 0) return

  ! near magnetic axis
  if (s < bplasma%cs(2)) then
     chi = datan2(x(2)-schimesh%zaxis, x(1)-schimesh%raxis)
     if (chi < 0.d0) chi = chi + pi2
  endif


  ! compute geometry parameters at (s, chi)
  call bplasma%find_s_index(s, i, ierr)
  if (ierr /= 0) return

  call eval_geometry(bplasma, s, i, chi, rtmp, ztmp, drds, drdc, dzds, dzdc)
  jac = (drds*dzdc - drdc*dzds) * x(1)
  if (type_field == X_FIELD) sg22 = sqrt(drdc**2 + dzdc**2)


  ! compute b10,b20,b30,x10 at (s0,chi0)
  ! b30 is computed only if type_field=1
  ! x10 is computed only if type_field=3
  ! go through all toroidal harmonics
  if (s < 1.d0  .and.  i > 1) i = i-1
  h1 = bplasma%cs(i)    - s
  h2 = bplasma%cs(i+1)  - s
  h3 = bplasma%cs(i+2)  - s
  h4 = bplasma%csm(i)   - s
  h5 = bplasma%csm(i+1) - s
  d  = (h1-h2)*(h2-h3)*(h3-h1)

  f   = 0.d0
  b10 = (0.d0,0.d0)
  b20 = (0.d0,0.d0)
  b30 = (0.d0,0.d0)
  x10 = (0.d0,0.d0)
  do j=1,bplasma%mmaxp
     ctmp = exp((0.d0,1.d0)*(j-1+bplasma%m1)*chi)

     if (type_field == B_FIELD  .or.  type_field == A_FIELD) then
        y1 = bplasma%b1(i,j)
        y2 = bplasma%b1(i+1,j)
        y3 = bplasma%b1(i+2,j)
        y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
        b10 = b10 + y0*ctmp

        y1 = (h5*bplasma%b2(i,j)-h4*bplasma%b2(i+1,j))/(h5-h4)
        b20 = b20 + y1*ctmp
     endif

     if (type_field == B_FIELD) then
        y1 = (h5*bplasma%b3(i,j)-h4*bplasma%b3(i+1,j))/(h5-h4)
        b30 = b30 + y1*ctmp
     endif

     if (type_field == X_FIELD  .and.  s <= 1.d0) then
        y1 = bplasma%x1(i,j)
        y2 = bplasma%x1(i+1,j)
        y3 = bplasma%x1(i+2,j)
        y0=(h3*h2*(h3-h2)*y1+h1*h3*(h1-h3)*y2+h2*h1*(h2-h1)*y3)/d
        x10 = x10 + y0*ctmp
     endif
  enddo


  ! compute toroidal harmonics of B- or A-field at (s0,chi0)
  if (type_field == B_FIELD) then     !B-field
         br = (b10*drds+b20*drdc)/jac
         bp =-b30*x(1)/jac
         bz = (b10*dzds+b20*dzdc)/jac
  elseif (type_field == A_FIELD) then !A-field
         br = (b10*dzds+b20*dzdc)*x(1)/jac/((0.d0, 1.d0)*bplasma%nn)
         bp = (0.d0, 0.d0)
         bz = (b10*drds+b20*drdc)*x(1)/jac/(-(0.d0, 1.d0)*bplasma%nn)
  elseif (type_field == X_FIELD) then !xi_n
         br = x10*jac*sg22/x(1)
         bp = (0.d0, 0.d0)
         bz = (0.d0, 0.d0)
  endif


  ! compute field in real space
  ctmp = exp(-(0.d0, 1.d0)*bplasma%nn * x(3))
  if (type_field == B_FIELD) then
     f(1) = real(br*ctmp)
     f(2) = real(bz*ctmp)
     f(3) = real(bp*ctmp)
  elseif (type_field == A_FIELD) then
     f(1) = real(br*ctmp)
     f(2) = real(bz*ctmp)
  elseif (type_field == X_FIELD) then
     f(1) = real(br*ctmp)
  endif

  end subroutine eval_field
  !-----------------------------------------------------------------------------

end module flare_marsf
