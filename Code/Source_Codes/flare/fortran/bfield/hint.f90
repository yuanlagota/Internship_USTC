module flare_hint
  use iso_fortran_env
  use flare_bfield
  implicit none
  private


  integer, parameter :: l3d = 4

  real(real64), parameter :: &
                         c411 = -5.0d0      / 2048.0d0,   &
                         c412 =  9611.0d0   / 737280.0d0, &
                         c413 =  259.0d0    / 23040.0d0,  &
                         c414 = -6629.0d0   / 46080.0d0,  &
                         c415 = -7.0d0      / 1152.0d0,   &
                         c416 =  819.0d0    / 1024.0d0,   &
                         c417 =  1.0d0      / 1440.0d0,   &
                         c418 = -1067.0d0   / 360.0d0,    &
                         c41a =  3941.0d0   / 576.0d0,    &
                         c41c = -1603.0d0   / 180.0d0,    &
                         c41e =  901.0d0    / 180.0d0,    &
                         c421 =  49.0d0     / 2048.0d0,   &
                         c422 = -70733.0d0  / 737280.0d0, &
                         c423 = -499.0d0    / 4608.0d0,   &
                         c424 =  47363.0d0  / 46080.0d0,  &
                         c425 =  59.0d0     / 1152.0d0,   &
                         c426 = -86123.0d0  / 15360.0d0,  &
                         c427 = -1.0d0      / 288.0d0,    &
                         c431 = -245.0d0    / 2048.0d0,   &
                         c432 =  27759.0d0  / 81920.0d0,  &
                         c433 =  1299.0d0   / 2560.0d0,   &
                         c434 = -50563.0d0  / 15360.0d0,  &
                         c435 = -15.0d0     / 128.0d0,    &
                         c436 =  51725.0d0  / 3072.0d0,   &
                         c437 =  1.0d0      / 160.0d0,    &
                         c441 =  1225.0d0   / 2048.0d0,   &
                         c442 = -240077.0d0 / 147456.0d0, &
                         c443 = -1891.0d0   / 4608.0d0,   &
                         c444 =  52931.0d0  / 9216.0d0,   &
                         c445 =  83.0d0     / 1152.0d0,   &
                         c446 = -86251.0d0  / 3072.0d0,   &
                         c447 = -1.0d0      / 288.0d0,    &
                         d413 =  c413 * 2.0d0,                &
                         d423 =  c423 * 2.0d0,                &
                         d433 =  c433 * 2.0d0,                &
                         d443 =  c443 * 2.0d0,                &
                         d414 =  c414 * 3.0d0,                &
                         d424 =  c424 * 3.0d0,                &
                         d434 =  c434 * 3.0d0,                &
                         d444 =  c444 * 3.0d0,                &
                         d415 =  c415 * 4.0d0,                &
                         d425 =  c425 * 4.0d0,                &
                         d435 =  c435 * 4.0d0,                &
                         d445 =  c445 * 4.0d0,                &
                         d416 =  c416 * 5.0d0,                &
                         d426 =  c426 * 5.0d0,                &
                         d436 =  c436 * 5.0d0,                &
                         d446 =  c446 * 5.0d0,                &
                         d417 =  c417 * 6.0d0,                &
                         d427 =  c427 * 6.0d0,                &
                         d437 =  c437 * 6.0d0,                &
                         d447 =  c447 * 6.0d0,                &
                         d418 =  c418 * 7.0d0,                &
                         d41a =  c41a * 9.0d0,                &
                         d41c =  c41c * 11.0d0,               &
                         d41e =  c41e * 13.0d0


  type, extends(magnetic_field), public :: hint_bfield
     real(real64), allocatable :: f3d(:,:,:,:)
     real(real64) :: h3x, h3y, h3z
     integer :: nx3d, ny3d, nz3d

     contains
     procedure :: broadcast
     procedure :: free

     procedure, private :: spl3df, spl3dd
     procedure :: eval
     procedure :: jac
  end type hint_bfield


  public :: &
     load_hint_bfield

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init(this, nr, nz, nt, mtor, rmin, rmax, zmin, zmax, bmax)
  class(hint_bfield), intent(inout) :: this
  integer,            intent(in   ) :: nr, nz, nt, mtor
  real(real64),       intent(in   ) :: rmin, rmax, zmin, zmax, bmax

  real(real64), allocatable :: rg(:), zg(:)
  real(real64) :: bb, eps, r, z
  integer :: i, j, k, l


  this%nx3d = nr
  this%ny3d = nz
  this%nz3d = nt + 1
  call init_magnetic_field(this, rmin, rmax, zmin, zmax, mtor)
  this%h3x = (this%ub(1) - this%lb(1)) / (this%nx3d - 1)
  this%h3y = (this%ub(2) - this%lb(2)) / (this%ny3d - 1)
  this%h3z = (this%ub(3) - this%lb(3)) / (this%nz3d - 1)
  allocate (rg(nr), zg(nz))
  do i=1,nr
     rg(i) = rmin + this%h3x * (i-1)
  enddo
  do j=1,nz
     zg(j) = zmin + this%h3y * (j-1)
  enddo


  ! enforce max. field value
  this%f3d(4,:,:,:) = sqrt(this%f3d(1,:,:,:)**2 + this%f3d(2,:,:,:)**2 + this%f3d(3,:,:,:)**2)
  do k=1,nt
     do j=1,this%ny3d
        do i=1,this%nx3d
           bb = this%f3d(4,i,j,k)
           if (bb > bmax) then
              this%f3d(1,i,j,k) =  bmax * this%f3d(1,i,j,k) / bb
              this%f3d(2,i,j,k) =  bmax * this%f3d(2,i,j,k) / bb
              this%f3d(3,i,j,k) =  bmax * this%f3d(3,i,j,k) / bb
           endif
        enddo
     enddo
  enddo


  ! boundary condition (radial direction)
  do k=1,nt
     do j=1,nz
        do i=1,3
           do l=1,l3d
              r = rmin - this%h3x * (4 - i)
              call polint(rg(1:4),     this%f3d(l,1:4,j,k),     4, r, this%f3d(l,i-3,j,k),  eps)
              r = rmax + this%h3x * i
              call polint(rg(nr-3:nr), this%f3d(l,nr-3:nr,j,k), 4, r, this%f3d(l,nr+i,j,k), eps)
           enddo
        enddo
     enddo
  enddo


  ! boundary condition (vertical direction)
  do k=1,nt
     do i=-2,nr+3
        do j=1,3
           do l=1,l3d
              z = zmin - this%h3y * (4 - j)
              call polint(zg(1:4),     this%f3d(l,i,1:4,k),     4, z, this%f3d(l,i,j-3,k),  eps)
              z = zmax + this%h3y * j
              call polint(zg(nz-3:nz), this%f3d(l,i,nz-3:nz,k), 4, z, this%f3d(l,i,nz+j,k), eps)
           enddo
        enddo
     enddo
  enddo
  deallocate (rg, zg)


  ! boundary condition (toroidal direction)
  if (nt == 1) then
     this%f3d(:,:,:,nt+1) =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,nt+2) =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,nt+3) =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,nt+4) =  this%f3d(:,:,:,1)

     this%f3d(:,:,:,-2)     =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,-1)     =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,0)      =  this%f3d(:,:,:,1)
  else
     this%f3d(:,:,:,nt+1) =  this%f3d(:,:,:,1)
     this%f3d(:,:,:,nt+2) =  this%f3d(:,:,:,2)
     this%f3d(:,:,:,nt+3) =  this%f3d(:,:,:,3)
     this%f3d(:,:,:,nt+4) =  this%f3d(:,:,:,4)

     this%f3d(:,:,:,-2)     =  this%f3d(:,:,:,nt-2)
     this%f3d(:,:,:,-1)     =  this%f3d(:,:,:,nt-1)
     this%f3d(:,:,:,0)      =  this%f3d(:,:,:,nt)
  endif

  end subroutine aux_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_hint_bfield(filename, group, bmax) result(this)
  use moose_error
  use moose_utils, only: basename, endswith, str
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: group
  real(real64),     intent(in) :: bmax
  type(hint_bfield)            :: this

  real(real64), allocatable :: br(:,:,:), bz(:,:,:), bp(:,:,:)
  real(real64) :: time, tmp, rmin, rmax, zmin, zmax
  integer :: i, itmp(4), iu, j, k, kstep, nr, nz, nt, mtor


  print *
  print 1000, basename(filename)
 1000 format(3x,"- Magnetic field from HINT: ",a)


  if (endswith(filename, ".h5")) then
     call loadnc()
  else
     call loadtxt()
  endif


  if (bmax == 0.d0) then
     call aux_init(this, nr, nz, nt, mtor, rmin, rmax, zmin, zmax, huge(1.d0))
  else
     call aux_init(this, nr, nz, nt, mtor, rmin, rmax, zmin, zmax, bmax)
  endif
  print 1001
  print 1002, mtor, rmin, rmax, zmin, zmax
  print *
  print 1003
  print 1004, nr, nz, nt
 1001 format(8x,"domain:    m   rmin     rmax     zmin     zmax ")
 1002 format(15x,i5,4f9.4)
 1003 format(8x,"resolution:    nr    nz    nt")
 1004 format(20x,3i6)

  contains
  !.............................................................................
  subroutine loadtxt()

  logical :: plasma


  select case(group)
  case (-1)
     plasma = .true.
  case (0)
     plasma = .false.
  case default
     call ERROR("invalid group '"//str(group)//"'")
  end select


  ! header -> domain + resolution
  open  (newunit=iu, file=filename, form='unformatted', status='old', convert='big_endian')
  if (plasma) then
     read (iu) kstep
     read (iu) time
  endif
  read (iu) nr, nz, nt, mtor
  read (iu) rmin, zmin, rmax, zmax
  allocate (this%f3d(l3d,-2:nr+3,-2:nz+3,-2:nt+4), source=0.d0)
  allocate (br(nr,nz,nt), bp(nr,nz,nt), bz(nr,nz,nt))


  ! magnetic field data
  do
     read (iu) (((br(i,j,k), bp(i,j,k), bz(i,j,k), i=1,nr), j=1,nz), k=1,nt)
     if (.not.plasma) exit

     read (iu) (((tmp, tmp, tmp, i=1,nr), j=1,nz), k=1,nt)
     read (iu) (((tmp, i=1,nr), j=1,nz), k=1,nt)

     read (iu, end=4) kstep
     read (iu) time
     read (iu) (itmp(i), i=1,4)
     read (iu) tmp, tmp, tmp, tmp
  enddo
4 close (iu)
  this%f3d(1,1:nr,1:nz,1:nt) = br
  this%f3d(2,1:nr,1:nz,1:nt) = bz
  this%f3d(3,1:nr,1:nz,1:nt) = bp
  deallocate (br, bz, bp)

  end subroutine loadtxt
  !.............................................................................
  subroutine loadnc()
  use moose_netcdf

  type(netcdf_dataset) :: N, G
  integer :: ngroups, group_ids(1024)


  N = netcdf_open(filename)

  ! mesh layout
  call N%get_var("nr", nr)
  call N%get_var("nz", nz)
  call N%get_var("ntor", nt)
  call N%get_var("mtor", mtor)

  call N%get_var("rminb", rmin)
  call N%get_var("rmaxb", rmax)
  call N%get_var("zminb", zmin)
  call N%get_var("zmaxb", zmax)


  ! select data group
  call N%inq_grps(ngroups, group_ids)
  if (group >= ngroups  .or.  group < -ngroups) then
     call ERROR("selected group out of bounds")
  elseif (group >= 0) then
     G%ncid = group_ids(group+1)
  else
     G%ncid = group_ids(ngroups+1+group)
  endif


  ! magnetic field data
  allocate (br(nr, nz, nt), bz(nr, nz, nt), bp(nr, nz, nt))
  call G%get_var("B_R", br)
  call G%get_var("B_Z", bz)
  call G%get_var("B_phi", bp)
  allocate (this%f3d(l3d,-2:nr+3,-2:nz+3,-2:nt+4), source=0.d0)
  this%f3d(1,1:nr,1:nz,1:nt) = br
  this%f3d(2,1:nr,1:nz,1:nt) = bz
  this%f3d(3,1:nr,1:nz,1:nt) = bp
  deallocate (br, bz, bp)

  call N%close()

  end subroutine loadnc
  !.............................................................................
  end function load_hint_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  use moose_error
  class(hint_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  call proc(0)%broadcast_allocatable(this%f3d)
  call proc(0)%broadcast(this%h3x)
  call proc(0)%broadcast(this%h3y)
  call proc(0)%broadcast(this%h3z)
  call proc(0)%broadcast(this%nx3d)
  call proc(0)%broadcast(this%ny3d)
  call proc(0)%broadcast(this%nz3d)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(hint_bfield), intent(inout) :: this


  call this%mfunc_free()
  deallocate (this%f3d)

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine spl3df(this, xd, w)
  class(hint_bfield), intent(in   ) :: this
  real(real64),       intent(in   ) :: xd(3)
  real(real64),       intent(  out) :: w(l3d)

  integer, parameter :: m2 = 8

  real(real64) :: x, y, z, ux, uy, uz, usx, usy, usz, x400, y400, z400, x41m, y41m, z41m, &
                  x41p, y41p, z41p, x42m, y42m, z42m, x42p, y42p, z42p, x43m, y43m, z43m, &
                  x43p, y43p, z43p, x44m, y44m, z44m, x44p, y44p, z44p, cyz, cy(8), cz(8)
  integer :: l, mz, lz, my, ly, ix, iy, iz


  w = 0.0d0


  x = xd(1)
  if (x < this%lb(1)) return
  if (x > this%ub(1)) return

  y = xd(2)   
  if (y < this%lb(2)) return
  if (y > this%ub(2)) return

  z = xd(3)
  if (z < this%lb(3)) z = this%lb(3)
  if (z > this%ub(3)) z = this%ub(3)

  ux = (x - this%lb(1)) / this%h3x
  ix = ux
  if (ix >= this%nx3d - 1) then
    ix = this%nx3d - 2
    ux = this%nx3d - 1
  endif
 
  uy = (y - this%lb(2)) / this%h3y
  iy = uy
  if (iy >= this%ny3d - 1) then
    iy = this%ny3d - 2
    uy = this%ny3d - 1
  endif

  uz = (z - this%lb(3)) / this%h3z
  iz = uz
  if (iz >= this%nz3d - 1) then
    iz = this%nz3d - 2
    uz = this%nz3d - 1
  endif

  ux  = ux - ix - 0.5d0
  uy  = uy - iy - 0.5d0
  uz  = uz - iz - 0.5d0
  usx = ux * ux
  usy = uy * uy
  usz = uz * uz


  x400 = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
  y400 = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
  z400 = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
  
  x41m = (((x400       + c416) * usx + c414) * usx + c412) * ux
  x41p =  ((c417 * usx + c415) * usx + c413) * usx + c411
  y41m = (((y400       + c416) * usy + c414) * usy + c412) * uy
  y41p =  ((c417 * usy + c415) * usy + c413) * usy + c411
  z41m = (((z400       + c416) * usz + c414) * usz + c412) * uz
  z41p =  ((c417 * usz + c415) * usz + c413) * usz + c411

  x42m = (((-x400 * 7.0d0 + c426) * usx + c424) * usx + c422) * ux
  x42p =  ((c427  * usx   + c425) * usx + c423) * usx + c421
  y42m = (((-y400 * 7.0d0 + c426) * usy + c424) * usy + c422) * uy
  y42p =  ((c427  * usy   + c425) * usy + c423) * usy + c421
  z42m = (((-z400 * 7.0d0 + c426) * usz + c424) * usz + c422) * uz
  z42p =  ((c427  * usz   + c425) * usz + c423) * usz + c421

  x43m = (((x400 * 21.0d0 + c436) * usx + c434) * usx + c432) * ux
  x43p =  ((c437 * usx    + c435) * usx + c433) * usx + c431
  y43m = (((y400 * 21.0d0 + c436) * usy + c434) * usy + c432) * uy
  y43p =  ((c437 * usy    + c435) * usy + c433) * usy + c431
  z43m = (((z400 * 21.0d0 + c436) * usz + c434) * usz + c432) * uz
  z43p =  ((c437 * usz    + c435) * usz + c433) * usz + c431

  x44m = (((-x400 * 35.0d0 + c446) * usx + c444) * usx + c442) * ux
  x44p =  ((c447  * usx    + c445) * usx + c443) * usx + c441
  y44m = (((-y400 * 35.0d0 + c446) * usy + c444) * usy + c442) * uy
  y44p =  ((c447  * usy    + c445) * usy + c443) * usy + c441
  z44m = (((-z400 * 35.0d0 + c446) * usz + c444) * usz + c442) * uz
  z44p =  ((c447  * usz    + c445) * usz + c443) * usz + c441

  cy(1) =  y41p + y41m
  cy(2) =  y42p + y42m
  cy(3) =  y43p + y43m
  cy(4) =  y44p + y44m
  cy(5) =  y44p - y44m
  cy(6) =  y43p - y43m
  cy(7) =  y42p - y42m
  cy(8) =  y41p - y41m
  cz(1) =  z41p + z41m
  cz(2) =  z42p + z42m
  cz(3) =  z43p + z43m
  cz(4) =  z44p + z44m
  cz(5) =  z44p - z44m
  cz(6) =  z43p - z43m
  cz(7) =  z42p - z42m
  cz(8) =  z41p - z41m

  loop010 : do l=1,l3d
     loop020 : do mz=1,m2
        lz = iz + mz - 3
        loop030 : do my=1,m2
           ly    = iy + my - 3
           cyz   = cy(my) * cz(mz)
           w(l) = ((x41p + x41m) * this%f3d(l,ix-2,ly,lz) &
                +  (x42p + x42m) * this%f3d(l,ix-1,ly,lz) &
                +  (x43p + x43m) * this%f3d(l,ix  ,ly,lz) &
                +  (x44p + x44m) * this%f3d(l,ix+1,ly,lz) &
                +  (x44p - x44m) * this%f3d(l,ix+2,ly,lz) & 
                +  (x43p - x43m) * this%f3d(l,ix+3,ly,lz) &
                +  (x42p - x42m) * this%f3d(l,ix+4,ly,lz) &
                +  (x41p - x41m) * this%f3d(l,ix+5,ly,lz)) * cyz + w(l)
        enddo loop030
     enddo loop020
  enddo loop010

  end subroutine spl3df
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine spl3dd(this, xd, w, wx, wy, wz)
  class(hint_bfield), intent(in   ) :: this
  real(real64),       intent(in   ) :: xd(3)
  real(real64),       intent(  out) :: w(l3d), wx(l3d), wy(l3d), wz(l3d)

  integer, parameter :: m2 = 8

  real(real64) :: x, y, z, ux, uy, uz, usx, usy, usz, x400, y400, z400, &
                  x41m, y41m, z41m, x41p, y41p, z41p, x42m, y42m, z42m, x42p, y42p, z42p, &
                  x43m, y43m, z43m, x43p, y43p, z43p, x44m, y44m, z44m, x44p, y44p, z44p, &
                  dx400, dy400, dz400, dx41m, dy41m, dz41m, dx41p, dy41p, dz41p, &
                  dx42m, dy42m, dz42m, dx42p, dy42p, dz42p, dx43m, dy43m, dz43m, &
                  dx43p, dy43p, dz43p, dx44m, dy44m, dz44m, dx44p, dy44p, dz44p, &
                  cx1, cx2, cx3, cx4, cx5, cx6, cx7, cx8, cy1, cy2, cy3, cy4, cy5, cy6, cy7, cy8, &
                  dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, &
                  w0l, wxl, wyl, wzl, w01, w02, w03, w04, w05, w06, w07, w08, w00, cz(8), dz(8)
  integer :: l, ix, iy, iz, mz, lz, ix1, ix2, ix4, ix5, ix6, ix7, ix8, &
             iy1, iy2, iy4, iy5, iy6, iy7, iy8

  w  = 0.d0
  wx = 0.d0
  wy = 0.d0
  wz = 0.d0

  x = xd(1)
  if (x < this%lb(1)) return
  if (x > this%ub(1)) return

  y = xd(2)
  if (y < this%lb(1)) return
  if (y > this%ub(2)) return

  z = xd(3)
  if (z < this%lb(3)) z = this%lb(3)
  if (z > this%ub(3)) z = this%ub(3)

  ux = (x - this%lb(1)) / this%h3x
  ix = ux
  if (ix >= this%nx3d - 1) then
    ix = this%nx3d - 2
    ux = this%nx3d - 1
  endif
 
  uy = (y - this%lb(2)) / this%h3y
  iy = uy
  if (iy >= this%ny3d - 1) then
    iy = this%ny3d - 2
    uy = this%ny3d - 1
  endif

  uz = (z - this%lb(3)) / this%h3z
  iz = uz
  if (iz >= this%nz3d - 1) then
    iz = this%nz3d - 2
    uz = this%nz3d - 1
  endif

  ux  = ux - ix - 0.5d0
  uy  = uy - iy - 0.5d0
  uz  = uz - iz - 0.5d0

  usx = ux * ux
  usy = uy * uy
  usz = uz * uz

  ix1 = ix - 2
  ix2 = ix - 1
  ix4 = ix + 1
  ix5 = ix + 2
  ix6 = ix + 3
  ix7 = ix + 4
  ix8 = ix + 5

  iy1 = iy - 2
  iy2 = iy - 1
  iy4 = iy + 1
  iy5 = iy + 2
  iy6 = iy + 3
  iy7 = iy + 4
  iy8 = iy + 5


  x400  = (((c41e * usx + c41c) * usx + c41a) * usx + c418) * usx
  y400  = (((c41e * usy + c41c) * usy + c41a) * usy + c418) * usy
  z400  = (((c41e * usz + c41c) * usz + c41a) * usz + c418) * usz
  dx400 = (((d41e * usx + d41c) * usx + d41a) * usx + d418) * usx
  dy400 = (((d41e * usy + d41c) * usy + d41a) * usy + d418) * usy
  dz400 = (((d41e * usz + d41c) * usz + d41a) * usz + d418) * usz

  x41m  = (((x400       + c416) * usx + c414) * usx + c412) * ux
  y41m  = (((y400       + c416) * usy + c414) * usy + c412) * uy
  z41m  = (((z400       + c416) * usz + c414) * usz + c412) * uz
  x41p  =  ((c417 * usx + c415) * usx + c413) * usx + c411
  y41p  =  ((c417 * usy + c415) * usy + c413) * usy + c411
  z41p  =  ((c417 * usz + c415) * usz + c413) * usz + c411
  dx41m =  ((dx400      + d416) * usx + d414) * usx + c412
  dy41m =  ((dy400      + d416) * usy + d414) * usy + c412
  dz41m =  ((dz400      + d416) * usz + d414) * usz + c412
  dx41p =  ((d417 * usx + d415) * usx + d413) * ux
  dy41p =  ((d417 * usy + d415) * usy + d413) * uy
  dz41p =  ((d417 * usz + d415) * usz + d413) * uz

  x42m  = (((-x400  * 7.0d0 + c426) * usx + c424) * usx + c422) * ux
  y42m  = (((-y400  * 7.0d0 + c426) * usy + c424) * usy + c422) * uy
  z42m  = (((-z400  * 7.0d0 + c426) * usz + c424) * usz + c422) * uz
  x42p  =  ((c427   * usx   + c425) * usx + c423) * usx + c421
  y42p  =  ((c427   * usy   + c425) * usy + c423) * usy + c421
  z42p  =  ((c427   * usz   + c425) * usz + c423) * usz + c421
  dx42m =  ((-dx400 * 7.0d0 + d426) * usx + d424) * usx + c422
  dy42m =  ((-dy400 * 7.0d0 + d426) * usy + d424) * usy + c422
  dz42m =  ((-dz400 * 7.0d0 + d426) * usz + d424) * usz + c422
  dx42p =  ((d427   * usx   + d425) * usx + d423) * ux
  dy42p =  ((d427   * usy   + d425) * usy + d423) * uy
  dz42p =  ((d427   * usz   + d425) * usz + d423) * uz

  x43m  = (((x400  * 21.0d0 + c436) * usx + c434) * usx + c432) * ux
  y43m  = (((y400  * 21.0d0 + c436) * usy + c434) * usy + c432) * uy
  z43m  = (((z400  * 21.0d0 + c436) * usz + c434) * usz + c432) * uz
  x43p  =  ((c437  * usx    + c435) * usx + c433) * usx + c431
  y43p  =  ((c437  * usy    + c435) * usy + c433) * usy + c431
  z43p  =  ((c437  * usz    + c435) * usz + c433) * usz + c431
  dx43m =  ((dx400 * 21.0d0 + d436) * usx + d434) * usx + c432
  dy43m =  ((dy400 * 21.0d0 + d436) * usy + d434) * usy + c432
  dz43m =  ((dz400 * 21.0d0 + d436) * usz + d434) * usz + c432
  dx43p =  ((d437  * usx    + d435) * usx + d433) * ux
  dy43p =  ((d437  * usy    + d435) * usy + d433) * uy
  dz43p =  ((d437  * usz    + d435) * usz + d433) * uz

  x44m  = (((-x400  * 35.0d0 + c446) * usx + c444) * usx + c442) * ux
  y44m  = (((-y400  * 35.0d0 + c446) * usy + c444) * usy + c442) * uy
  z44m  = (((-z400  * 35.0d0 + c446) * usz + c444) * usz + c442) * uz
  x44p  =  ((c447   * usx    + c445) * usx + c443) * usx + c441
  y44p  =  ((c447   * usy    + c445) * usy + c443) * usy + c441
  z44p  =  ((c447   * usz    + c445) * usz + c443) * usz + c441
  dx44m =  ((-dx400 * 35.0d0 + d446) * usx + d444) * usx + c442
  dy44m =  ((-dy400 * 35.0d0 + d446) * usy + d444) * usy + c442
  dz44m =  ((-dz400 * 35.0d0 + d446) * usz + d444) * usz + c442
  dx44p =  ((d447   * usx    + d445) * usx + d443) * ux
  dy44p =  ((d447   * usy    + d445) * usy + d443) * uy
  dz44p =  ((d447   * usz    + d445) * usz + d443) * uz

  cx1   = x41p +  x41m
  cx2   = x42p +  x42m
  cx3   = x43p +  x43m
  cx4   = x44p +  x44m
  cx5   = x44p -  x44m
  cx6   = x43p -  x43m
  cx7   = x42p -  x42m
  cx8   = x41p -  x41m
  cy1   = y41p +  y41m
  cy2   = y42p +  y42m
  cy3   = y43p +  y43m
  cy4   = y44p +  y44m
  cy5   = y44p -  y44m
  cy6   = y43p -  y43m
  cy7   = y42p -  y42m
  cy8   = y41p -  y41m
  cz(1) = z41p +  z41m
  cz(2) = z42p +  z42m
  cz(3) = z43p +  z43m
  cz(4) = z44p +  z44m
  cz(5) = z44p -  z44m
  cz(6) = z43p -  z43m
  cz(7) = z42p -  z42m
  cz(8) = z41p -  z41m
  dx1   = dx41p + dx41m
  dx2   = dx42p + dx42m
  dx3   = dx43p + dx43m
  dx4   = dx44p + dx44m
  dx5   = dx44p - dx44m
  dx6   = dx43p - dx43m
  dx7   = dx42p - dx42m
  dx8   = dx41p - dx41m
  dy1   = dy41p + dy41m
  dy2   = dy42p + dy42m
  dy3   = dy43p + dy43m
  dy4   = dy44p + dy44m
  dy5   = dy44p - dy44m
  dy6   = dy43p - dy43m
  dy7   = dy42p - dy42m
  dy8   = dy41p - dy41m
  dz(1) = dz41p + dz41m
  dz(2) = dz42p + dz42m
  dz(3) = dz43p + dz43m
  dz(4) = dz44p + dz44m
  dz(5) = dz44p - dz44m
  dz(6) = dz43p - dz43m
  dz(7) = dz42p - dz42m
  dz(8) = dz41p - dz41m


  loop010 : do l=1,l3d
     w0l = 0.d0
     wxl = 0.d0
     wyl = 0.d0
     wzl = 0.d0

     loop020 : do mz=1,m2
        lz = iz + mz - 3
        w01 = cx1 * this%f3d(l,ix1,iy1,lz) + cx2 * this%f3d(l,ix2,iy1,lz) &
            + cx3 * this%f3d(l,ix, iy1,lz) + cx4 * this%f3d(l,ix4,iy1,lz) &
            + cx5 * this%f3d(l,ix5,iy1,lz) + cx6 * this%f3d(l,ix6,iy1,lz) &
            + cx7 * this%f3d(l,ix7,iy1,lz) + cx8 * this%f3d(l,ix8,iy1,lz)
        w02 = cx1 * this%f3d(l,ix1,iy2,lz) + cx2 * this%f3d(l,ix2,iy2,lz) &
            + cx3 * this%f3d(l,ix, iy2,lz) + cx4 * this%f3d(l,ix4,iy2,lz) &
            + cx5 * this%f3d(l,ix5,iy2,lz) + cx6 * this%f3d(l,ix6,iy2,lz) &
            + cx7 * this%f3d(l,ix7,iy2,lz) + cx8 * this%f3d(l,ix8,iy2,lz)
        w03 = cx1 * this%f3d(l,ix1,iy, lz) + cx2 * this%f3d(l,ix2,iy, lz) &
            + cx3 * this%f3d(l,ix, iy, lz) + cx4 * this%f3d(l,ix4,iy, lz) &
            + cx5 * this%f3d(l,ix5,iy, lz) + cx6 * this%f3d(l,ix6,iy, lz) &
            + cx7 * this%f3d(l,ix7,iy, lz) + cx8 * this%f3d(l,ix8,iy, lz)
        w04 = cx1 * this%f3d(l,ix1,iy4,lz) + cx2 * this%f3d(l,ix2,iy4,lz) &
            + cx3 * this%f3d(l,ix, iy4,lz) + cx4 * this%f3d(l,ix4,iy4,lz) &
            + cx5 * this%f3d(l,ix5,iy4,lz) + cx6 * this%f3d(l,ix6,iy4,lz) &
            + cx7 * this%f3d(l,ix7,iy4,lz) + cx8 * this%f3d(l,ix8,iy4,lz)
        w05 = cx1 * this%f3d(l,ix1,iy5,lz) + cx2 * this%f3d(l,ix2,iy5,lz) &
            + cx3 * this%f3d(l,ix, iy5,lz) + cx4 * this%f3d(l,ix4,iy5,lz) &
            + cx5 * this%f3d(l,ix5,iy5,lz) + cx6 * this%f3d(l,ix6,iy5,lz) &
            + cx7 * this%f3d(l,ix7,iy5,lz) + cx8 * this%f3d(l,ix8,iy5,lz)
        w06 = cx1 * this%f3d(l,ix1,iy6,lz) + cx2 * this%f3d(l,ix2,iy6,lz) &
            + cx3 * this%f3d(l,ix, iy6,lz) + cx4 * this%f3d(l,ix4,iy6,lz) &
            + cx5 * this%f3d(l,ix5,iy6,lz) + cx6 * this%f3d(l,ix6,iy6,lz) &
            + cx7 * this%f3d(l,ix7,iy6,lz) + cx8 * this%f3d(l,ix8,iy6,lz)
        w07 = cx1 * this%f3d(l,ix1,iy7,lz) + cx2 * this%f3d(l,ix2,iy7,lz) &
            + cx3 * this%f3d(l,ix, iy7,lz) + cx4 * this%f3d(l,ix4,iy7,lz) &
            + cx5 * this%f3d(l,ix5,iy7,lz) + cx6 * this%f3d(l,ix6,iy7,lz) &
            + cx7 * this%f3d(l,ix7,iy7,lz) + cx8 * this%f3d(l,ix8,iy7,lz)
        w08 = cx1 * this%f3d(l,ix1,iy8,lz) + cx2 * this%f3d(l,ix2,iy8,lz) &
            + cx3 * this%f3d(l,ix, iy8,lz) + cx4 * this%f3d(l,ix4,iy8,lz) &
            + cx5 * this%f3d(l,ix5,iy8,lz) + cx6 * this%f3d(l,ix6,iy8,lz) &
            + cx7 * this%f3d(l,ix7,iy8,lz) + cx8 * this%f3d(l,ix8,iy8,lz)

        wxl = ((dx1 * this%f3d(l,ix1,iy1,lz) + dx2 * this%f3d(l,ix2,iy1,lz)         &
          & +   dx3 * this%f3d(l,ix, iy1,lz) + dx4 * this%f3d(l,ix4,iy1,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy1,lz) + dx6 * this%f3d(l,ix6,iy1,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy1,lz) + dx8 * this%f3d(l,ix8,iy1,lz)) * cy1  &
          & +  (dx1 * this%f3d(l,ix1,iy2,lz) + dx2 * this%f3d(l,ix2,iy2,lz)         &
          & +   dx3 * this%f3d(l,ix, iy2,lz) + dx4 * this%f3d(l,ix4,iy2,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy2,lz) + dx6 * this%f3d(l,ix6,iy2,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy2,lz) + dx8 * this%f3d(l,ix8,iy2,lz)) * cy2  &
          & +  (dx1 * this%f3d(l,ix1,iy, lz) + dx2 * this%f3d(l,ix2,iy, lz)         &
          & +   dx3 * this%f3d(l,ix, iy, lz) + dx4 * this%f3d(l,ix4,iy, lz)         &
          & +   dx5 * this%f3d(l,ix5,iy, lz) + dx6 * this%f3d(l,ix6,iy, lz)         &
          & +   dx7 * this%f3d(l,ix7,iy, lz) + dx8 * this%f3d(l,ix8,iy, lz)) * cy3  &
          & +  (dx1 * this%f3d(l,ix1,iy4,lz) + dx2 * this%f3d(l,ix2,iy4,lz)         &
          & +   dx3 * this%f3d(l,ix, iy4,lz) + dx4 * this%f3d(l,ix4,iy4,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy4,lz) + dx6 * this%f3d(l,ix6,iy4,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy4,lz) + dx8 * this%f3d(l,ix8,iy4,lz)) * cy4  &
          & +  (dx1 * this%f3d(l,ix1,iy5,lz) + dx2 * this%f3d(l,ix2,iy5,lz)         &
          & +   dx3 * this%f3d(l,ix, iy5,lz) + dx4 * this%f3d(l,ix4,iy5,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy5,lz) + dx6 * this%f3d(l,ix6,iy5,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy5,lz) + dx8 * this%f3d(l,ix8,iy5,lz)) * cy5  &
          & +  (dx1 * this%f3d(l,ix1,iy6,lz) + dx2 * this%f3d(l,ix2,iy6,lz)         &
          & +   dx3 * this%f3d(l,ix, iy6,lz) + dx4 * this%f3d(l,ix4,iy6,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy6,lz) + dx6 * this%f3d(l,ix6,iy6,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy6,lz) + dx8 * this%f3d(l,ix8,iy6,lz)) * cy6  &
          & +  (dx1 * this%f3d(l,ix1,iy7,lz) + dx2 * this%f3d(l,ix2,iy7,lz)         &
          & +   dx3 * this%f3d(l,ix, iy7,lz) + dx4 * this%f3d(l,ix4,iy7,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy7,lz) + dx6 * this%f3d(l,ix6,iy7,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy7,lz) + dx8 * this%f3d(l,ix8,iy7,lz)) * cy7  &
          & +  (dx1 * this%f3d(l,ix1,iy8,lz) + dx2 * this%f3d(l,ix2,iy8,lz)         &
          & +   dx3 * this%f3d(l,ix, iy8,lz) + dx4 * this%f3d(l,ix4,iy8,lz)         &
          & +   dx5 * this%f3d(l,ix5,iy8,lz) + dx6 * this%f3d(l,ix6,iy8,lz)         &
          & +   dx7 * this%f3d(l,ix7,iy8,lz) + dx8 * this%f3d(l,ix8,iy8,lz)) * cy8) &
          & *   cz(mz) + wxl

        w00 =  w01 * cy1 + w02 * cy2 + w03 * cy3 + w04 * cy4 &
            +  w05 * cy5 + w06 * cy6 + w07 * cy7 + w08 * cy8
        w0l =  w00 * cz(mz) + w0l
        wzl =  w00 * dz(mz) + wzl
        wyl = (w01 * dy1 + w02 * dy2 + w03 * dy3 + w04 * dy4           &
            +  w05 * dy5 + w06 * dy6 + w07 * dy7 + w08 * dy8) * cz(mz) &
            +  wyl
     enddo loop020
         
     w(l)  =  w0l
     wx(l) =  wxl / this%h3x
     wy(l) =  wyl / this%h3y
     wz(l) =  wzl / this%h3z
  enddo loop010

  end subroutine spl3dd
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(hint_bfield), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: B(this%ndim)

  real(real64) :: xmod(3), w(l3d)
  integer :: iphi


  xmod(1:2) = x(1:2)
  xmod(3) = modulo(x(3), this%ub(3))
  call this%spl3df(xmod, w)
  B = w(1:3)

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(hint_bfield), intent(in) :: this
  real(real64),       intent(in) :: x(this%ndim)
  real(real64)                   :: jac(this%mdim, this%ndim)

  real(real64) :: xmod(3), w(l3d), wx(l3d), wy(l3d), wz(l3d)
  integer :: iphi


  xmod(1:2) = x(1:2)
  xmod(3) = modulo(x(3), this%ub(3))
  call this%spl3dd(xmod, w, wx, wy, wz)
  jac(:,1) = wx(1:3)
  jac(:,2) = wy(1:3)
  jac(:,3) = wz(1:3)

  end function jac
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine polint(xa, ya, n, x, y, dy)
  use moose_error
  integer,      intent(in   ) :: n
  real(real64), intent(in   ) :: x, xa(n), ya(n)
  real(real64), intent(  out) :: y, dy

  real(real64) :: den, dif, dift, ho, hp, w, c(n), d(n)
  integer :: i, ns, m
  

  ns  = 1
  dif = abs(x - xa(1))

  loop100: do i=1,n
     dift = abs(x - xa(i))
     if (dift < dif) then
        ns  = i
        dif = dift
     endif
     c(i) = ya(i)
     d(i) = ya(i)
  enddo loop100

  y  = ya(ns)
  ns = ns - 1

  loop200: do m=1,n-1
      loop210: do i=1,n-m
         ho  = xa(i)   - x
         hp  = xa(i+m) - x
         w   = c(i+1)  - d(i)
         den = ho      - hp
         if (den == 0.d0) call ERROR("failure in polint")
         den  = w  / den
         d(i) = hp * den
         c(i) = ho * den
      enddo loop210
      if (2 * ns < n - m) then
         dy = c(ns+1)
      else
         dy = d(ns)
         ns = ns - 1
      endif
      y = y + dy
  enddo loop200

  end subroutine polint
  !-----------------------------------------------------------------------------

end module flare_hint
