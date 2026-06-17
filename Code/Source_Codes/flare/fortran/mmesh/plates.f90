module flare_mmesh_plates
  use iso_fortran_env
  use flare_mmesh, only: mmesh_workspace => workspace
  implicit none


  character(len=*), parameter :: PLATES_FILENAME = "plates.dat"


  type, private :: zone
     integer, allocatable :: plate(:,:,:)
  end type zone
  type(zone), allocatable :: workspace(:)


  contains
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine init_workspace()

  integer :: iz, nr, np, nt, nz


  nz = size(mmesh_workspace)
  allocate (workspace(0:nz-1))
  do iz=0,nz-1
     nr = mmesh_workspace(iz)%n(1)-1
     np = mmesh_workspace(iz)%n(2)-1
     nt = mmesh_workspace(iz)%n(3)-1
     allocate (workspace(iz)%plate(0:nr-1, 0:np-1, 0:nt-1), source=0)
  enddo

  end subroutine init_workspace
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine load_plates(plate_format)
  use moose_error
  integer, intent(in) :: plate_format

  integer, allocatable :: jtcell(:), jtbundle(:)
  integer :: iostat, iu, ir, ip, it, iz, k, k1, k2, kbundle, nk, nt


  select case (plate_format)
  case (1)
     print *, "Loading plate representation (bundled format)"

  case (3)
     print *, "Loading plate representation"

  case default 
     call ERROR("invalid plate format")
  end select


  nt = 0
  do iz=0,size(workspace)-1
     nt = max(nt, mmesh_workspace(iz)%n(3)-1)
  enddo
  allocate (jtcell(nt), jtbundle(nt))


  open  (newunit=iu, file=PLATES_FILENAME)
  select case (plate_format)
  ! 1. bundled format
  case (1)
     do
        read  (iu, *, iostat=iostat) iz, ir, ip, nk, jtbundle(1:nk);   if (iostat /= 0) exit
        do kbundle=1,nk/2
           k1 = (kbundle-1)*2 + 1
           k2 =  kbundle   *2

           do it=jtbundle(k1),jtbundle(k2)
              workspace(iz)%plate(ir,ip,it) = 1
           enddo
        enddo
     enddo

  ! 3. explicit format
  case (3)
     do
        read  (iu, *, iostat=iostat) iz, ir, ip, nk, jtcell(1:nk);   if (iostat /= 0) exit
        do k=1,nk
           it = jtcell(k)
           workspace(iz)%plate(ir,ip,it) = 1
        enddo
     enddo

  end select
  close (iu)
  deallocate (jtcell, jtbundle)

  end subroutine load_plates
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine broadcast_plates()
  use moose_mpi

  integer :: iz


  do iz=0,size(workspace)-1
     call proc(0)%broadcast(workspace(iz)%plate)
  enddo

  end subroutine broadcast_plates
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine save_plates(plate_format)
  use moose_error
  integer, intent(in) :: plate_format

  integer, allocatable :: jtcell(:), jtbundle(:)
  integer :: iu, ir, ip, it, iz, nr, np, nt, nplate, nbundle


  select case (plate_format)
  case (1)
     print *, "Saving plate representation (bundled format)"

  case (3)
     print *, "Saving plate representation"

  case default
     call ERROR("invalid plate format")
  end select


  open  (newunit=iu, file=PLATES_FILENAME)
  do iz=0,size(workspace)-1
     nr = mmesh_workspace(iz)%n(1)-1
     np = mmesh_workspace(iz)%n(2)-1
     nt = mmesh_workspace(iz)%n(3)-1
     allocate (jtcell(nt), jtbundle(nt))

     do ir=0,nr-1
     do ip=0,np-1
        ! get number of plate cells in flux tube
        nplate = 0
        jtcell = 0
        do it=0,nt-1
           if (workspace(iz)%plate(ir,ip,it) > 0) then
              nplate = nplate + 1
              jtcell(nplate) = it
           endif
        enddo

        ! no plate cell in flux tube?
        if (nplate == 0) cycle

        select case (plate_format)
        ! 1. bundle plate cells
        case (1)
           jtbundle = 0
           nbundle = 0
           it = 0
           scan_flux_tube: do
              ! find 1st plate cell in bundle
              if (workspace(iz)%plate(ir,ip,it) > 0) then
                 nbundle = nbundle + 1
                 jtbundle(nbundle) = it

                 ! find last plate cell in bundle
                 scan_bundle: do
                    ! reached end of flux tube
                    if (it == nt-1) exit

                    if (workspace(iz)%plate(ir,ip,it+1) == 0) exit
                    it = it + 1
                 enddo scan_bundle

                 nbundle = nbundle + 1
                 jtbundle(nbundle) = it
              endif

              ! move to next cell
              it = it + 1
              if (it >= nt) exit
           enddo scan_flux_tube
           write (iu, 1001) iz, ir, ip, nbundle, jtbundle(1:nbundle)

        ! 3. explicit definition
        case (3)
           write (iu, *) iz, ir, ip, nplate, jtcell(1:nplate)

        end select
     enddo
     enddo

     deallocate (jtcell, jtbundle)
  enddo
  close (iu)
 1001 format(1x,i0,1x,i4,1x,i4,1x,i4,*(1x,i5,1x,i5))

  end subroutine save_plates
  !-------------------------------------------------------------------


  !-------------------------------------------------------------------
  subroutine free_workspace()

  integer :: iz


  do iz=0,size(workspace)-1
     deallocate (workspace(iz)%plate)
  enddo
  deallocate (workspace)

  end subroutine free_workspace
  !-------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine exclude_cells_from_flux_tubes(iz, kmin)
  !
  ! exclude cells from flux tubes in zone *iz* if fewer than *kmin*
  ! cells are between two target plates.
  !
  integer, intent(in) :: iz, kmin

  integer, dimension(:), allocatable :: kbeg, kend
  integer :: ir, ip, it, iplate, ns, nr, np, nt, k, k1, k2


  nr = mmesh_workspace(iz)%n(1)-1
  np = mmesh_workspace(iz)%n(2)-1
  nt = mmesh_workspace(iz)%n(3)-1
  allocate (kbeg(nt), kend(nt))

  do ir=0,nr-1
  do ip=0,np-1
     ! 1. set up kbeg and kend from workspace(iz)%plate
     ns = 0
     do it=0,nt-2
        iplate = workspace(iz)%plate(ir,ip,it)

        ! first cell is already plate cell
        if (it == 0  .and.  iplate == 1) then
           ns       = ns + 1
           kbeg(ns) = it
        endif

        ! transition between plasma and plate cells
        if (iplate /= workspace(iz)%plate(ir,ip,it+1)) then
           ! plasma -> target
           if (iplate == 0) then
              ns       = ns + 1
              kbeg(ns) = it + 1

           ! target -> plasma
           else
              kend(ns) = it
           endif
        endif
     enddo
     ! last cell is plate cell
     if (workspace(iz)%plate(ir,ip,nt-1) == 1) then
        kend(ns) = it
     endif


     ! 2. adjust definition for kmin
     if (ns > 0) then
        ! update kbeg, kend
        k2 = ns
        k1 = 1
        do
           if (k1 >= k2) exit

           if (kbeg(k1+1) - kend(k1) < kmin) then
              kend(k1) = kend(k1+1)
              k2       = k2 - 1
              do k=k1+2,k2
                 kbeg(k) = kbeg(k+1)
                 kend(k) = kend(k+1)
              enddo
           else
              k1 = k1 + 1
           endif
        enddo
        ns = k2

        ! update workspace(iz)%plate
        do k=1,ns
           do it=kbeg(k),kend(k)
              workspace(iz)%plate(ir,ip,it) = 1
           enddo
        enddo
     endif
  enddo
  enddo
  deallocate (kbeg, kend)

  end subroutine exclude_cells_from_flux_tubes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine scan_cells(iz, match, samples)
  use moose_error
  use moose_mpi
  use moose_geometry, only: hypersurf2d
  use flare_model,    only: bfield, boundary
  use flare_control
  use flare_tasks
  integer,          intent(in) :: iz, samples
  character(len=*), intent(in) :: match

  type(hypersurf2d), allocatable :: slice(:, :)
  character(len=256) :: msg
  logical :: partial_match
  real(real64) :: r, s, t, phi, x(2)
  integer :: i, ir, ir1, ir2, ip, ip1, ip2, it, jr, jp, jt, n, nt


  call begin_task()
  if (rank == 0) then
     select case (match)
     case ('partial')
        print 1001, iz
        partial_match = .true.

     case ('complete')
        print 1002, iz
        partial_match = .false.

     case default
        write (msg, 9001) match
        call ERROR(msg)
     end select
     print 1010, samples**3
  endif
  call proc(0)%broadcast(partial_match)
 1001 format(1x,"Collecting cells in zone ",i0," which are partially outside of boundary"/)
 1002 format(1x,"Collecting cells in zone ",i0," which are completely outside of boundary"/)
 1010 format(3x,"- Number of equidistant samples per cell: ",i0/)
 9001 format("invalid match keyword '",a,"' for SCAN_CELLS plate generator")


  ! cell resolution and index range
  ir1 = mmesh_workspace(iz)%r_surf_pl_trans_range(1)
  ir2 = mmesh_workspace(iz)%r_surf_pl_trans_range(2)-1
  ip1 = mmesh_workspace(iz)%p_surf_pl_trans_range(1)
  ip2 = mmesh_workspace(iz)%p_surf_pl_trans_range(2)-1
  nt  = mmesh_workspace(iz)%n(3)-1


  ! construct boundary slices
  allocate (slice(0:nt-1, samples))
  do it=0,nt-1
     do jt=1,samples
        t = it + (jt-0.5d0) / samples
        phi = mmesh_workspace(iz)%toroidal_angle(t)
        x = bfield%equi%magnetic_axis(phi)
        slice(it, jt) = boundary%rzslice(phi)
        call slice(it, jt)%set_inside_flag(x)
     enddo
  enddo


  ! scan cells
  i = 0
  n = (ir2-ir1+1) * (ip2-ip1+1) * nt
  call progress_bar(0,n)
  do ir=ir1,ir2
  do ip=ip1,ip2
  do it=0,nt-1
     i = i + 1;   if (mod(i,nproc) /= rank) cycle

     ! sample points in cell (ir,ip,it)
     do jr=1,samples
     do jp=1,samples
     do jt=1,samples
        r = -1.d0 + 2.d0 * (jr-0.5d0) / samples
        s = -1.d0 + 2.d0 * (jp-0.5d0) / samples
        t = it + (jt-0.5d0) / samples

        x = mmesh_workspace(iz)%rz_real_coordinates(ir, ip, r, s, t)
        ! sample point is inside boundary
        if (slice(it, jt)%includes(x)) then
           ! restore value if complete match is required
           if (.not.partial_match) workspace(iz)%plate(ir, ip, it) = 0

        ! sample point is outside boundary
        else
           ! mark cell as out of bounds
           workspace(iz)%plate(ir, ip, it) = 1
        endif
     enddo
     enddo
     enddo
     call progress_bar(i, n)
  enddo
  enddo
  enddo
  call moose_mpi_sum(workspace(iz)%plate)
  call finalize_progress_bar()
  call finalize_task()

  end subroutine scan_cells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine scan_rpaths(iz, samples, ircut, filter)
  use moose_error
  use moose_mpi
  use moose_geometry, only: hypersurf2d
  use flare_model,    only: boundary
  use flare_tasks
  integer, intent(in) :: iz, samples
  character(len=*), intent(in) :: ircut, filter

  type(hypersurf2d) :: slice
  real(real64), dimension(:,:), allocatable :: rtouch
  real(real64) :: phi, s, t, rx, x(2), x1(2), x2(2), u(1)
  integer :: i, ir, ir1, ir2, ip, ip1, ip2, it, jp, jt, nb, nt, m


  call begin_task()
  if (rank == 0) print 1000, iz
 1000 format(1x,"Scanning mesh along radial direction in zone ",i0," for intersection with plates"/)


  ir1 = mmesh_workspace(iz)%r_surf_pl_trans_range(1)
  ir2 = mmesh_workspace(iz)%r_surf_pl_trans_range(2)-1
  ip1 = mmesh_workspace(iz)%p_surf_pl_trans_range(1)
  ip2 = mmesh_workspace(iz)%p_surf_pl_trans_range(2)-1
  nt  = mmesh_workspace(iz)%n(3)-1
  allocate (rtouch(ip1:ip2,0:nt-1), source=0.d0)


  ! loop over toroidal direction
  i = 0
  m = nt * (ip2 - ip1 + 1)
  call progress_bar(0, m)
  do it=0,nt-1
     ! loop over toroidal samples per cell
     do jt=1,samples
        t = it + (jt-0.5d0) / samples
        phi = mmesh_workspace(iz)%toroidal_angle(t)

        ! construct boundary slice at this location
        slice = boundary%rzslice(phi, filter=filter)

        ! loop over poloidal direction
        do ip=ip1+rank,ip2,nproc
           i = i + nproc
           rtouch(ip,it) = ir2+1
           do jp=1,samples
              s = -1.d0 + 2.d0 * (jp-0.5d0) / samples
              x1 = mmesh_workspace(iz)%rz_real_coordinates(ir1, ip, -1.d0, s, t)

              ! scan path in radial direction for intersection with boundary
              do ir=ir1,ir2
                 x2 = mmesh_workspace(iz)%rz_real_coordinates(ir, ip, 1.d0, s, t)
                 if (slice%intersect(x1, x2, x, rx, nb, u)) then
                    rtouch(ip,it) = min(rtouch(ip,it), ir+rx)
                 endif
                 x1 = x2
              enddo
           enddo
           call progress_bar(i, m)
        enddo
     enddo
  enddo
  call moose_mpi_sum(rtouch)
  call finalize_progress_bar()


  ! save plate cells
  do it=0,nt-1
  do ip=ip1,ip2
     select case(ircut)
     case("nint")
        ir = nint(rtouch(ip,it))

     case("floor")
        ir = floor(rtouch(ip,it))

     case default
        call ERROR("invalid ircut = '"//trim(ircut)//"'")
     end select

     workspace(iz)%plate(ir:,ip,it) = 1
  enddo
  enddo
  call finalize_task()

  end subroutine scan_rpaths
  !-----------------------------------------------------------------------------

end module flare_mmesh_plates
