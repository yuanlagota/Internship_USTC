subroutine flare_task_fluxsurf3d_grid(r0, nsym, nphi, ntheta, endpoints, output)
  !
  ! Create 3D grid on flux surface from Poincare map for `r0`. The flux surface must be closed, i.e. must not intersect the boundary.
  !
  ! **Parameters:**
  !
  ! :r0:         Reference point on flux surface (r[m], z[m], phi[deg]).
  !
  ! :nsym:       Toroidal symmetry of flux surface.
  !
  ! :nphi:       Number of steps in toroidal direction.
  !
  ! :ntheta:     Number of steps in poloidal direction.
  !
  ! :endpoints:  Include/exclude periodic endpoints of grid.
  !
  ! :output:     Filename for output of grid.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math, only: pi2, CYLINDRICAL_COORDINATES
  use moose_grids
  use moose_units
  use moose_r3grid
  use flare_control
  use flare_fluxsurf3d
  use flare_tasks
  implicit none
  real(real64),     intent(in) :: r0(3)
  integer,          intent(in) :: nsym, nphi, ntheta
  logical,          intent(in) :: endpoints
  character(len=*), intent(in) :: output

  type(fluxsurf3d) :: F
  type(tpzmesh3d)  :: M
  type(r3grid)     :: R3
  real(real64)     :: r0_rad(3), theta
  integer :: i, j, nu, nv


  call begin_task()
  if (rank == 0) then
     print *, "Constructing grid on flux surface ..."
     print *


     ! construct flux surface
     r0_rad(1:2) = r0(1:2);   r0_rad(3) = r0(3) / 360 * pi2
     F = fluxsurf3d(r0_rad, nsym, nphi, report=.true.)
     print *
     if (verbose) call F%save("fluxsurf3d")


     ! sample points on flux surface for grid
     print 1001, nphi, ntheta, nsym
     nu = nphi
     nv = ntheta
     if (endpoints) then
        nu = nu + 1
        nv = nv + 1
     endif
     M = tpzmesh3d(nu, nv, "Toroidal angle [deg]", "Poloidal angle [deg]", "Major radius [m]", "Z [m]")
     do i=0,nphi-1
        M%domain%u(i) = F%section(i)%phiX / pi2 * 360.d0
        do j=0,ntheta-1
           theta = 360.d0 * j / ntheta
           M%domain%v(i,j) = theta
           M%x(:,i,j)      = F%slice(i)%eval(j, ntheta-1)
        enddo
        if (endpoints) then
           M%domain%v(i,ntheta) = 360.d0
           M%x(:,i,ntheta)      = M%x(:,i,0)
        endif
     enddo
     if (endpoints) then
        M%domain%u(nphi)   = 360.d0 / nsym
        M%domain%v(nphi,:) = M%domain%v(0,:)
        M%x(:,nphi,:)      = M%x(:,0,:)
     endif

     ! wrap Tpzmesh3d in R3grid
     R3 = cylindrical_r3grid(M, METER, DEGREE, 1, 2, 3, .false.)
     call R3%savetxt(output)
  endif
 1001 format(3x,"- Sampling ",i0," x ",i0," mesh nodes with toroidal symmetry ",i0/)


  call finalize_task()

end subroutine flare_task_fluxsurf3d_grid
