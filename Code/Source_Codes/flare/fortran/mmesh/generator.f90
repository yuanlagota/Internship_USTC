module flare_mmesh_generator
  use iso_fortran_env
  use flare_mmesh_layout
  implicit none


  integer :: npoints = 1024, nctrl = 0, fit_method = 0

  real(real64) :: epsabs = 1.d-7, &
     lambda1 = 0.d0,  & ! regularization parameter for non-linear fit
     lambda2 = 1.e-5


  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_flux_tubes(iz, it_base, phi)
  !
  ! Construct 3-D finite flux tubes in zone *iz* by tracing field lines from base nodes.
  !
  use moose_mpi
  use moose_utils, only: str
  use moose_grids
  use flare_tasks
  use flare_mmesh, make => construct_flux_tubes
  integer,      intent(in) :: iz
  integer,      intent(in) :: it_base
  real(real64), intent(in) :: phi(0:)

  type(mmesh) :: F
  type(qmesh) :: B
  integer :: n


  call begin_task()
  if (rank == 0) then
     print 1000, iz
     print *

     ! load base mesh
     B = qmesh("base"//str(iz)//".dat")
     print 1001, B%n, phi(it_base)
     print *

     n = ubound(phi,1)
     write (6, 1002, advance='no')
     if (it_base > 0) write (6, 1003, advance='no') phi(0)
     write (6, 1004, advance='no') phi(it_base)
     if (it_base < n) write (6, 1005, advance='no') phi(n)
     write (6, 1006)
     print *
  endif
  call B%broadcast()
 1000 format(1x,"Constructing flux tubes in zone ",i0, " ...")
 1001 format(3x,"- Base mesh with ",i0," x ",i0," nodes at ",f0.1," deg")
 1002 format(3x,"- Tracing field lines from mesh nodes:   [")
 1003 format(f0.1," <- ")
 1004 format(f0.1)
 1005 format(" -> ",f0.1)
 1006 format("]")


  F = make(B, it_base, phi)
  if (rank == 0) call F%savetxt(mmesh_filename(iz))
  call finalize_task()

  end subroutine construct_flux_tubes
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine construct_n0_domain(iz, iblock, ilayer)
  use moose_mpi
  use flare_mmesh
  use flare_mmesh_parameters
  use flare_mmesh_core, core_mesh_generator => mesh_generator
  use flare_mmesh_vacuum
  integer, intent(in) :: iz, iblock, ilayer

  integer :: iside, nr1, nr2



  ! lower radial boundary
  nr1 = 0
  select case(connectZ(ilayer)%radial(1))
  case(CORE)
     nr1 = nr_core
  case(VACUUM)
     nr1 = nr_vac(ilayer)
  end select

  ! upper radial boundary
  nr2 = 0
  select case(connectZ(ilayer)%radial(2))
  case(VACUUM)
     nr2 = nr_vac(ilayer)
  end select


  ! load flux tubes
  if (rank == 0) then
     workspace(iz) = load_flux_tubes(mmesh_filename(iz), [nr1, nr2])
  endif
  call workspace(iz)%broadcast()


  ! construct core domain (if applicable)
  if (nr_core > 0  .and.  connectZ(ilayer)%radial(1) == CORE) then
     call core_mesh_generator(iz, core_domain)
  endif


  ! construct vacuum domain (if applicable)
  if (nr_vac(ilayer) > 0) then
  do iside=1,2
     if (connectZ(ilayer)%radial(iside) == VACUUM) then
        call mesh_generator(iz, iside, vacuum_domain(layer_index(iz)), Bmod_in_vacuum_domain)
     endif
  enddo
  do iside=1,2
     if (connectZ(ilayer)%poloidal(iside) == PLATE) then
        if (iside == 1  .and.  closure_R(ilayer) /= "") then
           call close_mesh_domain_usr(iz, iside, closure_R(ilayer))
        elseif (iside == 2  .and.  closure_L(ilayer) /= "") then
           call close_mesh_domain_usr(iz, iside, closure_L(ilayer))
        else
           call close_mesh_domain(iz, iside)
        endif
     endif
  enddo
  endif


  end subroutine construct_n0_domain
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_input_files()
  use flare_mmesh
  use flare_mmesh_parameters
  use flare_mmesh_plates, only: PLATES_FILENAME

  integer, allocatable :: &
     ZON_RADI(:), ZON_POLO(:), ZON_TORO(:), &
     SRF_RADI(:), SRF_POLO(:), SRF_TORO(:)
  integer :: iu, iz


  allocate (ZON_RADI(0:nz-1), ZON_POLO(0:nz-1), ZON_TORO(0:nz-1))
  allocate (SRF_RADI(0:nz-1), SRF_POLO(0:nz-1), SRF_TORO(0:nz-1))
  do iz=0,nz-1
     SRF_RADI(iz) = workspace(iz)%n(1)
     SRF_POLO(iz) = workspace(iz)%n(2)
     SRF_TORO(iz) = workspace(iz)%n(3)
  enddo
  ZON_RADI = SRF_RADI - 1
  ZON_POLO = SRF_POLO - 1
  ZON_TORO = SRF_TORO - 1


  write (6, 4000, advance='no') 'input.geo';   call write_input_geo()
  write (6, 4000, advance='no') 'input.N0G';   call write_input_n0g()
  if (generate_add_sf_n0) call write_ADD_SF_N0()
  deallocate (ZON_RADI, ZON_POLO, ZON_TORO, SRF_RADI, SRF_POLO, SRF_TORO)
  print *
 4000 format(3x,a)


  open  (newunit=iu, file='.geometry')
  write (iu, 1000)
  write (iu, 1001) "layout", LAYOUT_FILENAME
  write (iu, 1001) "grid3d", GRID3D_FILENAME
  write (iu, 1001) "n0g",    "input.N0G"
  write (iu, 1001) "plates", PLATES_FILENAME
  write (iu, 1001) "bfield", BFIELD_FILENAME
  close (iu)
 1000 format("[geometry]")
 1001 format(a,': "',a,'"')

  contains
  !.............................................................................
  subroutine write_input_geo()

  real(real64) :: Atot, A
  integer :: ir, ip, it, iz, irun, iside, iu, itype, n


  open  (newunit=iu, file='input.geo')
  write (iu, 1000)
  write (iu, 1010) nsymmetry(symmetry, stellarator_symmetry)
 1000 format ('* geometry information for EMC3')
 1010 format ('* SYMMETRY ',i0)

  ! 1. geometry, mesh resolution
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003) nz
  write (iu, 1004)
  do iz=0,nz-1
     write (iu, *) SRF_RADI(iz), SRF_POLO(iz), SRF_TORO(iz)
  enddo
 1001 format ('*** 1. grid resolution')
 1002 format ('* number of zones/blocks')
 1003 format (i0)
 1004 format ('* number of radial, poloidal and toroidal grid points')


  ! 2. surface definitions
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  ! 2.1 non default surfaces (periodic, mapping, ...)
  write (iu, 2001)
  ! 2.1.a - radial
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
        ir = 0
        do iside=1,2
           itype = connectZ(layer_index(iz))%radial(iside)
           if (itype > 0) then
              if (irun == 0) then
                 n = n + 1
              else
                 write (iu, *) ir, iz, itype
                 write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
              endif
           endif
           ir = SRF_RADI(iz)-1
        enddo
     enddo
  enddo
  ! 2.1.b - poloidal
  write (iu, 2003)
  n = 0
  do irun=0,1
     ! write number of non default poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
        ip = 0
        do iside=1,2
           itype = connectZ(layer_index(iz))%poloidal(iside)
           if (itype > 0) then
              if (irun == 0) then
                 n = n + 1
              else
                 write (iu, *) ip, iz, itype
                 write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
              endif
           endif
           ip = SRF_POLO(iz)-1
        enddo
     enddo
  enddo
  ! 2.1.c - toroidal
  write (iu, 2004)
  n = 0
  do irun=0,1
     ! write number of non default toroidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
        it = 0
        do iside=1,2
           itype = connectZ(layer_index(iz))%toroidal(iside)
           if (itype > 0) then
              if (itype == MAPPING) itype = toroidal_mapping_type(block_index(iz), iside)
              if (irun == 0) then
                 n = n + 1
              else
                 write (iu, *) it, iz, itype
                 write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_POLO(iz)-1
              endif
           endif
           it = SRF_TORO(iz)-1
        enddo
     enddo
  enddo

  ! 2.2 non transparent surfaces (boundary conditions)
  Atot = 0.d0
  write (iu, 2005)
  ! 2.2.a - radial
  write (iu, 2002)
  do irun=0,1
     ! write number of non transparent radial surfaces
     if (irun == 1) write (iu, *) n
     n = 0

     do iz=0,nz-1
     do iside=1,2
        itype = connectZ(layer_index(iz))%radial(iside)
        if (itype == CORE) then
           A = workspace(iz)%phi(ZON_TORO(iz)) - workspace(iz)%phi(0)
           n = n + 1
           if (irun == 0) then
              Atot = Atot + A
           else
              write (iu, 2010) n, A/Atot
           endif
        endif
        if (itype == VACUUM) then
           n = n + 1
           if (irun == 1) write (iu, 2011) n
        endif


        if (irun == 1  .and.  itype < 0) then
           write (iu, *) workspace(iz)%r_surf_pl_trans_range(iside), iz, (-1)**(iside+1)
           write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
        endif
     enddo
     enddo
  enddo
  ! 2.2.b - poloidal
  write (iu, 2003)
  n = 0
  do irun=0,1
     ! write number of non transparent poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
     do iside=1,2
        itype = connectZ(layer_index(iz))%poloidal(iside)
        if (itype < 0) then
           if (irun == 0) then
              n = n + 1
           else
              write (iu, *) workspace(iz)%p_surf_pl_trans_range(iside), iz, (-1)**(iside+1)
              write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
           endif
        endif
     enddo
     enddo
  enddo
  ! 2.2.c - toroidal
  write (iu, 2004)
  write (iu, *) 0

  ! 2.3 plate surfaces
  write (iu, 2006)
  write (iu, 2002)
  write (iu, *) -plate_format ! user defined
  write (iu, 2003)
  write (iu, *) -plate_format ! user defined
  write (iu, 2004)
  write (iu, *) -plate_format ! user defined
 2000 format ('*** 2. surface definitions')
 2001 format ('*** 2.1 non default surface')
 2002 format ('* radial')
 2003 format ('* poloidal')
 2004 format ('* toroidal')
 2005 format ('*** 2.2 non transparent surface (Boundary condition must be defined)')
 2006 format ('*** 2.3 plate surface (Bohm Boundary condition)')
 2010 format ('* ',i0,': CORE BOUNDARY 'e11.5)
 2011 format ('* ',i0,': VACUUM BOUNDARY ')


  ! 3. physical cell definition
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu,    *) cell_def
  write (iu, 3003) (cell_param, iz=1,nz)
  write (iu, 3001)
  write (iu, 3002) .true.
 3000 format ('*** 3. physical cell definition')
 3001 format ('* run cell check?')
 3002 format (L1)
 3003 format (99(2x,i0))
  close (iu)
 9999 format ('*',32('-'))

  end subroutine write_input_geo
  !.............................................................................
  subroutine write_input_n0g()

  integer :: N0_DENS(4), DIA_SFS
  integer :: dp, dt, il, ir, ip, iz, irun, iside, itype, iu, n


  open  (newunit=iu, file='input.N0G')
  write (iu, 1000)
  write (iu, 9999)
  write (iu, 1001)
  write (iu, 9999)
  write (iu, 1002)
  write (iu, 1003)
  write (iu, 1004)
 1000 format ('******** additional geometry and parameters for EIRENE ****')
 1001 format ('*** 1. non-transparent surfaces for neutral particles')
 1002 format ('*  non-transparent surfaces with informations about')
 1003 format ('*  this surface being defined in EIRENE. The surface')
 1004 format ('*  number must be indicated here.')

  ! 1. non-transparent surfaces
  ! 1.1 radial
  write (iu, 1012)
  n = 0
  do irun=0,1
     ! write number of non default radial surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
        ir = 0
        do iside=1,2
           itype = connectZ(layer_index(iz))%radial(iside)
           if (itype < 0) then
              if (irun == 0) then
                 n = n + 1
              else
                 write (iu, 1015) ir, iz, -EIRENE_SF_NUM(-itype)
                 write (iu, *) 0, ZON_POLO(iz)-1, 0, ZON_TORO(iz)-1
              endif
           endif
           ir = SRF_RADI(iz)-1
        enddo
     enddo
  enddo
  ! 1.2. poloidal
  write (iu, 1013)
  n = 0
  do irun=0,1
     ! write number of non transparent poloidal surfaces
     if (irun == 1) write (iu, *) n

     do iz=0,nz-1
        ip = 0
        do iside=1,2
           itype = connectZ(layer_index(iz))%poloidal(iside)
           if (itype < 0) then
              if (irun == 0) then
                 n = n + 1
              else
                 write (iu, 1015) ip, iz, -EIRENE_SF_NUM(-itype)
                 write (iu, *) 0, ZON_RADI(iz)-1, 0, ZON_TORO(iz)-1
              endif
           endif
           ip = SRF_POLO(iz)-1
        enddo
     enddo
  enddo
  ! 1.3. toroidal
  write (iu, 1014)
  write (iu, *) 0
 1012 format ('* radial')
 1013 format ('* poloidal')
 1014 format ('* toroidal')
 1015 format (3i8)


  ! 2. additional physical cells for neutrals
  write (iu, 9999)
  write (iu, 2000)
  write (iu, 9999)
  write (iu, 2001)
  write (iu, 2002)
  n = 0
  do irun=0,1
     ! write number of additional cell blocks
     if (irun == 1) write (iu, *) n, 70

     ! confined region
     do iz=0,nz-1
        if (connectZ(layer_index(iz))%radial(1) == CORE) then
           if (irun == 0) then
              n = n + workspace(iz)%r_surf_pl_trans_range(1)
           else
              dp = dp_core;   if (dp < 0) dp = ZON_POLO(iz)
              dt = dt_core;   if (dt < 0) dt = ZON_TORO(iz)
              do ir=0,workspace(iz)%r_surf_pl_trans_range(1)-1
                 write (iu, *) 2, 'EIRENE_CORE_MODEL'
                 write (iu, 2003) iz, ir, ir+1, 1, 0, ZON_POLO(iz), dp, 0, ZON_TORO(iz), dt
                 write (iu, 2004) iz, ir
              enddo
           endif
        endif
     enddo

     ! vacuum region (collect what's left)
     do iz=0,nz-1
        if (irun == 0) then
           n = n + 1
        else
           il = layer_index(iz)
           dp = dp_vac(il);   if (dp < 0) dp = ZON_POLO(iz)
           dt = dt_vac(il);   if (dt < 0) dt = ZON_TORO(iz)
           write (iu, *) 2, 0
           write (iu, 2003) iz, 0, ZON_RADI(iz), 1, 0, ZON_POLO(iz), dp, 0, ZON_TORO(iz), dt
           write (iu, 2005) iz
        endif
     enddo
  enddo
 2000 format ('*** 2. DEFINE ADDITIONAL PHYSICAL CELLS FOR NEUTRALS')
 2001 format ('*   ZONE  R1    R2    DR    P1    P2    DP    T1    T2    DT')
 2002 format ('* ne       Te      Ti        M')
 2003 format (10i6)
 2004 format ('EIRENE_CORE_DATA_',i0,'_',i0)
 2005 format ('EIRENE_VACUUM_DATA_',i0)


  ! 3. Neutral sources
  write (iu, 9999)
  write (iu, 3000)
  write (iu, 9999)
  write (iu, 3001)
  write (iu, *) 0, 0, nsside
 3000 format ('*** 3. Neutral Source distribution')
 3001 format ('* N0S NS_PLACE  NSSIDE')


  ! 4. Additional surfaces
  write (iu, 9999)
  write (iu, 4000)
  write (iu, 9999)
  write (iu, 4001)
 4000 format ('*** 4. Additional surfaces')
 4001 format ('BOUNDARY/ADD_SF_N0')


  ! 5. Diagnostics
  write (iu, 9999)
  write (iu, 5000)
  write (iu, 5010)
  N0_DENS = 0
  write (iu, *) N0_DENS
  write (iu, 5020)
  DIA_SFS = 0
  write (iu, *) DIA_SFS
  close (iu)
 5000 format ('*** 5. Neutral gas diagnostics')
 5010 format ('*** 5.1. Particle and energy densities for atoms and molecules')
 5020 format ('*** 5.2. Flux and spectrum on a given surface')
 9999 format ('*',32('-'))

  end subroutine write_input_n0g
  !.............................................................................
  subroutine write_ADD_SF_N0()
  use moose_utils,    only: isdir, join, mkdir, make_filename, endswith
  use moose_geometry
  use flare_model,    only: boundary, bfield

  type(polygon2d)    :: P
  character(len=256) :: filename
  real(real64) :: phi, x0(2)
  integer :: i, isgn, iu, n


  ! check if path for BOUNDARY already exists
  if (isdir("BOUNDARY")) then
     print *
     print *, "BOUNDARY directory already exists, skipping setup of boundary files"
     return
  endif
  write (6, 4000, advance='no') "BOUNDARY/ADD_SF_N0"
  print *
  call mkdir("BOUNDARY")
 4000 format(3x,a)


  ! initialize layout file
  n = 0
  do i=1,boundary%nsurfaces
     if (endswith(boundary%surfaces(i)%key, ".T") .and. stellarator_symmetry) exit
     n = n + 1
  enddo
  open  (newunit=iu, file="BOUNDARY/ADD_SF_N0")
  write (iu, *) n


  ! loop over surfaces
  do i=1,boundary%nsurfaces
  select type (S => boundary%surfaces(i)%geometry)
  type is (axisurf)
     ! axisymmetric surfaces
     filename = make_filename(boundary%surfaces(i)%key)
     write (iu, 1000) trim(filename)
     call export_axisurf(S, join("BOUNDARY", filename))

     ! check if orientation is consistent with nsside
     if (S%P%is_closed()) then
        isgn = 1
        x0 = bfield%equi%magnetic_axis(0.d0)
        if (S%P%winding_number(x0) == 0) isgn = -1
        if (S%P%orientation() * nsside * isgn > 0) print 9001, filename
     endif

  type is (torosurf)
     ! non-axisymmetric surfaces
     filename = make_filename(boundary%surfaces(i)%key)
     write (iu, 1000) trim(filename)
     call S%savetxt(join("BOUNDARY", filename), iotype='legacy', scale=1.d2)

     ! check if orientation is consistent with nsside
     P = S%polygon2d(S%phi(0))
     if (P%is_closed()) then
        isgn = 1
        x0 = bfield%equi%magnetic_axis(S%phi(0))
        if (P%winding_number(x0) == 0) isgn = -1
        if (P%orientation() * nsside * isgn < 0) print 9001, trim(filename)
     endif
  end select
  enddo
  close (iu)
 1000 format("   0  -4   1",/,a)
 9001 format("WARNING: orientation of ",a," appears to be inconsistent with choice of 'nsside'")

  end subroutine write_ADD_SF_N0
  !.............................................................................
  subroutine export_axisurf(A, filename)
  use moose_geometry
  class(axisurf),   intent(in) :: A
  character(len=*), intent(in) :: filename

  integer :: i, iu, n


  n = A%P%segments()

  open  (newunit=iu, file=filename)
  write (iu, 1000) A%description()
  write (iu, *) 1, n+1, 1, 0.d0, 0.d0
  write (iu, *) -180.d0
  do i=0,n
     write (iu, *) A%P%node(i) * 1.d2 ! m -> cm
  enddo
  close (iu)
 1000 format("# ",a)

  end subroutine export_axisurf
  !.............................................................................
  end subroutine write_input_files
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine evaluate_psiN()
  use moose_mpi
  use flare_equi2d, only: equi2d
  use flare_model,  only: bfield
  use flare_mmesh
  use flare_mmesh_parameters, only: Bmod_in_vacuum_domain

  character(len=*), parameter :: PSIN_FILENAME   = "psiN.dat"

  real(real64), allocatable :: psiN(:,:,:)
  real(real64) :: r(3)
  integer :: i, iz, ir, ir1, ir2, ip, ip1, ip2, it, iu, n(3)


  if (rank == 0) open  (newunit=iu, file=PSIN_FILENAME)
  do iz=0,size(workspace)-1
     n = workspace(iz)%n
     allocate (psiN(0:n(1)-1, 0:n(2)-1, 0:n(3)-1), source=0.d0)

     ir1 = workspace(iz)%r_surf_pl_trans_range(1)
     ir2 = workspace(iz)%r_surf_pl_trans_range(2)
     ip1 = workspace(iz)%p_surf_pl_trans_range(1)
     ip2 = workspace(iz)%p_surf_pl_trans_range(2)
     if (Bmod_in_vacuum_domain) then
        ir1 = 0
        ir2 = n(1) - 1
        ip1 = 0
        ip2 = n(2) - 1
     endif

     i = 0
     do ir=ir1,ir2
     do ip=ip1,ip2
     do it=0,n(3)-1
        i = i + 1;   if (mod(i, nproc) /= rank) cycle
        r = [workspace(iz)%r(ir,ip,it), workspace(iz)%z(ir,ip,it), workspace(iz)%phi(it)]
        psiN(ir,ip,it) = bfield%equi%psiN(r)
     enddo
     enddo
     enddo
     call moose_mpi_sum(psiN)

     if (rank == 0) write (iu, *) psiN
     deallocate (psiN)
  enddo
  if (rank == 0) close (iu)

  end subroutine evaluate_psiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine apply_mmesh_updown_symmetry(iz)
  use moose_error
  use moose_grids, only: qmesh
  use flare_model, only: bfield
  use flare_mmesh_parameters
  use flare_mmesh_layout
  use flare_mmesh
  use flare_mmesh_base_generator, only: updown_symmetry_tolerance
  integer, intent(in) :: iz

  type(qmesh)  :: M
  real(real64) :: x0(2), b0(3), r, z, b, dr, dz, db
  integer :: ir, ip, jp, it0


  it0 = -1
  if (updown_symmetry < 0  .and.  block_index(iz) == 0) it0 = 0
  if (updown_symmetry > 0  .and.  block_index(iz) == blocks-1) it0 = workspace(iz)%n(3)-1
  if (it0 == -1) return


  associate (this => workspace(iz))
  x0 = bfield%equi%magnetic_axis(this%phi(it0))
  b0 = bfield%eval([x0, this%phi(it0)])
  do ir=0,this%n(1)-1
  do ip=0,(this%n(2)-1)/2
     jp = this%n(2)-1-ip

     r  = (this%r(ir,ip,it0) + this%r(ir,jp,it0)) / 2
     z  = (this%z(ir,ip,it0) - this%z(ir,jp,it0)) / 2
     b  = (this%b(ir,ip,it0) + this%b(ir,jp,it0)) / 2

     dr = this%r(ir,ip,it0) - r
     dz = this%z(ir,ip,it0) - z
     db = this%b(ir,ip,it0) - b
     if (max(sqrt(dr**2 + dz**2) / x0(1), abs(db) / b0(3)) > updown_symmetry_tolerance) then
        print *, "iz         = ", iz
        print *, "ir, ip, jp = ", ir, ip, jp
        print *, "x(ir,ip)   = ", this%r(ir,ip,it0), this%z(ir,ip,it0)
        print *, "x(ir,jp)   = ", this%r(ir,jp,it0), this%z(ir,jp,it0)
        print *, "dr         = ", dr
        print *, "dz         = ", dz
        print *, "db         = ", db
        M = qmesh(this%r(:,:,it0), this%z(:,:,it0))
        call M%savetxt("ERROR_UPDOWN_SYMMETRY")
        call ERROR("up/down symmetry exceeds tolerance")
     endif

     this%r(ir,ip,it0) = r;   this%z(ir,ip,it0) =  z;   this%b(ir,ip,it0) = b
     this%r(ir,jp,it0) = r;   this%z(ir,jp,it0) = -z;   this%b(ir,jp,it0) = b
  enddo
  enddo
  end associate

  end subroutine apply_mmesh_updown_symmetry
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_plates()
  use moose_mpi
  use moose_error
  use moose_utils, only: str, substring, nsubstrings
  use flare_mmesh_parameters
  use flare_mmesh_plates

  character(len=256) :: instruction, cmd, match, ircut, filter, msg
  integer :: i, iostat, iz, n, samples


  if (plate_generator == "") call ERROR("undefined plate generator")
  call init_workspace()


  n = nsubstrings(plate_generator)
  do i=1,n
     instruction = substring(plate_generator, i)
     if (instruction == "") cycle

     read (instruction, *) cmd
     select case (cmd)
     ! scan cells for points outside of poloidally closed boundary
     case ('SCAN_CELLS')
        read (instruction, *, iostat=iostat) cmd, match, samples
        if (iostat /= 0) call ERROR("invalid format for SCAN_CELLS parameters")
        if (samples <= 0) call ERROR("samples > 0 required for SCAN_CELLS plate generator")

        do iz=0,size(workspace)-1
           call scan_cells(iz, match, samples)
        enddo


     ! scan radial paths through mesh
     case ('SCAN_RPATHS')
        read (instruction, *, iostat=iostat) cmd, samples, ircut
        if (iostat /= 0) call ERROR("invalid format for SCAN_RPATHS parameters")
        if (samples <= 0) call ERROR("samples > 0 required for SCAN_RPATHS plate generator")

        read (instruction, *, iostat=iostat) cmd, samples, ircut, filter
        if (iostat /= 0) filter = ""

        do iz=0,size(workspace)-1
           call scan_rpaths(iz, samples, ircut, filter)
           call exclude_cells_from_flux_tubes(iz, 2)
        enddo


     ! invalid plate generator
     case default
        write (msg, 9001) trim(instruction)
        call ERROR(msg)
     end select
  enddo
 9001 format("invalid plate generator ",a)


  if (rank == 0) call save_plates(plate_format)
  call free_workspace()

  end subroutine generate_plates
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine generate_reservoirs()
  use moose_mpi
  use moose_error
  use flare_model, only: equi2d
  use flare_tasks
  use flare_mmesh_parameters, only: layout, tblock
  use flare_mmesh
  use base_mesh

  real(real64) :: r, z, x0(2), x1(2), x95(2)
  integer, allocatable :: id(:,:,:)
  integer :: ir, ip, it, iu, iz, nr, np, npX, nt


  call begin_task()
  if (rank == 0) then
     print 1000
     x0 = equi2d%r0
     x1 = equi2d%xpoint(1)
     x95 = x1 + (x0 - x1) * 0.05d0
     open  (newunit=iu, file="DIVERTOR_RESERVOIRS")
  endif
  call proc(0)%broadcast(x95)
  if (.not.layout == TOPO_LSN) call ERROR("lsn")
 1000 format(1x,"Generating reservoirs for divertor analysis ...")


  do iz=0,size(workspace)-1
     nr = size(workspace(iz)%r, 1) - 1
     np = size(workspace(iz)%r, 2) - 1
     nt = size(workspace(iz)%r, 3) - 1
     print *, iz, nr, np, nt
     allocate (id(0:nr-1, 0:np-1, 0:nt-1), source=0)

     select case(mod(iz,3))
     ! core
     case(0)

     ! SOL
     case(1)
        do it=0,nt-1
        do ir=0,nr-1
           ! z = z95 line
           do ip=0,np/2
              z = sum(workspace(iz)%z(ir:ir+1, ip:ip+1, it:it+1)) / 8
              id(ir,ip,it) = 1
              if (z > x95(2)) exit
           enddo

           ! r - r95 = z95 - z line (diagonal through r95, z95)
           do ip=np-1,np/2,-1
              r = sum(workspace(iz)%r(ir:ir+1, ip-1:ip, it:it+1)) / 8
              z = sum(workspace(iz)%z(ir:ir+1, ip-1:ip, it:it+1)) / 8
              id(ir,ip,it) = 2
              if (r - x95(1) + z - x95(2) > 0.d0) exit
           enddo
        enddo
        enddo

     ! private flux region
     case(2)
        npX = tblock(iz/3)%npR(1)
        do it=0,nt-1
        do ir=0,nr-1
           id(ir,0:npX-1,it) = 1
           id(ir,npX:np-1,it) = 2
        enddo
        enddo

     end select

     call moose_mpi_sum(id)
     if (rank == 0) write (iu, '(i8)') id
     deallocate (id)
  enddo


  if (rank == 0) close (iu)
  call finalize_task()

  end subroutine generate_reservoirs
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine exec(mmesh_parameters, subtasks)
  !
  ! Magnetic mesh generator for field line reconstruction.
  !
  use moose_error
  use moose_mpi
  use flare_tasks
  use flare_mmesh, make => construct_flux_tubes
  use flare_mmesh_parameters
  use flare_mmesh_inner_boundary
  use flare_mmesh_base_generator
  use flare_mmesh_unstructured_generator, generate_unstructured_mmesh => generate_mmesh
  character(len=*), intent(in) :: mmesh_parameters
  logical,          intent(in) :: subtasks(6)

  real(real64), allocatable :: phi_base(:)
  integer :: iblock, ilayer, iz


  if (rank == 0) then
     print '(80("="))'
     print *, "Starting magnetic mesh generator for field line reconstruction..."
     print *
     call load_mmesh_parameters(mmesh_parameters)
  endif
  call broadcast_mmesh_parameters()
  call setup_mmesh_parameters()
  call init_layout()


  ! Task 1: generate innermost boundaries
  if (subtasks(1)) then
     ! 1.1. construct Poincare maps
     if (default_layout) then
        call make_points_default(p1, p2, nsymmetry(symmetry, stellarator_symmetry), phi0, blocks, npoints)

     else
        allocate (phi_base(0:blocks-1))
        do iblock=0,blocks-1
           phi_base(iblock) = T(iblock)%phi_base
        enddo
        call make_points_usr(p1, p2, symmetry, phi_base, npoints)
        deallocate (phi_base)
     endif
     call mpi_barrier_world()

     ! 1.2. fit curve
     call bspline_multifit(blocks, 4, nctrl, epsabs, fit_method, lambda1, lambda2)
  endif
  call mpi_barrier_world()


  ! Task 2: construct base mesh(s)
  if (subtasks(2)) then
     if (rank == 0) then
        call init_base_mesh_generator()
        do iblock=0,blocks-1
           call generate_base_mesh(iblock)
        enddo
     endif
  endif
  call mpi_barrier_world()


  ! Task 3: construct 3-D finite flux tubes
  if (subtasks(3)) then
     if (layout == LAYOUT_UNSTRUCTURED) then
        if (nproc > 1) call ERROR("parallelization of unstructured mesh generation not implemented yet")
        call generate_unstructured_mmesh(curve_filename(0, 0), curve_filename(1, 0), &
           nsymmetry(symmetry, stellarator_symmetry), T(0)%phi, T(0)%it_base, tblock(0)%np(0), &
           delta_r, rinc, divmax, divavg)

     else
        iz = 0
        do iblock=0,blocks-1
           do ilayer=0,layers-1
              call construct_flux_tubes(iz, T(iblock)%it_base, T(iblock)%phi)
              iz = iz + 1
           enddo
        enddo
     endif
  endif
  call mpi_barrier_world()


  ! Task 4: construct neutral gas domain (+ extension into core)
  if (subtasks(4)) then
     call init_mmesh_workspace(nz)

     iz = 0
     do iblock=0,blocks-1
        do ilayer=0,layers-1
           call construct_n0_domain(iz, iblock, ilayer)
           if (updown_symmetry /= 0) call apply_mmesh_updown_symmetry(iz)
           iz = iz + 1
        enddo
     enddo

     if (rank == 0) then
        call save_mmesh()
        call write_input_files()
     endif
     call evaluate_psiN()
     if (.not.(subtasks(5).or.subtasks(6))) call free_mmesh_workspace()
  endif
  call mpi_barrier_world()


  ! Task 5: construct plate surface representation in mmesh
  if (subtasks(5)) then
     if (.not.subtasks(4)) call load_mmesh()

     call generate_plates()

     if (.not.subtasks(6)) call free_mmesh_workspace()
  endif


  ! Task 6: reservoirs for divertor analysis
  if (subtasks(6)) then
     if (.not.subtasks(4)) call load_mmesh()

     call generate_reservoirs()

     call free_mmesh_workspace()
   endif

  ! cleanup
  call free_layout()

  end subroutine exec
  !-----------------------------------------------------------------------------

end module flare_mmesh_generator
