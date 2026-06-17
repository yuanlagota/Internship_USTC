module flare_mmesh_parameters
  use iso_fortran_env
  implicit none


  character(len=*), parameter :: &
     LAYOUT_GENERIC      =  "generic", &
     LAYOUT_LSN          =  "lsn", &
     LAYOUT_DDN          =  "ddn", &
     LAYOUT_CDN          =  "cdn", &
     LAYOUT_UNSTRUCTURED =  "unstructured"


  ! internal parameters
  integer, parameter :: &
     max_blocks = 128, &
     max_layers = 128



  ! user defined parameters ....................................................
  character(len=256) :: &
     layout                              = "", &  ! zone layout (lsn, ddn, ...)
     pcoordinates                        = "", &
     radial_spacing(-1:max_layers-1)     = "", &
     poloidal_spacing(0:max_layers-1)    = "", &
     poloidal_spacing_L(0:max_layers-1)  = "", &
     poloidal_spacing_R(0:max_layers-1)  = "", &
     guiding_contour(-1:max_blocks-1)    = "", &
     guiding_contour_L(-1:max_blocks-1)  = "", &
     guiding_contour_R(-1:max_blocks-1)  = "", &
     upstream_adjust_L(0:max_blocks-1)   = "", &
     upstream_adjust_R(0:max_blocks-1)   = "", &
     qmesh_generator                     = "", &
     core_domain                         = "", &
     vacuum_domain(0:max_layers-1)       = "", &
     closure_L(0:max_layers-1)           = "", &
     closure_R(0:max_layers-1)           = "", &
     plate_generator                     = ""


  integer :: &
     symmetry                       =   1, &  ! toroidal symmetry (i.e. 5 => 72 deg)
     updown_symmetry                =   0, &  ! symmetric mesh at upper or lower toroidal boundary
     blocks                         =   1, &  ! number of toroidal blocks in 360 deg / symmetry
     nt                             =  12, &  ! default toroidal resolution
     np(0:max_layers-1)             = 360, &  ! default poloidal resolution
     npL(0:max_layers-1)            =  30, &  ! default poloidal resolution in divertor legs
     npR(0:max_layers-1)            =  30, &  !    (L)eft and (R)ight segments
     nr(0:max_layers-1)             =  32, &  ! default radial resolution
     n_interpolate                  =   4, &  ! number of interpolated flux surfaces (for the transition between the perturbed flux surfaces at the inner simulation boundary and unperturbed flux surfaces further outside)
     npXqo                          =  -1, &  ! number of quasi-orthogonal surfaces below X-point
     npSP_subres                    =   1, &  ! sub-resolution of strike point alignment
     npSP_extend                    =   0, &  ! extension of strike point alignment beyond zone
     nr_core                        =   1, &  ! radial resolution in core (EIRENE only)
     dp_core                        =  -1, &  ! cell increment in poloidal direction
     dt_core                        =  -1, &  ! cell increment in toroidal direction
     nr_vac(0:max_layers-1)         =   1, &  ! radial resolution in vacuum (EIRENE only)
     dp_vac(0:max_layers-1)         =  -1, &  ! cell increment in poloidal direction
     dt_vac(0:max_layers-1)         =  -1, &  ! cell increment in toroidal direction
     nr_perturbed                   =   2, &  ! number of perturbed flux surfaces at the inner boundary
     cell_def                       =  -1, &  ! physical cell definition method
     cell_param                     =   1, &  ! parameter for physical cell definition
     nsside                         =   1, &  ! flag for plasma side of boundary surface
     plate_format                   =   1, &  ! format for plate definition file
     EIRENE_SF_NUM(3)               =   0, &  ! non-default std. surface # in input.eir (block 3A)
     EIRENE_CORE_SF_NUM             =   3, &  ! for innermost (core) surface
     EIRENE_VAC_SF_NUM              =   2     ! for outermost (vacuum) surface


  real(real64) :: &
     phi0                           = -360.d0, &  ! lower boundary of simulation domain [deg]
     p1(3)                          = [120.d0, 0.d0, 0.d0], &  ! reference points
     p2(3)                          = [119.d0, 0.d0, 0.d0], &  ! ... on 1st and 2nd grid surface
     d_SOL(2)                       =   24.d0, &  ! radial width of scrape-off layer
     d_PFR(4)                       =   15.d0, &  ! radial width of private flux region
     polo_ext                       =    0.d0, &  ! poloidal extension for divertor leg closure
     delta_r                        =    1.d-2, & ! radial increment [cm] for unstructured mesh
     rinc                           =    0.d0, &  ! relative increase or delta_r per layer
     divmax                         =    0.3d0, & ! max. allowed flux conservation error (for unstructured mmesh)
     divavg                         =    0.3d0    ! avg. allowed flux conservation error (for unstructured mmesh)


  logical :: &
     auto_adjust_mesh               =  .false., &
     Bmod_in_vacuum_domain          =  .false., &
     stellarator_symmetry           =  .false., &
     generate_add_sf_n0             =  .false.


  ! block specific parameters
  type :: tblock_parameters
     integer :: &
        nr(0:max_layers-1)   =  -1, &  ! radial resolution
        np(0:max_layers-1)   =  -1, &  ! poloidal resolution
        npL(0:max_layers-1)  =  -1, &  ! poloidal resolution in divertor legs
        npR(0:max_layers-1)  =  -1, &
        nt      = -1, &                ! number of cells in toroidal direction
        it_base = -1                   ! 0 <= index of base grid position <= nt

     real(real64) :: &
        tsize   = -360.d0              ! toroidal width of block [deg]
  end type tblock_parameters
  type(tblock_parameters) :: tblock(0:max_blocks-1)
  ! user defined parameters ....................................................



  ! toroidal discretization ....................................................
  type :: toroidal_discretization
     integer :: nt, it_base
     real(real64), pointer :: phi(:), phi_base => null()
  end type toroidal_discretization
  ! toroidal discretization ....................................................



  ! dependent parameters
  type(toroidal_discretization), allocatable :: T(:)
  real(real64) :: dphi
  logical :: default_layout

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_mmesh_parameters(mmesh_parameters)
  use moose_error
  use moose_math,  only: pi
  use flare_model, only: equi2d, assert_equi2d
  character(len=*), intent(in) :: mmesh_parameters

  character(len=256) :: iomsg, err
  integer :: iu, ios, ierr


  namelist /MmeshParameters/ &
     layout, pcoordinates, tblock, updown_symmetry, stellarator_symmetry, &
     radial_spacing, poloidal_spacing, poloidal_spacing_L, poloidal_spacing_R, &
     guiding_contour, guiding_contour_L, guiding_contour_R, &
     upstream_adjust_L, upstream_adjust_R, auto_adjust_mesh, &
     qmesh_generator, core_domain, vacuum_domain, plate_generator, closure_L, closure_R, &
     symmetry, blocks, nt, np, npL, npR, nr, n_interpolate, npXqo, npSP_subres, npSP_extend, &
     nr_core, dp_core, dt_core, nr_vac, dp_vac, dt_vac, delta_r, rinc, divmax, divavg, &
     cell_def, cell_param, nsside, plate_format, EIRENE_CORE_SF_NUM, EIRENE_VAC_SF_NUM, &
     phi0, p1, p2, d_SOL, d_PFR, polo_ext, &
     Bmod_in_vacuum_domain, generate_add_sf_n0



  ! 0. automatic range limits
  if (n_interpolate < 0) n_interpolate = 0


  ! 1. read user defined parameters
!  open  (newunit=iu, file="mmesh_parameters.txt", action='read')
!  read  (iu, MmeshParameters, iostat=ios, iomsg=iomsg)
!  close (iu)
  read  (mmesh_parameters, MmeshParameters, iostat=ios, iomsg=iomsg)
  if (ios /= 0) then
     print *, "while reading mmesh_parameters from input file:"
     call ERROR(trim(iomsg))
  endif
  EIRENE_SF_NUM(1) = EIRENE_CORE_SF_NUM
  EIRENE_SF_NUM(2) = EIRENE_VAC_SF_NUM
  EIRENE_SF_NUM(3) = EIRENE_VAC_SF_NUM ! grid closure at plates in vacuum domain
  if (updown_symmetry /= 0) stellarator_symmetry = .true.
  call assert_layout()


  ! 2. convert p1 and p2 to cylindrical coordinates
  select case(pcoordinates)
  case("", "cylindrical")

  case("magnetic")
     call assert_equi2d("load_mmesh_parameters")
     p1(1:2) = equi2d%rzcoords(p1(1), p1(2) / 180.d0 * pi, ierr=ierr)
     if (ierr /= 0) then
        print 9002, 'p1'
        stop
     endif

     p2(1:2) = equi2d%rzcoords(p2(1), p2(2) / 180.d0 * pi, ierr=ierr)
     if (ierr /= 0) then
        print 9002, 'p2'
        stop
     endif

  case default
     write (err, 9003) trim(pcoordinates)
     call ERROR(err)
  end select
 9002 format("failed to convert ",a," to cylindrical coordinates")
 9003 format("invalid choice '",a,"' for pcoordinates")


  ! 3. convert cm to m
  d_SOL = d_SOL / 1.d2
  d_PFR = d_PFR / 1.d2


  contains
  !.............................................................................
  subroutine assert_layout()
  use moose_error


  select case(layout)
  case(LAYOUT_GENERIC, "", "default")
     layout = LAYOUT_GENERIC

  case(LAYOUT_LSN, "lower_single_null")
     layout = LAYOUT_LSN

  case(LAYOUT_DDN, "disconnected_double_null")
     layout = LAYOUT_DDN

  case(LAYOUT_CDN, "connected_double_null")
     layout = LAYOUT_CDN

  case(LAYOUT_UNSTRUCTURED)
     layout = LAYOUT_UNSTRUCTURED

  case default
     call VALUE_ERROR("invalid layout '" // trim(layout) // "'")
  end select

  end subroutine assert_layout
  !.............................................................................
  end subroutine load_mmesh_parameters
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast_mmesh_parameters()
  use moose_mpi


  integer :: iblock

  call proc(0)%broadcast(layout)
  call proc(0)%broadcast(pcoordinates)
  call proc(0)%broadcast(updown_symmetry)
  call proc(0)%broadcast(stellarator_symmetry)
  call proc(0)%broadcast(radial_spacing)
  call proc(0)%broadcast(poloidal_spacing)
  call proc(0)%broadcast(poloidal_spacing_L)
  call proc(0)%broadcast(poloidal_spacing_R)
  call proc(0)%broadcast(guiding_contour)
  call proc(0)%broadcast(guiding_contour_L)
  call proc(0)%broadcast(guiding_contour_R)
  call proc(0)%broadcast(upstream_adjust_L)
  call proc(0)%broadcast(upstream_adjust_R)
  call proc(0)%broadcast(auto_adjust_mesh)
  call proc(0)%broadcast(qmesh_generator)
  call proc(0)%broadcast(core_domain)
  call proc(0)%broadcast(vacuum_domain)
  call proc(0)%broadcast(closure_L)
  call proc(0)%broadcast(closure_R)
  call proc(0)%broadcast(plate_generator)
  call proc(0)%broadcast(symmetry)
  call proc(0)%broadcast(blocks)
  call proc(0)%broadcast(nt)
  call proc(0)%broadcast(np)
  call proc(0)%broadcast(npL)
  call proc(0)%broadcast(npR)
  call proc(0)%broadcast(nr)
  call proc(0)%broadcast(n_interpolate)
  call proc(0)%broadcast(npXqo)
  call proc(0)%broadcast(npSP_subres)
  call proc(0)%broadcast(npSP_extend)
  call proc(0)%broadcast(nr_core)
  call proc(0)%broadcast(dp_core)
  call proc(0)%broadcast(dt_core)
  call proc(0)%broadcast(nr_vac)
  call proc(0)%broadcast(dp_vac)
  call proc(0)%broadcast(dt_vac)
  call proc(0)%broadcast(delta_r)
  call proc(0)%broadcast(rinc)
  call proc(0)%broadcast(divmax)
  call proc(0)%broadcast(divavg)
  call proc(0)%broadcast(cell_def)
  call proc(0)%broadcast(cell_param)
  call proc(0)%broadcast(nsside)
  call proc(0)%broadcast(plate_format)
  call proc(0)%broadcast(EIRENE_CORE_SF_NUM)
  call proc(0)%broadcast(EIRENE_VAC_SF_NUM)
  call proc(0)%broadcast(EIRENE_SF_NUM)
  call proc(0)%broadcast(phi0)
  call proc(0)%broadcast(p1)
  call proc(0)%broadcast(p2)
  call proc(0)%broadcast(d_SOL)
  call proc(0)%broadcast(d_PFR)
  call proc(0)%broadcast(polo_ext)
  call proc(0)%broadcast(Bmod_in_vacuum_domain)
  call proc(0)%broadcast(generate_add_sf_n0)
  do iblock=0,max_blocks-1
     call proc(0)%broadcast(tblock(iblock)%nr)
     call proc(0)%broadcast(tblock(iblock)%np)
     call proc(0)%broadcast(tblock(iblock)%npL)
     call proc(0)%broadcast(tblock(iblock)%npR)
     call proc(0)%broadcast(tblock(iblock)%nt)
     call proc(0)%broadcast(tblock(iblock)%it_base)
     call proc(0)%broadcast(tblock(iblock)%tsize)
  enddo

  end subroutine broadcast_mmesh_parameters
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine setup_mmesh_parameters()
  use moose_error
  use moose_math,  only: linspace
  use flare_control
  use flare_mmesh, only: sector, nsymmetry

  character(len=128) :: err
  real(real64) :: dphi_sum, phi1, phi2
  integer :: i, j


  ! 1. resolution and index of base mesh
  default_layout = .true.
  do i=0,blocks-1
     ! toroidal resolution
     call set_value(tblock(i)%nt, nt, "nt")

     ! radial and poloidal resolution
     do j=0,max_layers-1;   call set_value(tblock(i)%nr(j),  nr(j),  "nr");   enddo
     do j=0,max_layers-1;   call set_value(tblock(i)%np(j),  np(j),  "np");   enddo
     do j=0,max_layers-1;   call set_value(tblock(i)%npL(j), npL(j), "npL");   enddo
     do j=0,max_layers-1;   call set_value(tblock(i)%npR(j), npR(j), "npR");   enddo

     ! base mesh index
     if (tblock(i)%it_base == -1) then
        ! set default index
        tblock(i)%it_base = tblock(i)%nt / 2
        if (mod(tblock(i)%nt, 2) == 1) default_layout = .false.

        ! up-down symmetric slice at lower boundary
        if (updown_symmetry < 0  .and.  i == 0) then
           tblock(i)%it_base = 0
           default_layout = .false.
        ! up-down symmetric slice at upper boundary
        elseif (updown_symmetry > 0  .and.  i == blocks-1) then
           tblock(i)%it_base = tblock(i)%nt
           default_layout = .false.
        endif
     else
        default_layout = .false.
     endif

     ! check input
     if (tblock(i)%it_base < 0  .or.  tblock(i)%it_base > tblock(i)%nt) then
        write (err, 9001) tblock(i)%it_base, i
        call INDEX_ERROR(err)
     endif
  enddo
 9001 format("invalid index ",i0," for base mesh position in block ",i0)


  ! 2. toroidal layout of simulation domain
  dphi = sector(nsymmetry(symmetry, stellarator_symmetry))
  dphi_sum = 0.d0
  do i=0,blocks-1
     ! default size
     if (tblock(i)%tsize < 0.d0) then
        tblock(i)%tsize = dphi / blocks

     ! non-default size
     else
        default_layout = .false.
     endif
     dphi_sum = dphi_sum + tblock(i)%tsize
  enddo
  if (abs(dphi_sum - dphi) > 1.d-10) then
     call ERROR("block sizes are incompatible with total size")
  endif

  ! 2.3. lower toroidal boundary of simulation domain
  if (phi0 == -360.d0) then
     phi0 = -dphi / 2
     if (stellarator_symmetry) phi0 = 0.d0
  else
     default_layout = .false.
  endif


  ! 3. toroidal discretization
  allocate (T(0:blocks-1))
  phi1 = phi0
  do i=0,blocks-1
     T(i)%nt = tblock(i)%nt
     T(i)%it_base = tblock(i)%it_base

     phi2 = phi1 + tblock(i)%tsize
     allocate (T(i)%phi(0:T(i)%nt), source=linspace(phi1, phi2, T(i)%nt+1))
     T(i)%phi_base => T(i)%phi(T(i)%it_base)

     phi1 = phi2
  enddo


  ! SCREEN OUTPUT
  if (report) then
     print *
     if (layout == LAYOUT_UNSTRUCTURED) then
        print 1010, phi0, phi0 + dphi
        print 1011, T(0)%phi(T(0)%it_base)

     else
        print 1000, phi0, phi0 + dphi
        if (default_layout) print 1004
        print 1001
        do i=0,blocks-1
           print 1002, i, T(i)%phi(T(i)%it_base), T(i)%phi(0), T(i)%phi(T(i)%nt)
        enddo
     endif
     print *
  endif
 1000 format (3x,'- Layout of simulation domain (',f8.3,' -> ',f8.3,' deg):')
 1001 format (8x,'block #,  base mesh [deg],  domain [deg]')
 1002 format (8x,      i7,7x,f8.3,':',5x,f8.3,' -> ',f8.3)
 1004 format (8x,'using default layout')
 1010 format (3x,'- Simulation domain: ',f8.3,' -> ',f8.3,' deg')
 1011 format (8x,'base mesh at ',f8.3,' deg')

  contains
  !.............................................................................
  subroutine set_value(n, n_default, varname)
  integer, intent(inout) :: n
  integer,          intent(in   ) :: n_default
  character(len=*), intent(in   ) :: varname


  if (n == -1) n = n_default
  if (n <= 0) call VALUE_ERROR(varname//" > 0 required")

  end subroutine set_value
  !.............................................................................
  end subroutine setup_mmesh_parameters
  !-----------------------------------------------------------------------------

end module flare_mmesh_parameters
