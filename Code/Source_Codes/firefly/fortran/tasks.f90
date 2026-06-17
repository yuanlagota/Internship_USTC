module firefly_tasks
  use iso_fortran_env
  use moose_geometry, only: hypermesh3d
  use moose_data,     only: dataset
  implicit none


  ! results from strike point density and heat load proxy
  type(hypermesh3d), target :: pfc_mesh
  type(dataset) :: results

  real(real64), allocatable :: summary(:,:)
  logical :: results_initialized = .false.


  ! results from particle exhaust proxy: ionization in core and edge/SOL
  real(real64) :: sp(2), emin_fast, lost_particles

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_results(dphi, dl, kind)
  use moose_data,       only: CELL_DATA
  use firefly_geometry, only: pfc, nsurfaces
  use firefly_version
  real(real64),     intent(in) :: dphi, dl
  character(len=*), intent(in), optional :: kind


  if (results_initialized) then
     call pfc_mesh%free()
     call results%free()
     deallocate (summary)
  endif


  allocate (summary(2, nsurfaces), source = 0.d0)
  pfc_mesh = hypermesh3d(pfc, dphi, dl)
  results = dataset(2, pfc_mesh, CELL_DATA, kind=kind)
  results_initialized = .true.
  call results%annotations%set("firefly_version", version)

  end subroutine init_results
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine connection_length(filename, max_lc, output)
  !
  ! Compute field line connection length for set of initial points. This is the
  ! equivalent of the :ref:`fieldline_connection` task but using field line
  ! reconstruction instead of numerical integration.
  !
  ! **Parameters:**
  !
  ! :filename:  Name of :class:`moose:moose.grids.R3grid` file for initial points.
  !
  ! :max_lc:    Truncation for field line length [m].
  !
  use moose_mpi
  use moose_grids, only: r3grid
  use moose_dataset
  use flare_control, only: report, progress_bar, finalize_progress_bar
  use firefly_geometry
  character(len=*), intent(in) :: filename, output
  real(real64),     intent(in) :: max_lc

  type(r3grid)  :: grid
  type(dataset) :: D
  type(mcoords) :: c, c0
  real(real64), pointer :: element(:)
  real(real64) :: dl, lc(-1:1), t, t1, t2, tx, p1(3), p2(3), x(3), x0(3), u(2)
  integer :: i, idt, idt0, ierr, iphi, iphi1, iphi_tube(0:1), n, nnodes


  call assert_mmesh_and_pfc("connection_length")


  if (report) then
     print *, "Loading grid for sampling of connection length ..."
     grid = r3grid(filename)
     print *, "... done"
     print *

     print *, "Computing connection length ..."
     print *
  endif
  call grid%broadcast()


  call cpu_time(t1)
  nnodes = grid%nnodes()
  D = dataset(3, nnodes)
  call D%set_geometry(filename, output)
  call D%set_metadata(1, 'Lc_neg')
  call D%set_metadata(2, 'Lc_pos')
  call D%set_metadata(3, 'ierr')
  call D%set_expression("Lc", "Lc_neg + Lc_pos", label="Connection length", units="m")
  call progress_bar(0, nnodes)
  do i=rank,nnodes-1,nproc
     x0 = grid%node(i)
     call mesh%mcoords(x0, c0, 0, ierr)
     if (ierr /= 0) then
        !print *, x0, i, ierr
        D%values(3,i) = 1
        cycle
     endif

     do idt0=-1,1,2
        idt = idt0
        c = c0
        p1(1:2) = mesh%rzcoords(c)
        p1(3) = mesh%phi(c%iphi) + c%t * (mesh%phi(c%iphi+1) - mesh%phi(c%iphi))
        lc(idt) = tracing%arclength(c%itube, c%iphi)
        if (idt == 1) then
           iphi1 = c%iphi + 1
           t = 1.d0 - c%t
        else
           iphi1 = c%iphi
           t = c%t
        endif
        lc(idt) = lc(idt) * t

        p2(1:2) = mesh%rzcoords(c%itube, iphi1, c%xi)
        p2(3) = mesh%phi(iphi1)
        if (tracing%boundary_flag(c%itube) == 1 .and. pfc%intersect(p1, p2, x, tx, n, u)) then
           lc(idt) = lc(idt) * tx
           cycle
        endif


        trace_loop: do
           iphi_tube = mesh%iphi_zone(:, mesh%izone_tube(c%itube))
           p1(1:2) = mesh%rzcoords(c%itube, iphi1, c%xi)
           p1(3) = mesh%phi(iphi1)
           do iphi=iphi1,iphi_tube(torosf_side(idt))-idt,idt
              dl = tracing%arclength(c%itube, c%iphi)
              if (tracing%boundary_flag(c%itube) == 1) then
              p2(1:2) = mesh%rzcoords(c%itube, iphi+idt, c%xi)
              p2(3) = mesh%phi(iphi+idt)

              if (pfc%intersect(p1, p2, x, tx, n, u)) then
                 lc(idt0) = lc(idt0) + tx * dl
                 exit trace_loop
              endif
              p1 = p2
              endif
              lc(idt0) = lc(idt0) + dl
           enddo
           if (lc(idt0) > max_lc) exit trace_loop

           call mesh%torosf_map(c, idt, ierr)
           if (ierr /= 0) then
              !print *, i, c, ierr
              D%values(3,i) = 2
              exit trace_loop
           endif
           iphi1 = mesh%iphi_zone(torosf_side(-idt), mesh%izone_tube(c%itube))
        enddo trace_loop
     enddo

     element => D%element(i)
     element(1) = lc(-1)
     element(2) = lc(1)
     call progress_bar(i+1, nnodes)
  enddo
  call finalize_progress_bar()
  call D%allreduce()
  call cpu_time(t2)
  if (rank == 0) print *, "... finished, total time: ", t2 - t1, " s"

  if (rank == 0) call D%savetxt(output)
  call D%free()

  end subroutine connection_length
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine strike_point_density(dcoeff, nsamples, bstep, dphi, dl, output)
  !
  ! Compute the strike point density of magnetic field lines on divertor targets,
  ! baffles and the first wall. This is the equivalent of the :ref:`fieldline_connection`
  ! task with artifical cross-field diffusion and processing the output with the
  ! :ref:`strike_point_density` task. Field lines are initialized along the
  ! inner mesh boundary here. The output is a :class:`moose:moose.data.Dataset`
  ! that can be visualized with the :ref:`moose:mplot` and/or :ref:`moose:mplot3d` commands.
  !
  !
  ! **Parameters:**
  !
  ! :dcoeff:     Field line diffusion coefficient [m**2 / m].
  !
  ! :nsamples:   Number of field line samples.
  !
  ! :bstep:      Step size along field lines [m].
  !
  ! :dphi:       Toroidal resolution for output mesh [deg].
  !
  ! :dl:         Output mesh resolution along pfc [m].
  !
  use moose_mpi
  use moose_math
  use moose_geometry, only: hypermesh3d
  use moose_dataset
  use flare_control
  use flare_mmesh_unstructured_sampling
  use flare_mmesh_unstructured_tracing
  use firefly_geometry
  real(real64),     intent(in) :: dcoeff, bstep, dphi, dl
  integer,          intent(in) :: nsamples
  character(len=*), intent(in) :: output

  type(mcoords_workspace) :: c
  type(source) :: S
  real(real64) :: w, xstep, dx(2), ds, u(2), loss(4), r, total_src, t1, t2
  integer :: i, icell, ids, istat, k1, k2, n


  call assert_mmesh_and_pfc("strike_point_density")


  total_src = 1.d0 / abs(mesh%symmetry);   if (mesh%symmetry < 0) total_src = total_src / 2
  if (report) then
     print *, "Computing strike point density ..."
     print *
     print *, "dcoeff = ", dcoeff
     print *, "nsamples = ", nsamples
     print *, "dphi, dl = ", dphi, dl
     print *, "total source = ", total_src
     print *
  endif


  ! initialize pfc mesh with fine resolution
  call init_results(dphi, dl, "strike_point_density")
  associate (M => pfc_mesh, D => results)
  call cpu_time(t1)
  call D%set_metadata(1, "p", "Strike point density", "m**(-2)")
  call D%set_metadata(2, "area", "Surface area per cell", "m**2")
  call D%set_parameter("dcoeff", dcoeff)
  call D%set_parameter("bstep", bstep)
  call D%set_parameter("nsamples", nsamples)


  ! jump steps [m]
  xstep = sqrt(2 * bstep * dcoeff)
  if (report) then
     print *, "xstep, bstep = ", xstep, bstep
     print *
  endif


  loss = 0.d0
  S = surface_source(mesh, ISB_TAG, 1, total_src, nsamples)
  call progress_bar(0, nsamples)
  do i=1+rank,nsamples,nproc
     ! sample new particle
     call S%sample(c%mcoords, w)
     if (tracing%boundary_flag(c%itube) == 1) c%p = mesh%rzphicoords(c%mcoords)
     call random_number(r)
     ids = 2 * int(2*r) - 1

     ! trace particle until it reaches divertor targets
     trace_loop: do
        ! 1. cross-field step
        dx = random_number_stdnorm2d() * xstep
        call tracing%xstep(c, dx, istat, n, u)
        ! - at targets
        if (istat == -1) then
           icell = M%cell_index(n, u)
           D%values(1,icell) = D%values(1,icell) + w
           exit trace_loop

        ! - reflection at inner boundary
        elseif (mod(istat, 10) == ISB_TAG) then
           c%xi(2) = -0.9999876543210d0
           istat = 0

        ! - out of bounds
        elseif (istat /= 0) then
           !print *, "cross-field step exited with istat = ", istat
           loss(1) = loss(1) + w
           exit trace_loop
        endif


        ! 2. step along field line
        ds = ids * bstep
        bstep_loop: do
           call tracing%bstep(c, ds, istat, n, u)
           ! - finished bstep
           if (istat == 0) then
               ids = ids * n
               exit bstep_loop

           ! - at targets
           elseif (istat == -1) then
              icell = M%cell_index(n, u)
              D%values(1,icell) = D%values(1,icell) + w
              exit trace_loop

           ! - out of bounds
           else
              loss(2) = loss(2) + w
              exit trace_loop
           endif
        enddo bstep_loop
     enddo trace_loop
     call progress_bar(i, nsamples)
  enddo
  call finalize_progress_bar()
  call D%allreduce()
  call cpu_time(t2)


  ! prepare summary and save results
  ! 1. collect integral fluxes
  do i=1,size(M%refined_tpzmesh3d)
     k1 = M%cell_offset(i)
     k2 = M%cell_offset(i+1) - 1
     summary(1, i) = sum(D%values(1,k1:k2)) / total_src * 100.d0
  enddo

  ! 2. evaluate strike point density
  D%values(2,:) = M%area()
  where (D%values(2,:) > 0.d0) D%values(1,:) = D%values(1,:) / D%values(2,:)
  do i=1,size(M%refined_tpzmesh3d)
     k1 = M%cell_offset(i)
     k2 = M%cell_offset(i+1) - 1
     summary(2, i) = maxval(D%values(1,k1:k2))
  enddo
  if (rank == 0) call D%savenc(output)


  ! print summary
  if (rank == 0) then
     print *
     print *, "results:                            contribution             max. value [m**(-2)]"
     do i=1,size(M%refined_tpzmesh3d)
        print 4000, M%refined_tpzmesh3d(i)%title(), summary(:, i)
     enddo
     print *
     if (any(loss > 0.d0)) print *, "loss = ", loss
     print *, "computation time: ", t2 - t1, " s"
  endif
 4000 format(1x,a32,2x,f8.3," %",8x,e12.5)

  end associate
  end subroutine strike_point_density
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine heat_load_proxy(n0, T0, chi, nparticles, tau, dphi, dl, output)
  !
  ! Compute heat loads onto divertor targets, baffles and the first wall based
  ! on fixed background plasma (linearized heat conduction). Heat loads are
  ! normalized to the total power exhaust.
  ! The output is a :class:`moose:moose.data.Dataset` that can be visualized
  ! with the :ref:`moose:mplot` and/or :ref:`moose:mplot3d` commands.
  !
  !
  ! **Parameters:**
  !
  ! :n0:         Background plasma density [m**(-3)].
  !
  ! :T0:         Background plasma temperature [eV].
  !
  ! :chi:        Cross-field heat diffusion coefficient [m**2 / m].
  !
  ! :nparticles: Number of test particles.
  !
  ! :tau:        Time step for particle tracing [s].
  !
  ! :dphi:       Toroidal resolution for output mesh [deg].
  !
  ! :dl:         Output mesh resolution along pfc [m].
  !
  use moose_mpi
  use moose_math
  use moose_geometry, only: hypermesh3d
  use moose_dataset
  use flare_control
  use flare_mmesh_unstructured_sampling
  use flare_mmesh_unstructured_tracing
  use firefly_geometry
  real(real64),     intent(in) :: n0, T0, chi
  integer,          intent(in) :: nparticles
  real(real64),     intent(in) :: tau, dphi, dl
  character(len=*), intent(in) :: output

  real(real64), parameter :: &
     echarg        = 1.602176634d-19, &
     mp            = 1.67262192369d-27, &
     cs_prefac     = sqrt(echarg / mp), &
     kappa_prefac  = 2060 / echarg, &
     sheath_factor = 7.d0

  type(mcoords_workspace) :: c
  type(source) :: S
  real(real64) :: cs, kappa, lt_para, a, pload
  real(real64) :: w, wload, xstep, bstep, dx(2), ds, u(2), loss(4), total_src, t1, t2
  integer :: i, icell, istat, k1, k2, n


  call assert_mmesh_and_pfc("head_load_proxy")


  total_src = 1.d0 / abs(mesh%symmetry);   if (mesh%symmetry < 0) total_src = total_src / 2
  if (report) then
     print *, "Computing heat load proxy ..."
     print *
     print *, "n0 = ", n0
     print *, "T0 = ", T0
     print *, "chi = ", chi
     print *, "nparticles = ", nparticles
     print *, "total source = ", total_src
     print *
  endif


  ! initialize pfc mesh with fine resolution
  call init_results(dphi, dl, "heat_load_proxy")
  associate (M => pfc_mesh, D => results)
  call cpu_time(t1)
  call D%set_parameter("n0", n0)
  call D%set_parameter("T0", T0)
  call D%set_parameter("chi", chi)
  call D%set_parameter("nparticles", nparticles)
  call D%set_parameter("tau", tau)
  call D%set_metadata(1, "hload", "Heat load / power", "m**(-2)")
  call D%set_metadata(2, "area", "Surface area per cell", "m**2")


  ! pre-compute transport coefficients
  ! sound speed of hydrogen [m/s] with Te = Ti
  cs = cs_prefac * sqrt(2 * T0)
  ! parallel heat conductivity (normalized to n0) [m**2 / s]
  kappa = kappa_prefac * T0**2.5d0 / n0
  ! parallel T-decay length [m]
  lt_para = kappa / sheath_factor / cs
  ! jump steps [m]
  xstep = sqrt(2 * tau * chi)
  bstep = sqrt(2 * tau * kappa)
  ! loss probability (target loading)
  a = sqrt(pi/2) * bstep / lt_para
  pload = a / (1.d0 + 0.5d0 * a)
  if (report) then
     print *, "xstep, bstep, pload = ", xstep, bstep, pload
     print *
  endif


  loss = 0.d0
  S = surface_source(mesh, ISB_TAG, 1, total_src, nparticles)
  call progress_bar(0, nparticles)
  do i=1+rank,nparticles,nproc
     ! sample new particle
     call S%sample(c%mcoords, w)
     if (tracing%boundary_flag(c%itube) == 1) c%p = mesh%rzphicoords(c%mcoords)

     ! trace particle until it reaches divertor targets
     trace_loop: do
        ! 1. cross-field step
        dx = random_number_stdnorm2d() * xstep
        call tracing%xstep(c, dx, istat, n, u)
        ! - at targets
        if (istat == -1) then
           icell = M%cell_index(n, u)
           D%values(1,icell) = D%values(1,icell) + w
           exit trace_loop

        ! - reflection at inner boundary
        elseif (mod(istat, 10) == ISB_TAG) then
           c%xi(2) = -0.9999876543210d0
!           dx = -dx
!           call tracing%xstep(c, dx, istat, n, u)
!           if (istat /= 0) then
!              loss(3) = loss(3) + w
!              exit trace_loop
!           endif
           istat = 0

        ! - out of bounds
        elseif (istat /= 0) then
           !print *, "cross-field step exited with istat = ", istat
           loss(1) = loss(1) + w
           exit trace_loop
        endif


        ! 2. step along field line
        ds = random_number_stdnorm() * bstep
        bstep_loop: do
           call tracing%bstep(c, ds, istat, n, u)
           ! - finished bstep
           if (istat == 0) then
               exit bstep_loop

           ! - at targets
           elseif (istat == -1) then
              call S%boundary_event(w, pload, wload)
              if (wload > 0.d0) then
                 icell = M%cell_index(n, u)
                 D%values(1,icell) = D%values(1,icell) + wload
              endif
              if (w == 0.d0) exit trace_loop
              ds = -ds

           ! - out of bounds
           else
              loss(2) = loss(2) + w
              exit trace_loop
           endif
        enddo bstep_loop
     enddo trace_loop
     call progress_bar(i, nparticles)
  enddo
  call finalize_progress_bar()
  call D%allreduce()
  call cpu_time(t2)


  ! prepare summary and save results
  ! 1. collect integral fluxes
  do i=1,size(M%refined_tpzmesh3d)
     k1 = M%cell_offset(i)
     k2 = M%cell_offset(i+1) - 1
     summary(1, i) = sum(D%values(1,k1:k2)) / total_src * 100.d0
  enddo

  ! 2. evaluate heat loads
  D%values(2,:) = M%area()
  where (D%values(2,:) > 0.d0) D%values(1,:) = D%values(1,:) / D%values(2,:)
  do i=1,size(M%refined_tpzmesh3d)
     k1 = M%cell_offset(i)
     k2 = M%cell_offset(i+1) - 1
     summary(2, i) = maxval(D%values(1,k1:k2))
  enddo
  if (rank == 0) call D%savenc(output)


  ! print summary
  if (report) then
     print *
     print *, "results:                            contribution             max. value [m**(-2)]"
     do i=1,size(M%refined_tpzmesh3d)
        print 4000, M%refined_tpzmesh3d(i)%title(), summary(:, i)
     enddo
     print *
     if (any(loss > 0.d0)) print *, "loss = ", loss
     print *, "computation time: ", t2 - t1, " s"
  endif
 4000 format(1x,a32,2x,f8.3," %",8x,e12.5)

  end associate
  end subroutine heat_load_proxy
  !-----------------------------------------------------------------------------

end module firefly_tasks
