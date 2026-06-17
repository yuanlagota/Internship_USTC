subroutine flare_task_fieldline_connection(grid, lcmax, lptmax, lttmax, output, min_psiN, &
  xbwd, xfwd, ubwd, ufwd, alpha, final_psiN, ierr)
  !
  ! Compute connection length and radial connection for selected field lines.
  !
  ! **Parameters:**
  !
  ! :grid:    Filename for grid with initial points for field line tracing.
  !
  ! :lcmax:   Max. length [m] for field line tracing.
  !
  ! :lptmax:  Max. length [poloidal turns] for field line tracing.
  !
  ! :lttmax:  Max. length [toroidal turns] for field line tracing.
  !
  ! :output:  Filename for output of data set.
  !
  ! :min_psiN:     Stop field line trace when :math:`\psi_N` < *min_psiN*.
  !
  ! :xbwd, xfwd:   Store (R[m], Z[m], phi[rad]) coordinates of intersection point with boundary in backward/forward direction.
  !
  ! :ubwd, ufwd:   Store id and surface coordinates (u1, u2) of intersection point with boundary in backward/forward direction.
  !
  ! :alpha:        Store grazing angle on target in shortest and longest direction.
  !
  ! :final_psiN:   Store psiN at end point of field line.
  !
  ! :ierr:         Store status code of field line trace
  !
  use iso_fortran_env
  use moose_mpi
  use moose_error
  use moose_math,  only: pi2, deg3
  use moose_dataset
  use moose_grids, only: r3grid
  use flare_control
  use flare_model, only: bfield, boundary
  use flare_fieldline
  use flare_tasks
  implicit none
  character(len=*), intent(in   ) :: grid, output
  real(real64),     intent(in   ) :: lcmax, lptmax, lttmax, min_psiN
  logical,          intent(in   ) :: xbwd, xfwd, ubwd, ufwd, alpha, final_psiN, ierr

  real(real64), pointer :: element(:)
  real(real64), allocatable :: y(:)
  type(r3grid)  :: G
  type(dataset) :: D
  type(fdriver) :: F
  real(real64)  :: x(3), lc(-1:1), lpt(-1:1), ltt(-1:1), minPsiN(-1:1), psiN(-1:1), xb(6,-1:1), a(-1:1)
  real(real64)  :: vb(3), vn(3), aopt(2)
  integer       :: i, idir, istat(-1:1), ixbwd, ixfwd, iubwd, iufwd, ipsiN, ialpha, iierr, jdir, n, m


  call begin_task()
  if (rank == 0) then
     print *, "Sampling connection length on grid ", trim(grid)
     print *

     G = r3grid(grid)
     print 1000, G%nnodes()
     print *
  endif
  call G%broadcast()
 1000 format(3x,"- Grid resolution: ",i0," nodes")


  ! initialize optional output & storage
  m = 8
  call init_optional_output(xbwd, 3, m ,ixbwd)
  call init_optional_output(xfwd, 3, m ,ixfwd)
  call init_optional_output(ubwd, 3, m ,iubwd)
  call init_optional_output(ufwd, 3, m ,iufwd)
  call init_optional_output(alpha, 2, m ,ialpha)
  call init_optional_output(final_psiN, 2, m ,ipsiN)
  call init_optional_output(ierr, 2, m ,iierr)
  n = G%nnodes()
  D = dataset(m, n)
  allocate (y(m))

  ! set output metadata
  call D%set_metadata(1, 'Lc_bwd',  "Backward connection length", "m")
  call D%set_metadata(2, 'Lc_fwd',  "Forward connection length",  "m")
  call D%set_metadata(3, 'Lpt_bwd', "Backward connection length", "{poloidal turns}")
  call D%set_metadata(4, 'Lpt_fwd', "Forward connection length",  "{poloidal turns}")
  call D%set_metadata(5, 'Ltt_bwd', "Backward connection length", "{toroidal turns}")
  call D%set_metadata(6, 'Ltt_fwd', "Forward connection length",  "{toroidal turns}")
  call D%set_metadata(7, 'minPsiN_bwd', "Minimum psiN in backward direction")
  call D%set_metadata(8, 'minPsiN_fwd', "Minimum psiN in forward direction")
  call D%set_expression("Lc",      "Lc_bwd + Lc_fwd", "Connection length", "m")
  call D%set_expression("Lcs",     "min(Lc_bwd, Lc_fwd)", "Shortest connection length", "m")
  call D%set_expression("Lpt",     "abs(Lpt_bwd - Lpt_fwd)", "Connection length", "{poloidal turns}")
  call D%set_expression("Lpts",    "min(abs(Lpt_bwd), abs(Lpt_fwd))", "Shortest connection length", "{poloidal turns}")
  call D%set_expression("Ltt",     "abs(Ltt_bwd - Ltt_fwd)", "Connection length", "{toroidal turns}")
  call D%set_expression("Ltts",    "min(abs(Ltt_bwd), abs(Ltt_fwd))", "Shortest connection length", "{toroidal turns}")
  call D%set_expression("minPsiN", "min(minPsiN_bwd, minPsiN_fwd)", "Deepest radial incursion", "psiN")
  if (xbwd) then
     call D%set_metadata(ixbwd+1, 'R_bwd',   "Backward strike point", "m")
     call D%set_metadata(ixbwd+2, 'Z_bwd',   "Backward strike point", "m")
     call D%set_metadata(ixbwd+3, 'phi_bwd', "Backward strike point", "{rad}")
  endif
  if (xfwd) then
     call D%set_metadata(ixfwd+1, 'R_fwd',   "Forward strike point", "m")
     call D%set_metadata(ixfwd+2, 'Z_fwd',   "Forward strike point", "m")
     call D%set_metadata(ixfwd+3, 'phi_fwd', "Forward strike point", "{rad}")
  endif
  if (ubwd) then
     call D%set_metadata(iubwd+1, 'n_bwd',   "Backward strike point boundary id")
     call D%set_metadata(iubwd+2, 'u1_bwd',  "Backward strike point")
     call D%set_metadata(iubwd+3, 'u2_bwd',  "Backward strike point")
  endif
  if (ufwd) then
     call D%set_metadata(iufwd+1, 'n_fwd',   "Forward strike point boundary id")
     call D%set_metadata(iufwd+2, 'u1_fwd',  "Forward strike point")
     call D%set_metadata(iufwd+3, 'u2_fwd',  "Forward strike point")
  endif
  if (alpha) then
     call D%set_metadata(ialpha+1, 'alphaS', "Incident angle [deg] in shortest direction")
     call D%set_metadata(ialpha+2, 'alphaL', "Incident angle [deg] in longest direction")
  endif
  if (final_psiN) then
     call D%set_metadata(ipsiN+1, 'psiN_bwd', "Final psiN (backward direction)")
     call D%set_metadata(ipsiN+2, 'psiN_fwd', "Final psiN (forward direction)")
     call D%set_expression("dpsiN", "(psiN_fwd - psiN_bwd)/2", "Delta psiN")
  endif
  if (ierr) then
     call D%set_metadata(iierr+1, 'ierr_bwd', "Status code (backward direction)")
     call D%set_metadata(iierr+2, 'ierr_fwd', "Status code (forward direction)")
  endif


  ! main loop
  call progress_bar(0, n)
  F = fdriver(max_arclength=lcmax, max_theta=pi2*lptmax, min_psiN=min_psiN)
  do i=rank,n-1,nproc
     a = 0.d0
     do idir=-1,1,2
        call F%reset()
        jdir          = bfield%equi%Bt_sign * idir

        x = G%node(i)
        ltt(jdir) = - x(3) / pi2
        istat(idir) = F%evolve3(x, x(3) + idir*pi2*lttmax)
        if (.not.FIELDLINE_SUCCESS(istat(idir)) .and. .not.ierr) then
           print *, "initial point: ", deg3(G%node(i))
           call FIELDLINE_ERROR(F, istat(idir))
        endif

        lc(jdir)      = F%arclength
        lpt(jdir)     = F%theta / pi2
        ltt(jdir)     = ltt(jdir) + x(3) / pi2
        minPsiN(jdir) = F%minPsiN
        psiN(jdir)    = F%psiN
        xb(:,jdir) = 0.d0
        if (istat(idir) == INTERSECT_BOUNDARY) then
           xb(:,jdir) = [x, 1.d0*F%nb, F%ub]
           vb = bfield%eval(x)
           vn = boundary%normal_vector(F%nb, F%ub)
           a(jdir) = asin(abs(sum(vb * vn)) / norm2(vb)) / pi2 * 360.d0
        endif
     enddo
     y(1:8) = [lc(-1), lc(1), lpt(-1), lpt(1), ltt(-1), ltt(1), minPsiN(-1), minPsiN(1)]
     if (lc(-1) < lc(1)) then
        aopt = [a(-1), a(1)]
     else
        aopt = [a(1), a(-1)]
     endif
     call store_optional_output(xbwd, y, ixbwd, xb(1:3,-1))
     call store_optional_output(xfwd, y, ixfwd, xb(1:3, 1))
     call store_optional_output(ubwd, y, iubwd, xb(4:6,-1))
     call store_optional_output(ufwd, y, iufwd, xb(4:6, 1))
     call store_optional_output(alpha, y, ialpha, aopt)
     call store_optional_output(final_psiN, y, ipsiN, [psiN(-1), psiN(1)])
     call store_optional_output(ierr, y, iierr, 1.d0*[istat(-1), istat(1)])
     element => D%element(i);   element = y

     if (verbose) print 2000, i, lc(-1)+lc(1)
     call progress_bar(i+1, n)
  enddo
  call finalize_progress_bar()
  call F%free()
  call D%allreduce()
 2000 format (5x,i8,',',8x,'L_c = ',f9.2,' m')


  ! save output
  if (rank == 0) then
     call D%set_geometry(grid, output)
     call D%savetxt(output)
  endif
  deallocate (y)
  call D%free()
  call finalize_task()

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine init_optional_output(flag, increment, count, ioffset)
  logical, intent(in   ) :: flag
  integer, intent(in   ) :: increment
  integer, intent(inout) :: count
  integer, intent(  out) :: ioffset


  if (flag) then
     ioffset = count
     count   = count + increment
  endif

  end subroutine init_optional_output
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine store_optional_output(flag, y, ioffset, values)
  logical,      intent(in   ) :: flag
  real(real64), intent(inout) :: y(:)
  integer,      intent(in   ) :: ioffset
  real(real64), intent(in   ) :: values(:)


  if (flag) then
     y(ioffset+1:ioffset+size(values)) = values
  endif

  end subroutine store_optional_output
  !-----------------------------------------------------------------------------

end subroutine flare_task_fieldline_connection
