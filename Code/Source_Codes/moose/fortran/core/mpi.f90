!===============================================================================
! Supplemental subroutines for MPI programs
!===============================================================================
module moose_mpi
  use mpi_f08
  use moose_kinds
  implicit none
  private

  integer, parameter, private :: chunk_size = 2**16


  integer, public, save :: nproc = 1, rank = 0


  ! reference to process
  type, public :: process
     integer :: id ! process id

     contains
     ! send data to this process
     generic   :: send => send_int, send_real, send_real1, send_real2
     procedure :: send_int
     procedure :: send_real
     procedure :: send_real1
     procedure :: send_real2

     ! receive data from this process
     generic   :: recv => recv_int, recv_real, recv_real1, recv_real2
     procedure :: recv_int
     procedure :: recv_real
     procedure :: recv_real1
     procedure :: recv_real2

     ! allocate on rank /= id and broadcast data
     generic   :: broadcast_allocatable => allocate_int1
     generic   :: broadcast_allocatable => allocate_int2
     generic   :: broadcast_allocatable => allocate_int3
     generic   :: broadcast_allocatable => allocate_int4
     generic   :: broadcast_allocatable => allocate_real1
     generic   :: broadcast_allocatable => allocate_real2
     generic   :: broadcast_allocatable => allocate_real3
     generic   :: broadcast_allocatable => allocate_real4
     generic   :: broadcast_allocatable => allocate_complex1
     generic   :: broadcast_allocatable => allocate_complex2
     generic   :: broadcast_allocatable => allocate_complex3
     generic   :: broadcast_allocatable => allocate_complex4
     generic   :: broadcast_allocatable => allocate_logical1
     generic   :: broadcast_allocatable => allocate_logical2
     generic   :: broadcast_allocatable => allocate_logical3
     generic   :: broadcast_allocatable => allocate_logical4
     generic   :: broadcast_allocatable => allocate_string0
     generic   :: broadcast_allocatable => allocate_string1
     procedure :: allocate_int1
     procedure :: allocate_int2
     procedure :: allocate_int3
     procedure :: allocate_int4
     procedure :: allocate_real1
     procedure :: allocate_real2
     procedure :: allocate_real3
     procedure :: allocate_real4
     procedure :: allocate_complex1
     procedure :: allocate_complex2
     procedure :: allocate_complex3
     procedure :: allocate_complex4
     procedure :: allocate_logical1
     procedure :: allocate_logical2
     procedure :: allocate_logical3
     procedure :: allocate_logical4
     procedure :: allocate_string0
     procedure :: allocate_string1


     ! broadcast data from this process
     generic   :: broadcast => broadcast_integer
     generic   :: broadcast => broadcast_int1
     generic   :: broadcast => broadcast_int2
     generic   :: broadcast => broadcast_int3
     generic   :: broadcast => broadcast_int4
     generic   :: broadcast => broadcast_real
     generic   :: broadcast => broadcast_real1
     generic   :: broadcast => broadcast_real2
     generic   :: broadcast => broadcast_real3
     generic   :: broadcast => broadcast_real4
     generic   :: broadcast => broadcast_complex
     generic   :: broadcast => broadcast_complex1
     generic   :: broadcast => broadcast_complex2
     generic   :: broadcast => broadcast_complex3
     generic   :: broadcast => broadcast_complex4
     generic   :: broadcast => broadcast_logical
     generic   :: broadcast => broadcast_logical1
     generic   :: broadcast => broadcast_logical2
     generic   :: broadcast => broadcast_logical3
     generic   :: broadcast => broadcast_logical4
     generic   :: broadcast => broadcast_string
     generic   :: broadcast => broadcast_string1
     procedure :: broadcast_integer
     procedure :: broadcast_int1
     procedure :: broadcast_int2
     procedure :: broadcast_int3
     procedure :: broadcast_int4
     procedure :: broadcast_real
     procedure :: broadcast_real1
     procedure :: broadcast_real2
     procedure :: broadcast_real3
     procedure :: broadcast_real4
     procedure :: broadcast_complex
     procedure :: broadcast_complex1
     procedure :: broadcast_complex2
     procedure :: broadcast_complex3
     procedure :: broadcast_complex4
     procedure :: broadcast_logical
     procedure :: broadcast_logical1
     procedure :: broadcast_logical2
     procedure :: broadcast_logical3
     procedure :: broadcast_logical4
     procedure :: broadcast_string
     procedure :: broadcast_string1
  end type process
  type(process), allocatable, public :: proc(:)



  ! sum data
  interface moose_mpi_sum
     procedure sum_int
     procedure sum_integer_array1
     procedure sum_integer_array2
     procedure sum_integer_array3
     procedure sum_integer_array4
     procedure sum_real
     procedure sum_real_array1
     procedure sum_real_array2
     procedure sum_real_array3
     procedure sum_real_array4
  end interface



  public :: &
     moose_mpi_init, &
     moose_mpi_setup, &
     moose_mpi_finalize, &
     moose_mpi_sum, &
     moose_mpi_range, &
     mpi_barrier_world


  contains
  !---------------------------------------------------------------------


! type process =================================================================
! send -------------------------------------------------------------------------
  subroutine send_int(this, I)
  class(process), intent(in) :: this
  integer,        intent(in) :: I

  integer :: ierr


  call MPI_SEND(I, 1, MPI_INTEGER, this%id, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "send_int")

  end subroutine send_int
  !-----------------------------------------------------------------------------
  subroutine send_real(this, R)
  class(process), intent(in) :: this
  real(real64),   intent(in) :: R

  integer :: ierr


  call MPI_SEND(R, 1, MPI_REAL8, this%id, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "send_real")

  end subroutine send_real
  !-----------------------------------------------------------------------------
  subroutine send_real1(this, R)
  class(process), intent(in) :: this
  real(real64),   intent(in) :: R(:)

  integer :: ierr


  call MPI_SEND(R, size(R), MPI_REAL8, this%id, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "send_real1")

  end subroutine send_real1
  !-----------------------------------------------------------------------------
  subroutine send_real2(this, R)
  class(process), intent(in) :: this
  real(real64),   intent(in) :: R(:,:)

  integer :: ierr


  call MPI_SEND(R, size(R), MPI_REAL8, this%id, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "send_real2")

  end subroutine send_real2
! send -------------------------------------------------------------------------


! recv --------------------------------------------------------------------------
  subroutine recv_int(this, I)
  class(process), intent(in)  :: this
  integer,        intent(out) :: I

  type(MPI_Status) :: status
  integer :: ierr


  call MPI_RECV(I, 1, MPI_INTEGER, this%id, 0, MPI_COMM_WORLD, status, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "recv_int")

  end subroutine recv_int
  !-----------------------------------------------------------------------------
  subroutine recv_real(this, R)
  class(process), intent(in)  :: this
  real(real64),   intent(out) :: R

  type(MPI_Status) :: status
  integer :: ierr


  call MPI_RECV(R, 1, MPI_REAL8, this%id, 0, MPI_COMM_WORLD, status, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "recv_real")

  end subroutine recv_real
  !-----------------------------------------------------------------------------
  subroutine recv_real1(this, R)
  class(process), intent(in)  :: this
  real(real64),   intent(out) :: R(:)

  type(MPI_Status) :: status
  integer :: ierr


  call MPI_RECV(R, size(R), MPI_REAL8, this%id, 0, MPI_COMM_WORLD, status, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "recv_real2")

  end subroutine recv_real1
  !-----------------------------------------------------------------------------
  subroutine recv_real2(this, R)
  class(process), intent(in)  :: this
  real(real64),   intent(out) :: R(:,:)

  type(MPI_Status) :: status
  integer :: ierr


  call MPI_RECV(R, size(R), MPI_REAL8, this%id, 0, MPI_COMM_WORLD, status, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "recv_real2")

  end subroutine recv_real2
! recv -------------------------------------------------------------------------


! allocate_and_broadcast -------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine allocate_int1(this, I)
  class(process),       intent(in)    :: this
  integer, allocatable, intent(inout) :: I(:)

  integer :: n, lb, ub


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(I)
     lb = lbound(I,1)
     ub = ubound(I,1)
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(I(lb:ub))

  ! broadcast data
  call this%broadcast(I)

  end subroutine allocate_int1
  !---------------------------------------------------------------------
  subroutine allocate_int2(this, I)
  class(process),       intent(in)    :: this
  integer, allocatable, intent(inout) :: I(:,:)

  integer :: k, n, lb(2), ub(2)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(I)
     do k=1,2
        lb(k) = lbound(I,k)
        ub(k) = ubound(I,k)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(I(lb(1):ub(1), lb(2):ub(2)))

  ! broadcast data
  call this%broadcast(I)

  end subroutine allocate_int2
  !---------------------------------------------------------------------
  subroutine allocate_int3(this, I)
  class(process),       intent(in)    :: this
  integer, allocatable, intent(inout) :: I(:,:,:)

  integer :: k, n, lb(3), ub(3)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(I)
     do k=1,3
        lb(k) = lbound(I,k)
        ub(k) = ubound(I,k)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(I(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))

  ! broadcast data
  call this%broadcast(I)

  end subroutine allocate_int3
  !---------------------------------------------------------------------
  subroutine allocate_int4(this, I)
  class(process),       intent(in)    :: this
  integer, allocatable, intent(inout) :: I(:,:,:,:)

  integer :: k, n, lb(4), ub(4)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(I)
     do k=1,4
        lb(k) = lbound(I,k)
        ub(k) = ubound(I,k)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(I(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))

  ! broadcast data
  call this%broadcast(I)

  end subroutine allocate_int4
  !---------------------------------------------------------------------
  subroutine allocate_real1(this, R)
  class(process),            intent(in)    :: this
  real(real64), allocatable, intent(inout) :: R(:)

  integer :: n, lb, ub


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(R)
     lb = lbound(R,1)
     ub = ubound(R,1)
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(R(lb:ub))

  ! broadcast data
  call this%broadcast(R)

  end subroutine allocate_real1
  !---------------------------------------------------------------------
  subroutine allocate_real2(this, R)
  class(process),            intent(in)    :: this
  real(real64), allocatable, intent(inout) :: R(:,:)

  integer :: i, n, lb(2), ub(2)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(R)
     do i=1,2
        lb(i) = lbound(R,i)
        ub(i) = ubound(R,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(R(lb(1):ub(1), lb(2):ub(2)))

  ! broadcast data
  call this%broadcast(R)

  end subroutine allocate_real2
  !---------------------------------------------------------------------
  subroutine allocate_real3(this, R)
  class(process),            intent(in)    :: this
  real(real64), allocatable, intent(inout) :: R(:,:,:)

  integer :: i, n, lb(3), ub(3)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(R)
     do i=1,3
        lb(i) = lbound(R,i)
        ub(i) = ubound(R,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(R(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))

  ! broadcast data
  call this%broadcast(R)

  end subroutine allocate_real3
  !---------------------------------------------------------------------
  subroutine allocate_real4(this, R)
  class(process),            intent(in)    :: this
  real(real64), allocatable, intent(inout) :: R(:,:,:,:)

  integer :: i, n, lb(4), ub(4)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(R)
     do i=1,4
        lb(i) = lbound(R,i)
        ub(i) = ubound(R,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(R(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))

  ! broadcast data
  call this%broadcast(R)

  end subroutine allocate_real4
  !---------------------------------------------------------------------
  subroutine allocate_complex1(this, C)
  class(process),                intent(in)    :: this
  complex(kind=dp), allocatable, intent(inout) :: C(:)

  integer :: n, lb, ub


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(C)
     lb = lbound(C,1)
     ub = ubound(C,1)
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(C(lb:ub))

  ! broadcast data
  call this%broadcast(C)

  end subroutine allocate_complex1
  !---------------------------------------------------------------------
  subroutine allocate_complex2(this, C)
  class(process),                intent(in)    :: this
  complex(kind=dp), allocatable, intent(inout) :: C(:,:)

  integer :: i, n, lb(2), ub(2)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(C)
     do i=1,2
        lb(i) = lbound(C,i)
        ub(i) = ubound(C,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(C(lb(1):ub(1), lb(2):ub(2)))

  ! broadcast data
  call this%broadcast(C)

  end subroutine allocate_complex2
  !---------------------------------------------------------------------
  subroutine allocate_complex3(this, C)
  class(process),                intent(in)    :: this
  complex(kind=dp), allocatable, intent(inout) :: C(:,:,:)

  integer :: i, n, lb(3), ub(3)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(C)
     do i=1,3
        lb(i) = lbound(C,i)
        ub(i) = ubound(C,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(C(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))

  ! broadcast data
  call this%broadcast(C)

  end subroutine allocate_complex3
  !---------------------------------------------------------------------
  subroutine allocate_complex4(this, C)
  class(process),                intent(in)    :: this
  complex(kind=dp), allocatable, intent(inout) :: C(:,:,:,:)

  integer :: i, n, lb(4), ub(4)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(C)
     do i=1,4
        lb(i) = lbound(C,i)
        ub(i) = ubound(C,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memoray
  if (rank /= this%id) allocate(C(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))

  ! broadcast data
  call this%broadcast(C)

  end subroutine allocate_complex4
  !---------------------------------------------------------------------
  subroutine allocate_logical1(this, L)
  class(process),       intent(in)    :: this
  logical, allocatable, intent(inout) :: L(:)

  integer :: n, lb, ub


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(L)
     lb = lbound(L,1)
     ub = ubound(L,1)
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(L(lb:ub))

  ! broadcast data
  call this%broadcast(L)

  end subroutine allocate_logical1
  !---------------------------------------------------------------------
  subroutine allocate_logical2(this, L)
  class(process),       intent(in)    :: this
  logical, allocatable, intent(inout) :: L(:,:)

  integer :: i, n, lb(2), ub(2)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(L)
     do i=1,2
        lb(i) = lbound(L,i)
        ub(i) = ubound(L,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(L(lb(1):ub(1), lb(2):ub(2)))

  ! broadcast data
  call this%broadcast(L)

  end subroutine allocate_logical2
  !---------------------------------------------------------------------
  subroutine allocate_logical3(this, L)
  class(process),       intent(in)    :: this
  logical, allocatable, intent(inout) :: L(:,:,:)

  integer :: i, n, lb(3), ub(3)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(L)
     do i=1,3
        lb(i) = lbound(L,i)
        ub(i) = ubound(L,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(L(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))

  ! broadcast data
  call this%broadcast(L)

  end subroutine allocate_logical3
  !---------------------------------------------------------------------
  subroutine allocate_logical4(this, L)
  class(process),       intent(in)    :: this
  logical, allocatable, intent(inout) :: L(:,:,:,:)

  integer :: i, n, lb(4), ub(4)


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(L)
     do i=1,4
        lb(i) = lbound(L,i)
        ub(i) = ubound(L,i)
     enddo
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(L(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4)))

  ! broadcast data
  call this%broadcast(L)

  end subroutine allocate_logical4
  !---------------------------------------------------------------------
  subroutine allocate_string0(this, S)
  class(process),            intent(in)    :: this
  character(:), allocatable, intent(inout) :: S

  integer :: n


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = len(S)
  endif
  call this%broadcast(n)

  ! allocate memory
  if (rank /= this%id) allocate(character(len=n) :: S)

  ! broadcast data
  call this%broadcast(S)

  end subroutine allocate_string0
  !---------------------------------------------------------------------
  subroutine allocate_string1(this, S)
  class(process),                intent(in)    :: this
  character(len=*), allocatable, intent(inout) :: S(:)

  integer :: n, lb, ub


  ! array size and boundaries (this process)
  if (rank == this%id) then
     n  = size(S)
     lb = lbound(S,1)
     ub = ubound(S,1)
  endif
  call this%broadcast(n)
  call this%broadcast(lb)
  call this%broadcast(ub)

  ! allocate memory
  if (rank /= this%id) allocate(S(lb:ub))

  ! broadcast data
  call this%broadcast(S)

  end subroutine allocate_string1
  !---------------------------------------------------------------------
! allocate_and_broadcast -------------------------------------------------------


! broadcast --------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine broadcast_integer(this, I)
  class(process), intent(in)    :: this
  integer,        intent(inout) :: I

  integer :: ierr


  call MPI_BCAST(I, 1, MPI_INTEGER, this%id, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "broadcast_integer")

  end subroutine broadcast_integer
  !-----------------------------------------------------------------------------
  subroutine broadcast_real(this, R)
  class(process), intent(in)    :: this
  real(real64),   intent(inout) :: R

  integer :: ierr


  call MPI_BCAST(R, 1, MPI_REAL8, this%id, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "broadcast_real")

  end subroutine broadcast_real
  !-----------------------------------------------------------------------------
  subroutine broadcast_complex(this, C)
  class(process),   intent(in)    :: this
  complex(kind=dp), intent(inout) :: C

  integer :: ierr


  call MPI_BCAST(C, 1, MPI_DOUBLE_COMPLEX, this%id, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "broadcast_complex")

  end subroutine broadcast_complex
  !-----------------------------------------------------------------------------
  subroutine broadcast_logical(this, L)
  class(process),   intent(in)    :: this
  logical,          intent(inout) :: L

  integer :: ierr


  call MPI_BCAST(L, 1, MPI_LOGICAL, this%id, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "broadcast_logical")

  end subroutine broadcast_logical
  !-----------------------------------------------------------------------------
  subroutine broadcast_string(this, S)
  class(process),   intent(in)    :: this
  character(len=*), intent(inout) :: S

  integer :: ierr, n


  n = len(S)
  call MPI_BCAST(S, n, MPI_CHARACTER, this%id, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "broadcast_string")

  end subroutine broadcast_string
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine aux_broadcast_int(k, n, I)
  integer,      intent(in)    :: k, n
  integer,      intent(inout) :: I(*)

  integer :: ierr


  call MPI_BCAST(I, n, MPI_INTEGER, k, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "aux_broadcast_int")

  end subroutine aux_broadcast_int
  !-----------------------------------------------------------------------------
  subroutine broadcast_int1(this, I)
  class(process), intent(in)    :: this
  integer,        intent(inout) :: I(:)
  call aux_broadcast_int(this%id, size(I), I)
  end subroutine broadcast_int1
  !-----------------------------------------------------------------------------
  subroutine broadcast_int2(this, I)
  class(process), intent(in)    :: this
  integer,        intent(inout) :: I(:,:)
  call aux_broadcast_int(this%id, size(I), I)
  end subroutine broadcast_int2
  !-----------------------------------------------------------------------------
  subroutine broadcast_int3(this, I)
  class(process), intent(in)    :: this
  integer,        intent(inout) :: I(:,:,:)
  call aux_broadcast_int(this%id, size(I), I)
  end subroutine broadcast_int3
  !-----------------------------------------------------------------------------
  subroutine broadcast_int4(this, I)
  class(process), intent(in)    :: this
  integer,        intent(inout) :: I(:,:,:,:)
  call aux_broadcast_int(this%id, size(I), I)
  end subroutine broadcast_int4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine aux_broadcast_real(k, n, R)
  integer,      intent(in)    :: k, n
  real(real64), intent(inout) :: R(*)

  integer :: ierr


  call MPI_BCAST(R, n, MPI_REAL8, k, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "aux_broadcast_real")

  end subroutine aux_broadcast_real
  !-----------------------------------------------------------------------------
  subroutine broadcast_real1(this, R)
  class(process), intent(in)    :: this
  real(real64),   intent(inout) :: R(:)
  call aux_broadcast_real(this%id, size(R), R)
  end subroutine broadcast_real1
  !-----------------------------------------------------------------------------
  subroutine broadcast_real2(this, R)
  class(process), intent(in)    :: this
  real(real64),   intent(inout) :: R(:,:)
  call aux_broadcast_real(this%id, size(R), R)
  end subroutine broadcast_real2
  !-----------------------------------------------------------------------------
  subroutine broadcast_real3(this, R)
  class(process), intent(in)    :: this
  real(real64),   intent(inout) :: R(:,:,:)
  call aux_broadcast_real(this%id, size(R), R)
  end subroutine broadcast_real3
  !-----------------------------------------------------------------------------
  subroutine broadcast_real4(this, R)
  class(process), intent(in)    :: this
  real(real64),   intent(inout) :: R(:,:,:,:)
  call aux_broadcast_real(this%id, size(R), R)
  end subroutine broadcast_real4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine aux_broadcast_complex(k, n, C)
  integer,          intent(in)    :: k, n
  complex(kind=dp), intent(inout) :: C(*)

  integer :: ierr


  call MPI_BCAST(C, n, MPI_DOUBLE_COMPLEX, k, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "aux_broadcast_complex")

  end subroutine aux_broadcast_complex
  !-----------------------------------------------------------------------------
  subroutine broadcast_complex1(this, C)
  class(process),   intent(in)    :: this
  complex(kind=dp), intent(inout) :: C(:)
  call aux_broadcast_complex(this%id, size(C), C)
  end subroutine broadcast_complex1
  !-----------------------------------------------------------------------------
  subroutine broadcast_complex2(this, C)
  class(process),   intent(in)    :: this
  complex(kind=dp), intent(inout) :: C(:,:)
  call aux_broadcast_complex(this%id, size(C), C)
  end subroutine broadcast_complex2
  !-----------------------------------------------------------------------------
  subroutine broadcast_complex3(this, C)
  class(process),   intent(in)    :: this
  complex(kind=dp), intent(inout) :: C(:,:,:)
  call aux_broadcast_complex(this%id, size(C), C)
  end subroutine broadcast_complex3
  !-----------------------------------------------------------------------------
  subroutine broadcast_complex4(this, C)
  class(process),   intent(in)    :: this
  complex(kind=dp), intent(inout) :: C(:,:,:,:)
  call aux_broadcast_complex(this%id, size(C), C)
  end subroutine broadcast_complex4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine aux_broadcast_logical(k, n, L)
  integer,          intent(in)    :: k, n
  logical,          intent(inout) :: L(*)

  integer :: ierr


  call MPI_BCAST(L, n, MPI_LOGICAL, k, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "aux_broadcast_logical")

  end subroutine aux_broadcast_logical
  !-----------------------------------------------------------------------------
  subroutine broadcast_logical1(this, L)
  class(process),   intent(in)    :: this
  logical,          intent(inout) :: L(:)
  call aux_broadcast_logical(this%id, size(L), L)
  end subroutine broadcast_logical1
  !-----------------------------------------------------------------------------
  subroutine broadcast_logical2(this, L)
  class(process),   intent(in)    :: this
  logical,          intent(inout) :: L(:,:)
  call aux_broadcast_logical(this%id, size(L), L)
  end subroutine broadcast_logical2
  !-----------------------------------------------------------------------------
  subroutine broadcast_logical3(this, L)
  class(process),   intent(in)    :: this
  logical,          intent(inout) :: L(:,:,:)
  call aux_broadcast_logical(this%id, size(L), L)
  end subroutine broadcast_logical3
  !-----------------------------------------------------------------------------
  subroutine broadcast_logical4(this, L)
  class(process),   intent(in)    :: this
  logical,          intent(inout) :: L(:,:,:,:)
  call aux_broadcast_logical(this%id, size(L), L)
  end subroutine broadcast_logical4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine broadcast_string1(this, S)
  class(process),   intent(in)    :: this
  character(len=*), intent(inout) :: S(:)

  integer :: i


  do i=1,size(S)
     call this%broadcast_string(S(i))
  enddo

  end subroutine broadcast_string1
  !-----------------------------------------------------------------------------
! broadcast --------------------------------------------------------------------
! type process =================================================================


! module procedures:
  !---------------------------------------------------------------------
  subroutine moose_mpi_init()
  integer :: ierr


  call MPI_INIT(ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "moose_mpi_init")

  call moose_mpi_setup()

  end subroutine moose_mpi_init
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine moose_mpi_setup()
  integer :: ierr, i


  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "moose_mpi_init size")

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "moose_mpi_init rank")

  allocate (proc(0:nproc-1))
  do i=0,nproc-1
     proc(i)%id = i
  enddo

  end subroutine moose_mpi_setup
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine moose_mpi_finalize()
  integer :: ierr


  deallocate (proc)
  call mpi_barrier_world()
  call MPI_FINALIZE(ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "moose_mpi_finalize")

  end subroutine moose_mpi_finalize
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine moose_mpi_error(ierr, location)
  !
  ! display error message after MPI call
  !
  integer,          intent(in) :: ierr
  character(len=*), intent(in) :: location


  write (6, 9000) ierr, location;   stop
 9000 format("MPI error ",i0," in ",a,"!")

  end subroutine moose_mpi_error
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine mpi_barrier_world()
  !
  ! blocks until all processes have reached this routine
  !
  integer :: ierr


  call MPI_BARRIER (MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "mpi_barrier_world")

  end subroutine mpi_barrier_world
  !---------------------------------------------------------------------


! moose_mpi_sum (based on MPI_ALLREDUCE) ---------------------------------------
  !-----------------------------------------------------------------------------
  subroutine sum_int(I)
  integer, intent(inout) :: I

  integer :: ierr, isum


  call mpi_barrier_world()
  call MPI_ALLREDUCE(I, isum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "sum_int")
  I = isum

  end subroutine sum_int
  !-----------------------------------------------------------------------------
  subroutine aux_sum_integer_array(n, I)
  integer, intent(in   ) :: n
  integer, intent(inout) :: I(*)

  integer, allocatable :: itmp(:)
  integer :: i1, i2, ierr, m


  call mpi_barrier_world()
  m = min(n, chunk_size)
  allocate (itmp(m))
  do i1=1,n,chunk_size
     i2 = min(n, i1+chunk_size-1)
     itmp(1:i2-i1+1) = I(i1:i2)
     call MPI_ALLREDUCE(itmp, I(i1:i2), i2-i1+1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     if (ierr /= 0) call moose_mpi_error(ierr, "allreduce_integer_array")
  enddo
  deallocate (itmp)

  end subroutine aux_sum_integer_array
  !-----------------------------------------------------------------------------
  subroutine sum_integer_array1(I)
  integer, intent(inout) :: I(:)
  call aux_sum_integer_array(size(I), I)
  end subroutine sum_integer_array1
  !-----------------------------------------------------------------------------
  subroutine sum_integer_array2(I)
  integer, intent(inout) :: I(:,:)
  call aux_sum_integer_array(size(I), I)
  end subroutine sum_integer_array2
  !-----------------------------------------------------------------------------
  subroutine sum_integer_array3(I)
  integer, intent(inout) :: I(:,:,:)
  call aux_sum_integer_array(size(I), I)
  end subroutine sum_integer_array3
  !-----------------------------------------------------------------------------
  subroutine sum_integer_array4(I)
  integer, intent(inout) :: I(:,:,:,:)
  call aux_sum_integer_array(size(I), I)
  end subroutine sum_integer_array4
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine sum_real(R)
  real(real64), intent(inout) :: R

  real(real64) :: rsum
  integer :: ierr


  call mpi_barrier_world()
  call MPI_ALLREDUCE(R, rsum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) call moose_mpi_error(ierr, "sum_real")
  R = rsum

  end subroutine sum_real
  !-----------------------------------------------------------------------------
  subroutine aux_sum_real_array(n, R)
  integer,      intent(in)    :: n
  real(real64), intent(inout) :: R(*)

  real(real64), allocatable   :: rtmp(:)

  integer :: i1, i2, ierr, m


  call mpi_barrier_world()
  m = min(n, chunk_size)
  allocate (rtmp(m))
  do i1=1,n,chunk_size
     i2 = min(n, i1+chunk_size-1)
     rtmp(1:i2-i1+1) = R(i1:i2)
     call MPI_ALLREDUCE(rtmp, R(i1:i2), i2-i1+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     if (ierr /= 0) call moose_mpi_error(ierr, "allreduce_real_array")
  enddo
  deallocate (rtmp)

  end subroutine aux_sum_real_array
  !-----------------------------------------------------------------------------
  subroutine sum_real_array1(R)
  real(real64), intent(inout) :: R(:)
  call aux_sum_real_array(size(R), R)
  end subroutine sum_real_array1
  !-----------------------------------------------------------------------------
  subroutine sum_real_array2(R)
  real(real64), intent(inout) :: R(:,:)
  call aux_sum_real_array(size(R), R)
  end subroutine sum_real_array2
  !-----------------------------------------------------------------------------
  subroutine sum_real_array3(R)
  real(real64), intent(inout) :: R(:,:,:)
  call aux_sum_real_array(size(R), R)
  end subroutine sum_real_array3
  !-----------------------------------------------------------------------------
  subroutine sum_real_array4(R)
  real(real64), intent(inout) :: R(:,:,:,:)
  call aux_sum_real_array(size(R), R)
  end subroutine sum_real_array4
  !-----------------------------------------------------------------------------
! moose_mpi_sum ----------------------------------------------------------------


  !---------------------------------------------------------------------
  function moose_mpi_range(n) result(i)
  !
  ! partition range [0,n-1] into intervals [i(1),i(2)] over all processes
  !
  integer, intent(in) :: n
  integer             :: i(1:2)


  i(1) = rank * n / nproc
  i(2) = min((rank+1) * n / nproc, n) - 1

  end function moose_mpi_range
  !---------------------------------------------------------------------

end module moose_mpi
