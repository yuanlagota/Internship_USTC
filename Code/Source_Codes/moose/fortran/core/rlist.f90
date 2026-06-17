module moose_rlist
  use iso_fortran_env
  use moose_txtio
  use moose_table
  implicit none
  private


  ! list of linear array elements of size ndim .................................
  type, extends(txtio), public :: rlist
     ! implement list as resizable table
     type(table), pointer, private :: implementation

     ! internal parameter for table size increments
     integer, private :: chunk_size

     ! indices for first and last element
     integer, private :: first, last

     ! size of array elements (THIS VALUE MUST NOT BE CHANGED AFTER INITIALIZATION)
     integer :: ndim

     contains
     procedure :: allreduce
     procedure :: send
     procedure :: broadcast
     procedure :: free

     procedure :: element     ! return pointer to i-th element
     procedure :: nelements   ! return number of elements in list

     procedure :: values      ! return pointer to array with all elements
     procedure :: column      ! return pointer to j-th component of all elements
     procedure :: columns     ! return pointer to array with components j1-j2 of all elements

     generic :: append => append_element, append_list
     procedure :: append_element      ! append element to list
     procedure :: append_list         ! append another list to list
     generic :: prepend => prepend_element, prepend_list
     procedure :: prepend_element     ! prepend element to list
     procedure :: prepend_list        ! prepend another list to list
     procedure :: pop         ! return last element and remove it from list
     procedure :: drop        ! drop i-th element from list
     procedure :: clear       ! remove all elements from list

     procedure :: sort        ! sort list with respect to values in j-th component
     procedure :: reverse     ! reverse elements in list

     procedure :: write_formatted
  end type rlist


  interface rlist
     procedure :: new
     procedure :: load
  end interface



  ! extended types for arrays of size 2 and 3
  type, extends(rlist), public :: rlist2
  end type rlist2
  type, extends(rlist), public :: rlist3
  end type rlist3

  interface rlist2
     procedure :: new2
  end interface

  interface rlist3
     procedure :: new3
  end interface



  public :: &
     readnc_rlist, &
     recv_rlist


  contains
  !---------------------------------------------------------------------


! constructors:
  !---------------------------------------------------------------------
  function new(ndim, chunk_size) result(L)
  integer, intent(in) :: ndim
  integer, intent(in), optional :: chunk_size
  type(rlist)         :: L


  call init_txtio(L, "rlist")
  L%chunk_size = 1024;   if (present(chunk_size)) L%chunk_size = chunk_size

  allocate (L%implementation, source=table(ndim, L%chunk_size, lbounds=[1,0]))
  L%ndim  = ndim
  L%first = 0
  L%last  = -1

  end function new
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function new2(chunk_size) result(L)
  integer, intent(in), optional :: chunk_size
  type(rlist2) :: L


  L%rlist = rlist(2, chunk_size)

  end function new2
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function new3(chunk_size) result(L)
  integer, intent(in), optional :: chunk_size
  type(rlist3) :: L


  L%rlist = rlist(3, chunk_size)

  end function new3
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function load(filename, ndim, chunk_size) result(this)
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: ndim, chunk_size
  type(rlist)                  :: this

  type(table) :: T
  integer :: n, m


  T    = table(filename, ndim, .true., [1, 0])
  this = new(T%rows(), chunk_size)

  n = T%columns()
  if (this%implementation%columns() < n) then
     m = n / this%chunk_size;   if (mod(n, this%chunk_size) /= 0) m = m + 1
     call this%implementation%resize_columns(m * this%chunk_size)
  endif

  this%implementation%values(:,0:n-1) = T%values
  this%last = n-1

  end function load
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function readnc_rlist(nc, varname) result(this)
  use moose_error
  use moose_utils, only: str
  use moose_netcdf
  type(netcdf_dataset), intent(in) :: nc
  character(len=*),     intent(in) :: varname
  type(rlist)                      :: this

  integer :: varid, ndim, ndims, nelements, dimids(2)


  varid = nc%inq_varid(varname)
  call nc%inquire_variable(varid, ndims=ndims)
  if (ndims < 1  .or.  ndims > 2) then
     call ERROR("unexpected ndims = "//str(ndims)//" for variable '"//varname//"'")
  endif


  call nc%inquire_variable(varid, dimids=dimids)
  if (ndims == 1) then
     ndim = 1
     call nc%inquire_dimension(dimids(1), len=nelements)
  else
     call nc%inquire_dimension(dimids(1), len=ndim)
     call nc%inquire_dimension(dimids(2), len=nelements)
  endif


  this = new(ndim, nelements)
  this%last = nelements - 1
  call nc%get_var(varname, this%implementation%values)

  end function readnc_rlist
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function recv_rlist(from) result(L)
  use moose_mpi
  integer, intent(in) :: from
  type(rlist)         :: L

  integer :: ndim, chunk_size, nelements, chunks


  call proc(from)%recv(nelements)
  call proc(from)%recv(ndim)
  call proc(from)%recv(chunk_size)
  L = new(ndim, chunk_size)

  chunks = ceiling(1.d0 * nelements / chunk_size)
  call L%implementation%resize(ndim, chunks * chunk_size)
  if (nelements > 0) call proc(from)%recv(L%implementation%values(:,0:nelements-1))
  L%first = 0
  L%last  = nelements - 1

  end function recv_rlist
  !---------------------------------------------------------------------


! type-bound procedures:
  !---------------------------------------------------------------------
  subroutine allreduce(this)
  use moose_mpi
  class(rlist), intent(inout) :: this

  integer :: nelements(0:nproc-1), i, m, chunks


  if (nproc == 1) return

  ! collect number of elements from each process
  call proc(0)%send(this%nelements())
  if (rank == 0) then
     do i=0,nproc-1
        call proc(i)%recv(nelements(i))
     enddo

     ! prepare implementation on process 0 to store all elements
     this%last = this%last + sum(nelements) - nelements(0)
     m         = this%last - ubound(this%implementation%values, 2)
     if (m > 0) then
        chunks = (m + this%chunk_size - 1) / this%chunk_size
        call this%implementation%append_columns(chunks * this%chunk_size)
     endif
  endif


  ! send data to first process
  if (rank > 0) then
     call proc(0)%send(this%implementation%values(:,this%first:this%last))
  endif


  ! collect data from each process
  if (rank == 0) then
     m = this%first
     do i=1,nproc-1
        m = m + nelements(i-1)
        call proc(i)%recv(this%implementation%values(:,m:m+nelements(i)-1))
     enddo
  endif

  end subroutine allreduce
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine send(this, to)
  use moose_mpi
  class(rlist), intent(in) :: this
  integer,      intent(in) :: to


  call proc(to)%send(this%nelements())
  call proc(to)%send(this%ndim)
  call proc(to)%send(this%chunk_size)
  if (this%nelements() == 0) return
  call proc(to)%send(this%implementation%values(:,this%first:this%last))

  end subroutine send
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(rlist), intent(inout) :: this


  call this%txtio_broadcast()
  call proc(0)%broadcast(this%chunk_size)
  call proc(0)%broadcast(this%first)
  call proc(0)%broadcast(this%last)
  call proc(0)%broadcast(this%ndim)
  if (rank > 0) allocate (this%implementation)
  call this%implementation%broadcast()

  end subroutine broadcast
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine free(this)
  class(rlist), intent(inout) :: this


  call this%implementation%free()
  call this%txtio_free()

  end subroutine free
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(rlist),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%implementation%write_selection(unit, vlist, iostat, iomsg, &
     columns=[this%first, this%last], transposed=.true.)

  end subroutine write_formatted
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function element(this, i) result(x)
  !
  ! i-th element in list for 0 <= i < this%nelements
  !    or this%nelements+i-th element for -this%nelements <= i < 0
  ! returns:
  !    element(1:ndim)
  !
  class(rlist), intent(in) :: this
  integer,      intent(in) :: i
  real(real64), pointer    :: x(:)


  nullify(x)
  if (i >= 0  .and.  i <= this%last-this%first) then
     x => this%implementation%values(:,this%first+i)
  elseif (i < 0 .and. i >= this%first-this%last-1) then
     x => this%implementation%values(:,this%last+1+i)
  endif

  end function element
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  pure function nelements(this) result(n)
  class(rlist), intent(in) :: this
  integer                  :: n


  n = this%last - this%first + 1

  end function nelements
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function values(this)
  !
  ! pointer to array with all elements
  ! returns:
  !    values(1:ndim , 0:nelements-1)
  !
  class(rlist), intent(in) :: this
  real(real64), pointer    :: values(:,:)


  values(1:,0:) => this%implementation%values(:,this%first:this%last)

  end function values
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function column(this, j)
  !
  ! pointer to j-th component of all elements
  ! returns:
  !    column(0:nelements-1)
  !
  class(rlist), intent(in) :: this
  integer,      intent(in) :: j
  real(real64), pointer    :: column(:)


  nullify(column)
  if (j < 1  .or.  j > this%ndim) return
  column(0:) => this%implementation%values(j,this%first:this%last)

  end function column
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function columns(this, j1, j2)
  !
  ! pointer to array with components j1-j2 of all elements
  ! returns:
  !    columns(j1:j2, 0:nelements-1)
  !
  class(rlist), intent(in) :: this
  integer,      intent(in) :: j1, j2
  real(real64), pointer    :: columns(:,:)


  columns(1:,0:) => this%implementation%values(j1:j2,this%first:this%last)

  end function columns
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine append_element(this, x, direction)
  class(rlist), intent(inout) :: this
  real(real64), intent(in)    :: x(this%ndim)
  integer,      intent(in), optional :: direction


  ! prepend values instead, if requested
  if (present(direction)) then
     if (direction == -1) then
        call this%prepend(x)
        return
     endif
  endif


  ! resize table if necessary
  if (this%last == ubound(this%implementation%values, 2)) then
     call this%implementation%append_columns(this%chunk_size)
  endif


  ! add data to list
  this%last = this%last + 1
  this%implementation%values(:,this%last) = x

  end subroutine append_element
  !---------------------------------------------------------------------
  subroutine append_list(this, L)
  class(rlist), intent(inout) :: this
  class(rlist), intent(in   ) :: L

  integer :: n, nadd, missing


  ! verify compatibility
  if (this%ndim /= L%ndim) then
     print 9000, L%ndim, this%ndim;   stop
  endif
 9000 format("ERROR in rlist%prepend: incompatible ndim ",i0," /= ",i0)


  ! update index for new last element
  n = L%nelements()
  this%last = this%last + n


  ! resize table if necessary
  missing = this%last - ubound(this%implementation%values, 2)
  if (missing > 0) then
     nadd = (missing / this%chunk_size + 1) * this%chunk_size
     call this%implementation%append_columns(nadd)
  endif

  ! add data to list
  this%implementation%values(:,this%last-n+1:this%last) = L%implementation%values(:,L%first:L%last)

  end subroutine append_list
  !---------------------------------------------------------------------
  subroutine prepend_element(this, x)
  class(rlist), intent(inout) :: this
  real(real64), intent(in)    :: x(this%ndim)


  ! resize table if necessary
  if (this%first == lbound(this%implementation%values, 2)) then
     call this%implementation%prepend_columns(this%chunk_size)
  endif


  ! add data to list
  this%first = this%first - 1
  this%implementation%values(:,this%first) = x

  end subroutine prepend_element
  !---------------------------------------------------------------------
  subroutine prepend_list(this, L)
  class(rlist), intent(inout) :: this
  class(rlist), intent(in   ) :: L

  integer :: n, nadd, missing


  ! verify compatibility
  if (this%ndim /= L%ndim) then
     print 9000, L%ndim, this%ndim;   stop
  endif
 9000 format("ERROR in rlist%prepend: incompatible ndim ",i0," /= ",i0)


  ! update index for new first element
  n = L%nelements()
  this%first = this%first - n


  ! resize table if necessary
  missing = lbound(this%implementation%values, 2) - this%first
  if (missing > 0) then
     nadd = (missing / this%chunk_size + 1) * this%chunk_size
     call this%implementation%prepend_columns(nadd)
  endif


  ! add data to list
  this%implementation%values(:,this%first:this%first+n-1) = L%implementation%values(:,L%first:L%last)

  end subroutine prepend_list
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  function pop(this) result(x)
  !
  ! return last element and remove it from list
  !
  class(rlist), intent(inout) :: this
  real(real64)                :: x(this%ndim)


  x = this%implementation%values(:,this%last)
  call this%drop(this%nelements()-1)

  end function pop
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  pure subroutine drop(this, i)
  !
  ! drop i-th element from list
  !
  class(rlist), intent(inout) :: this
  integer,      intent(in)    :: i

  integer :: ii
  associate (x => this%implementation%values)


  ii = i + this%first
  if (ii < this%first  .or.  ii > this%last) return
  ! remove last remaining element from list
  if (this%first == this%last) then
     call this%implementation%resize(this%ndim, this%chunk_size, lcol=0)
     this%first = 0
     this%last  = -1

  ! drop i-th element by moving elements > i up in list
  else
     if (ii < this%last) x(:,ii:this%last-1) = x(:,ii+1:this%last)
     this%last = this%last - 1

     if (ubound(x, 2) - this%last > this%chunk_size) then
        call this%implementation%remove_columns(this%chunk_size)
     endif
  endif

  end associate
  end subroutine drop
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine clear(this)
  class(rlist), intent(inout) :: this


  call this%implementation%resize_columns(this%chunk_size)
  this%first = 0
  this%last  = -1

  end subroutine clear
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine sort(this, column)
  class(rlist), intent(inout) :: this
  integer,      intent(in   ) :: column


  call this%implementation%sort(column, axis=2, irange=[this%first, this%last])

  end subroutine sort
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  subroutine reverse(this)
  class(rlist), intent(inout) :: this


  call this%implementation%reverse_columns(this%first, this%last)

  end subroutine reverse
  !---------------------------------------------------------------------

end module moose_rlist
