module utils
  use mpi
  use vars

  implicit none

  integer :: myrank, nprocs

  contains

  subroutine init_mpi()
    integer :: ierr

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  end subroutine init_mpi

  subroutine decompose_domain()
    integer :: dims(3), mycoords(3), cart_comm, ierr
    logical :: periods(3)

    nx = ie_global - is_global + 1
    ny = je_global - js_global + 1
    nz = ke_global - ks_global + 1

    dims = [nprocs, 1, 1] ! Trivial decomposition, for now... TODO
    periods = [.true., .true., .true.] ! Always set to periodic and allow boundary conditions to override

    call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, .false., cart_comm, ierr)
    call MPI_CART_COORDS(cart_comm, myrank, 3, mycoords, ierr)

    ! Compute local domain, including offset if starting index is not 1
    is = mycoords(1) * (nx / dims(1)) + 1 - (is_global - 1)
    ie = (mycoords(1) + 1) * (nx / dims(1)) - (is_global - 1)
    js = mycoords(2) * (ny / dims(2)) + 1 - (js_global - 1)
    je = (mycoords(2) + 1) * (ny / dims(2)) - (js_global - 1)
    ks = mycoords(3) * (nz / dims(3)) + 1 - (ks_global - 1)
    ke = (mycoords(3) + 1) * (nz / dims(3)) - (ks_global - 1)

    ! Special treatment of final proc, in case nx/ny/nz is not divisible by nprocs
    if (mycoords(1) == dims(1) - 1) ie = nx - (is_global - 1)
    if (mycoords(2) == dims(2) - 1) je = ny - (js_global - 1)
    if (mycoords(3) == dims(3) - 1) ke = nz - (ks_global - 1)

    if (myrank == 0) then
      write(*, '(A, 6I6)') 'Global domain:     ', is_global, ie_global, js_global, je_global, ks_global, ke_global
    endif

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    write(*, '(A, I3, A, 6I6)') 'Rank', myrank, ' has domain ', is, ie, js, je, ks, ke

  end subroutine decompose_domain

  subroutine create_type_mpi_array
    integer :: ierr, sizes(3), subsizes(3), starts(3)

    sizes = [nx,ny,nz]
    subsizes = [ie-is+1,je-js+1,ke-ks+1]
    starts = [is-1,js-1,ks-1]

    call mpi_type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_mpi_array, ierr)
    call mpi_type_commit(type_mpi_array, ierr)

  end subroutine create_type_mpi_array

  ! subroutine to try reading the header in sequential format. If it succeeds it should be a legacy file
  subroutine is_file_legacy(filename, legacy)
    character(len=*), intent(in) :: filename
    logical, intent(out) :: legacy
    real(8) :: header
    integer :: ierr, un

    ! Try reading the header in sequential format
    open(newunit=un, file=filename, status='old', form='unformatted', action='read', access='sequential')
    read(un, iostat=ierr) header
    close(un)

    if (ierr /= 0) then
      legacy = .false.
    else
      legacy = .true.
    end if

    if (myrank == 0) print *, 'Legacy file =', legacy

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine is_file_legacy

  ! a function to determine the total number of elements of a 3D array
  function total_elements(array) result(n)
    integer :: n
    real(8), intent(in) :: array(:,:,:)

    n = size(array,1) * size(array,2) * size(array,3)

  end function total_elements

  subroutine get_file_end(fh, end_bytes)
    integer, intent(in) :: fh
    integer :: ierr
    integer(kind=MPI_OFFSET_KIND), intent(out) :: end_bytes
    ! integer(kind=MPI_OFFSET_KIND) :: disp

    ! This may or may not be better than the commented bit below...
    call mpi_file_sync(fh, ierr)
    call MPI_File_get_size(fh, end_bytes, ierr)

    ! ! Reset view to default ("absolute" view in bytes)
    ! call mpi_file_set_view(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)

    ! ! Move individual file positions to absolute end of file
    ! call mpi_file_seek(fh, 0_MPI_OFFSET_KIND, MPI_SEEK_END, ierr)
    ! call mpi_file_get_position(fh, disp, ierr)
    ! call MPI_File_get_byte_offset(fh, disp, end_bytes, ierr)

  end subroutine get_file_end

  subroutine open_file_read(filename, fh)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: fh
    integer :: ierr
    call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine open_file_read

  subroutine open_file_write(filename, fh)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: fh
    integer :: ierr
    call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR + MPI_MODE_APPEND, MPI_INFO_NULL, fh, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine open_file_write

  subroutine close_file(fh)
    integer, intent(inout) :: fh
    integer :: ierr
    call mpi_file_close(fh, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine close_file

end module utils
