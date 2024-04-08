program Main
  use mpi

  implicit none

  integer :: nx, ny, nz
  integer :: ierr
  integer :: dims(3)
  logical :: periods(3)
  integer :: mycoords(3)
  integer :: is, ie, js, je, ks, ke
  integer :: is_global, ie_global, js_global, je_global, ks_global, ke_global
  integer :: cart_comm, myrank, nprocs
  integer :: i,j,k
  integer :: sizes(3), subsizes(3), starts(3)
  integer :: count, file, filetype, status(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer :: ranki, un
  real(8) :: header, header_new
  integer, parameter :: record_marker = 4

  real(8), allocatable :: xyz(:,:,:)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  header = 88.0
  disp = 2*record_marker + sizeof(header)

  is_global = 1
  ie_global = 8
  js_global = 1
  je_global = 8
  ks_global = 1
  ke_global = 8

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
    print*, 'Global domain: ', is_global, ie_global, js_global, je_global, ks_global, ke_global
  endif

  print*, 'Rank ', myrank, ' has domain ', is, ie, js, je, ks, ke

  sizes = [nx,ny,nz]
  subsizes = [ie-is+1,je-js+1,ke-ks+1]
  count = subsizes(1)*subsizes(2)*subsizes(3)
  starts = [is-1,js-1,ks-1]

  call mpi_type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, filetype, ierr)
  call mpi_type_commit(filetype, ierr)

  allocate(xyz(is:ie,js:je,ks:ke))

  ! initialise array
  do i = is,ie
    do j = js,je
      do k = ks,ke
        xyz(i,j,k) = 100*i + 10*j + k
      end do
    end do
  end do

  ! Write header
  if (myrank == 0) then
    open(newunit=un, file='output.dat', status='replace', form='unformatted', action='write')
    write(un) header
    close(un)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Write the array
  call mpi_file_open(MPI_COMM_WORLD, 'output.dat', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file, ierr)
  call mpi_file_set_view(file, disp, MPI_DOUBLE_PRECISION, filetype, 'native', MPI_INFO_NULL, ierr)

  call mpi_file_write_all(file, xyz, count, MPI_DOUBLE_PRECISION, status, ierr)
  call mpi_file_close(file, ierr)

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! put junk in array
  xyz = -1.0

  ! -------------- Now read back the same file -----------------------------------------

  header_new = -1.0

  open(newunit=un, file='output.dat', status='old', form='unformatted', action='read')
  read(un) header_new
  close(un)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  if (abs(header - header_new) > 1e-14) then
    print*, 'rank', myrank, 'Header mismatch: ', header, header_new
  else
    print*, 'rank', myrank, 'Header matches: ', header, header_new
  end if

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call mpi_file_open(MPI_COMM_WORLD, 'output.dat', MPI_MODE_RDONLY, MPI_INFO_NULL, file, ierr)
  call mpi_file_set_view(file, disp, MPI_DOUBLE_PRECISION, filetype, 'native', MPI_INFO_NULL, ierr)
  call mpi_file_read_all(file, xyz, count, MPI_DOUBLE_PRECISION, status, ierr)
  call mpi_file_close(file, ierr)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  do ranki = 0, nprocs-1
    if (myrank == ranki) then
      do i = is,ie
        do j = js,je
          do k = ks,ke
            if (abs(xyz(i,j,k) - (100*i + 10*j + k)) > 1e-14) then
              print*, 'Rank ', myrank, ' has a problem at ', i, j, k, xyz(i,j,k)
            end if
          end do
        end do
      end do
    endif
    call mpi_barrier(cart_comm, ierr)
  end do

  call MPI_FINALIZE(ierr)

end program Main
