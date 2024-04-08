program Main
  use mpi
  use utils

  implicit none

  integer :: ierr, array_view
  integer :: i,j,k, ranki
  real(8) :: header, header_new(1)
  character(len=*), parameter :: filename = 'output.dat'

  real(8), allocatable :: xyz(:,:,:)

  call init_mpi

  header = 88.0
  is_global = 1 ; ie_global = 8
  js_global = 1 ; je_global = 8
  ks_global = 1 ; ke_global = 8

  call decompose_domain
  call set_array_view(array_view)

  ! allocate and initialise array
  allocate(xyz(is:ie,js:je,ks:ke))
  do i = is,ie
    do j = js,je
      do k = ks,ke
        xyz(i,j,k) = 100*i + 10*j + k
      end do
    end do
  end do

  call write_header(filename, [header])
  call write_array_mpi(filename, xyz, array_view)

  ! put junk in array
  xyz = -1.0

  ! -------------- Now read back the same file -----------------------------------------

  header_new = -1.0
  call read_header(filename, header_new)

  if (abs(header - header_new(1)) > 1e-14) print*, 'rank', myrank, 'Header mismatch: ', header, header_new(1)

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call read_array_mpi(filename, xyz, array_view)

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
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end do

  if (myrank==0) then
    print*, 'Great successsss'
    print*, 'bytes read', bytes_read, 'bytes written', bytes_written
  endif

  call MPI_FINALIZE(ierr)

end program Main
