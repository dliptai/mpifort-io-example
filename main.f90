program Main
  use mpi
  use utils

  implicit none

  integer :: ierr, array_view
  integer :: i,j,k, ranki
  real(8), dimension(3) :: header, header_old
  character(len=*), parameter :: filename = 'output.dat'

  real(8), allocatable, dimension(:,:,:) :: x1, x2

  call init_mpi

  header = [1.,2.,3.]
  is_global = 1 ; ie_global = 8
  js_global = 1 ; je_global = 8
  ks_global = 1 ; ke_global = 8

  call decompose_domain
  call set_array_view(array_view)

  ! allocate and initialise array
  allocate(x1(is:ie,js:je,ks:ke))
  allocate(x2,mold=x1)

  do i = is,ie
    do j = js,je
      do k = ks,ke
        x1(i,j,k) = 100*i + 10*j + k
        x2(i,j,k) = -(100*i + 10*j + k)
      end do
    end do
  end do

  call write_header(filename, header)
  call write_array_mpi(filename, x1, array_view)

  ! put junk in array
  x1 = -1.0

  ! -------------- Now read back the same file -----------------------------------------

  header_old = header
  call read_header(filename, header)

  do i = 1,3
    if (abs(header(i) - header_old(i)) > 1e-14) print*, 'rank', myrank, 'Header mismatch: ', header(i), header_old(i)
  end do

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call read_array_mpi(filename, x1, array_view)

  do ranki = 0, nprocs-1
    if (myrank == ranki) then
      do i = is,ie
        do j = js,je
          do k = ks,ke
            if (abs(x1(i,j,k) - (100*i + 10*j + k)) > 1e-14) then
              print*, 'Rank ', myrank, ' has a problem at ', i, j, k, x1(i,j,k)
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
