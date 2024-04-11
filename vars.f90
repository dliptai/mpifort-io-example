module vars
  implicit none

  integer :: nx, ny, nz
  integer :: is, ie, js, je, ks, ke
  integer :: is_global, ie_global, js_global, je_global, ks_global, ke_global

  integer :: tn
  real(8) :: time

  real(8), allocatable :: d(:,:,:)
  real(8), allocatable :: v1(:,:,:)
  real(8), allocatable :: v2(:,:,:)
  real(8), allocatable :: v3(:,:,:)
  real(8), allocatable :: e(:,:,:)

  integer :: type_mpi_array

end module vars
