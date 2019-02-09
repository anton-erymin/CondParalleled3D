module GlobalDefinitions
use MPI

implicit none

! ������������� ��������
integer id
! ���-�� ��������� � �������
integer numProcesses
! ��������� ���������� ������� MPI
integer ierr
integer status(MPI_STATUS_SIZE)

! �������� ��������
real(8), parameter :: eps = 1.D-6

! ���������� ��������� �����
integer nx, ny, nz, n
! �����
real(8), allocatable, dimension(:) :: x, y, z, xu, yu, zu

! ��������� �������
real(8) xl, xr, yl, yr, zl, zr

real(8), allocatable, dimension(:,:,:) :: u, u0, u1, gradu, gradu0
real(8), allocatable, dimension(:,:,:) :: rho, gami, gamj, gamk
real(8), allocatable, dimension(:,:,:) :: con, ap, aps, b
real(8), allocatable, dimension(:,:,:) :: aip, aim, ajm, ajp, akm, akp
real(8), allocatable, dimension(:,:,:) :: vu, vv, vw

real(8) t, norm, maxgrad

end