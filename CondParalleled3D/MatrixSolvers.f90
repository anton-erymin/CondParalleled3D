module MatrixSolvers
implicit none
contains

! –ешение системы уравнений методом √аусса-«ейдел¤
!subroutine Sor(u, aim, aip, ajm, ajp, ap, con)
!real(8), allocatable, dimension(:,:) :: u, aim, aip, ajm, ajp, ap, con
!integer i, j
!    do j = 2, size(u(:, 1)) - 1
!        do i = 2, size(u(1, :)) - 1
!            u(i, j) = (aim(i, j) * u(i - 1, j) + aip(i, j) * u(i + 1, j) + &
!                       ajm(i, j) * u(i, j - 1) + ajp(i, j) * u(i, j + 1) + con(i, j)) / ap(i, j)
!        end do
!    end do
!
!end subroutine Sor


! –ешение системы уравнений методом переменных направлений
subroutine CD(u, aim, aip, ajm, ajp, akm, akp, ap, con)
real(8), allocatable, dimension(:,:,:) :: u, aim, aip, ajm, ajp, akm, akp, ap, con
integer i, j, k, nx, ny, nz
real(8), allocatable, dimension(:) :: a, b, c, d, x, a2, b2, c2, d2, x2

    nx = size(u(:, 1, 1))
    ny = size(u(1, :, 1))
    nz = size(u(1, 1, :))
    
    !allocate(a(nx - 2), b(nx - 2), c(nx - 2), d(nx - 2), x(nx - 2))
    allocate(a2(ny - 2), b2(ny - 2), c2(ny - 2), d2(ny - 2), x2(ny - 2))
    
!    do j = 2, ny - 1
!        do k = 2, nz - 1
!        
!            do i = 1, nx - 2
!                a(i) = ap(i + 1, j, k)
!                d(i) = con(i + 1, j, k) + ajm(i + 1, j, k) * u(i + 1, j - 1, k) + ajp(i + 1, j, k) * u(i + 1, j + 1, k)
!                d(i) = d(i) + akm(i + 1, j, k) * u(i + 1, j, k - 1) + akp(i + 1, j, k) * u(i + 1, j, k + 1)
!            end do
!            do i = 2, nx - 2
!                b(i) = -aim(i + 1, j, k)
!            end do
!            do i = 1, nx - 3
!                c(i) = -aip(i + 1, j, k)
!            end do
!            
!            d(1) = d(1) + aim(2, j, k) * u(1, j, k)
!            d(nx - 2) = d(nx - 2) + aip(nx - 1, j, k) * u(nx, j, k)
!            
!            call TDMA(a, b, c, d, x)
!            
!            do i = 2, nx - 1
!                u(i, j, k) = x(i - 1)
!            end do
!            
!        end do
!    end do
    
    
    do i = 2, nx - 1
        do k = 2, nz - 1
        
            do j = 1, ny - 2
                a2(j) = ap(i, j + 1, k)
                d2(j) = con(i, j + 1, k) + aim(i, j + 1, k) * u(i - 1, j + 1, k) + aip(i, j + 1, k) * u(i + 1, j + 1, k)
                d2(j) = d2(j) + akm(i, j + 1, k) * u(i, j + 1, k - 1) + akp(i, j + 1, k) * u(i, j + 1, k + 1)
            end do
            do j = 2, ny - 2
                b2(j) = -ajm(i, j + 1, k)
            end do
            do j = 1, ny - 3
                c2(j) = -ajp(i, j + 1, k)
            end do
            
            d2(1) = d2(1) + ajm(i, 2, k) * u(i, 1, k)
            d2(ny - 2) = d2(ny - 2) + ajp(i, ny - 1, k) * u(i, ny, k)
            
            call TDMA(a2, b2, c2, d2, x2)
            
            do j = 2, ny - 1
                u(i, j, k) = x2(j - 1)
            end do
        end do
    end do
    
    
    do i = 2, nx - 1
        do j = 2, ny - 1
        
            do k = 1, nz - 2
                a2(k) = ap(i, j, k + 1)
                d2(k) = con(i, j, k + 1) + aim(i, j, k + 1) * u(i - 1, j, k + 1) + aip(i, j, k + 1) * u(i + 1, j, k + 1)
                d2(k) = d2(k) + ajm(i, j, k + 1) * u(i, j - 1, k + 1) + ajp(i, j, k + 1) * u(i, j + 1, k + 1)
            end do
            do k = 2, nz - 2
                b2(k) = -akm(i, j, k + 1)
            end do
            do k = 1, nz - 3
                c2(k) = -akp(i, j, k + 1)
            end do
            
            d2(1) = d2(1) + akm(i, j, 2) * u(i, j, 1)
            d2(nz - 2) = d2(nz - 2) + akp(i, j, nz - 1) * u(i, j, nz)
            
            call TDMA(a2, b2, c2, d2, x2)
            
            do k = 2, nz - 1
                u(i, j, k) = x2(k - 1)
            end do
        end do
    end do
    
    
    deallocate(a2, b2, c2, d2, x2)
    

end subroutine CD


! ћетод прогонки
subroutine TDMA(a, b, c, d, x)
real(8), allocatable, dimension(:) :: a, b, c, d, x
real(8), allocatable, dimension(:) :: p, q
integer i, n
real(8) r
    n = size(a)
    allocate(p(size(a)), q(size(a)))

    p(1) = -c(1) / a(1)
    q(1) = d(1) / a(1)
    
    do i = 2, n - 1
        r = a(i) + b(i) * p(i - 1)
        p(i) = -c(i) / r
        q(i) = (d(i) - b(i) * q(i - 1)) / r    
    end do
    
    q(n) = (d(n) - b(n) * q(n - 1)) / (a(n) + b(n) * p(n - 1))
    x(n) = q(n)
    
    i = n - 1
    do while (i >= 1)
        x(i) = p(i) * x(i + 1) + q(i)    
        i = i - 1
    end do
    
    deallocate(p, q)
    
end subroutine TDMA


subroutine CG(u, aim, aip, ajm, ajp, akm, akp, ap, con)
real(8), allocatable, dimension(:,:,:) :: u, aim, aip, ajm, ajp, akm, akp, ap, con
integer i, j, k, nx, ny, nz

real(8), allocatable :: A(:,:), b(:)

    nx = size(u(:, 1, 1))
    ny = size(u(1, :, 1))
    nz = size(u(1, 1, :))

    i = (nx - 2) * (ny - 2) * (nz - 2) 

    !allocate(A(i, i), b(i))


end subroutine CG



end