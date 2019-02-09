module Solver
use GlobalDefinitions
use MatrixSolvers
use GridDecomposition

implicit none

real(8) tmax, dt, grmax, rmax, relax
integer i, j, k, iter, numIterations, nx1, ny1, nz1, alliter

contains

subroutine Solve
    ! Выделяем память под массивы
    allocate(u(nx, ny, nz), u0(nx, ny, nz), u1(nx, ny, nz), rho(nx, ny, nz), gami(nx, ny, nz), gamj(nx, ny, nz), con(nx, ny, nz), b(nx, ny, nz))
    allocate(ap(nx, ny, nz), aps(nx, ny, nz), aip(nx, ny, nz), aim(nx, ny, nz), ajp(nx, ny, nz), ajm(nx, ny, nz))
    allocate(akm(nx, ny, nz), akp(nx, ny, nz), gamk(nx, ny, nz))
    allocate(vu(nx, ny, nz), vv(nx, ny, nz), vw(nx, ny, nz))
    
    call Start
    
    u0 = u
    do while (t <= (tmax - 0.5*dt))
        t = t + dt      
        iter = 0
        
        do       
            ! Проверяем условие выхода
             if (iter >= numIterations) exit
        
            if (id == 0) then         
                call Info
            end if
        
            u1 = u
            
            ! Считаем коэффициенты уравнения и дискретных аналогов
            call BuildMatrix  
            
            ! Решаем СЛАУ
            !do iter = 1, 30
                call CD(u, aim, aip, ajm, ajp, akm, akp, ap, b)          
            !end do
            
            ! Производим обмен приграничными значениями
            call ExchangeBoundary     
               
            
            
            !call ConvergenceCondition
            
            iter = iter + 1
            alliter = alliter + 1                  
        end do
        
        if (id == 0) then         
            call Takeoff
        end if
        
        ! Промежуточный файл пишет только процесс 0
        
        
        u0 = u
    end do

    !call EquationCoeffs

    ! Собираем данные в единый массив
    call GatherData

    if (id == 0) then
        call CalcNorm
        call Output
    end if

end subroutine Solve


subroutine Info
integer nx2, ny2, nz2

    nx2 = nx / 2
    ny2 = ny / 2
    nz2 = nz / 2
    
    write(*, '(1I4P3E15.6)') iter, t, u(5, ny2, nz2), AnalyticSolution(x(5), y(ny2), z(nz2), t)     
    write(4, '(1I5P3E15.6)') alliter, u(5, ny2, nz2), u(nx2, ny2, nz2), u(nx - 5, ny2, nz2)
    
end subroutine Info



subroutine Output
real(8) maxabs, maxrel, minabs, minrel, ret

    ! Расчет граничных точек
    
    maxabs = 0.
    maxrel = 0.
    minabs = 100000.
    minrel = 100000.

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                ret = dabs(u(i, j, k) - AnalyticSolution(x(i), y(j), z(k), t))
                if (ret == 0.0e+0) continue
                
                maxabs = dmax1(maxabs, ret)
                minabs = dmin1(minabs, ret)
                
                ret = 100.0 * ret / norm
                maxrel = dmax1(maxrel, ret)
                minrel = dmin1(minrel, ret)

!                ret = dabs(gradu(i, j, k) - gradu0(i, j, k))
!                if (ret == 0.0e+0) continue
!                
!                maxabs = dmax1(maxabs, ret)
!                minabs = dmin1(minabs, ret)
!                
!                ret = 100.0 * ret / maxgrad
!                maxrel = dmax1(maxrel, ret)
!                minrel = dmin1(minrel, ret)
            end do
        end do
    end do

    !write(*, *)
    !write(1, *)

   print *, "Min, max abs err:", minabs, maxabs
   print *, "Min, max rel err:", minrel, maxrel
   print *, "Iterations: ", iter
   
   write(1, *) "Min, max abs err:", minabs, maxabs
   write(1, *) "Min, max rel err:", minrel, maxrel
   write(1, *) "Iterations: ", iter

end subroutine Output


subroutine Takeoff
integer nx2, ny2, nz2
character*256 fname
character*8 strtime

    nx2 = nx / 2
    ny2 = ny / 2
    nz2 = nz / 2
    
    write(3, '(1P3E15.6)') t, u(nx2, ny2, nz2), AnalyticSolution(x(nx2), y(ny2), z(nz2), t)

    784 format(F8.4)    
    write(strtime, 784) t
    fname = 'all_' // strtime // '.dat'
    
    open(unit = 5, file = fname, status = 'unknown')

    write(5, *) 'VARIABLES = "X" "Y" "Z" "U" "U0"'
    write(5, *) "ZONE I=", nx, "J=", ny, "K=", nz
    write(5, *) "SOLUTIONTIME=", t

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                write(5, '(1P5E15.6)') x(i), y(j), z(k), u(i, j, k), AnalyticSolution(x(i), y(j), z(k), t)
            end do
        end do
    end do

    endfile 5
    close(5)

end subroutine Takeoff


subroutine CalcNorm

    norm = 0.
    do i = 2, nx - 1
        do j = 2, ny - 1
            do k = 2, nz - 1
                norm = norm + ((AnalyticSolution(x(i), y(j), z(k), t))**2) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k))
            end do
        end do
    end do        
    norm = dsqrt(norm / ((xr - xl) * (yr - yl) * (zr - zl)))

end subroutine CalcNorm


subroutine Grad
real(8) gx, gy, gz

    allocate(gradu(nx, ny, nz), gradu0(nx, ny, nz))
    
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                u0(i, j, k) = AnalyticSolution(x(i), y(j), z(k), t)
            end do
        end do
    end do      
    
    maxgrad = 0.
    
    do i = 2, nx - 1
        do j = 2, ny - 1
            do k = 2, nz - 1
                gx = (u(i + 1, j, k) - u(i - 1, j, k)) / (x(i + 1) - x(i - 1))
                gy = (u(i, j + 1, k) - u(i, j - 1, k)) / (y(j + 1) - y(j - 1))
                gz = (u(i, j, k + 1) - u(i, j, k - 1)) / (z(k + 1) - z(k - 1))
                
                gradu(i, j, k) = dsqrt(gx * gx + gy * gy + gz * gz)
                
                gx = (u0(i + 1, j, k) - u0(i - 1, j, k)) / (x(i + 1) - x(i - 1))
                gy = (u0(i, j + 1, k) - u0(i, j - 1, k)) / (y(j + 1) - y(j - 1))
                gz = (u0(i, j, k + 1) - u0(i, j, k - 1)) / (z(k + 1) - z(k - 1))
                
                gradu0(i, j, k) = dsqrt(gx * gx + gy * gy + gz * gz)
                
                maxgrad = dmax1(dabs(gradu0(i, j, k)), maxgrad)
                
            end do
        end do
    end do            

end subroutine Grad


real(8) function Acoeff(P)
real(8) P

    if (P == 0) then
        Acoeff = 0.
    else
        Acoeff = dabs(P) / (exp(dabs(P)) - 1.)
    end if

end function Acoeff


! Расчет коэффициентов дискретного аналога
subroutine BuildMatrix
real(8) ap0, flow, pekle, diff

    call EquationCoeffs

    do i = 2, nx1
        do j = 2, ny1
            do k = 2, nz1
                aim(i, j, k) = gami(i, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k)) / (x(i) - x(i - 1))
                aip(i, j, k) = gami(i + 1, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k)) / (x(i + 1) - x(i))

                ajm(i, j, k) = gamj(i, j, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k)) / (y(j) - y(j - 1))
                ajp(i, j, k) = gamj(i, j + 1, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k)) / (y(j + 1) - y(j))
                
                akm(i, j, k) = gamk(i, j, k) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) / (z(k) - z(k - 1))
                akp(i, j, k) = gamk(i, j, k + 1) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) / (z(k + 1) - z(k))
                
!                ! aE
!                diff = gami(i + 1, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k)) / (x(i + 1) - x(i))
!                flow = vu(i + 1, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k))
!                pekle = flow / diff
!                aip(i, j, k) = diff * Acoeff(pekle) + dmax1(-flow, 0.)
!            
!                ! aW
!                diff = gami(i, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k)) / (x(i) - x(i - 1))
!                flow = vu(i, j, k) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k))
!                pekle = flow / diff
!                aim(i, j, k) = diff * Acoeff(pekle) + dmax1(flow, 0.)
!
!                ! aN
!                diff = gamj(i, j + 1, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k)) / (y(j + 1) - y(j))
!                flow = vv(i, j + 1, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k))
!                pekle = flow / diff
!                ajp(i, j, k) = diff * Acoeff(pekle) + dmax1(-flow, 0.)
!            
!                ! aS
!                diff = gamj(i, j, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k)) / (y(j) - y(j - 1))
!                flow = vv(i, j, k) * (xu(i + 1) - xu(i)) * (zu(k + 1) - zu(k))
!                pekle = flow / diff
!                ajm(i, j, k) = diff * Acoeff(pekle) + dmax1(flow, 0.)
!                
!                ! aH
!                diff = gamk(i, j, k + 1) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) / (z(k + 1) - z(k))
!                flow = vw(i, j, k + 1) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j))
!                pekle = flow / diff
!                akp(i, j, k) = diff * Acoeff(pekle) + dmax1(-flow, 0.)
!            
!                ! aB
!                diff = gamk(i, j, k) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) / (z(k) - z(k - 1))
!                flow = vw(i, j, k) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j))
!                pekle = flow / diff
!                akm(i, j, k) = diff * Acoeff(pekle) + dmax1(flow, 0.)

            end do
        end do
    end do

    do i = 2, nx1
        do j = 2, ny1
            do k = 2, nz1
                ap0 = rho(i, j, k) * (xu(i + 1) - xu(i)) * (yu(j + 1) - yu(j)) * (zu(k + 1) - zu(k)) / dt
                b(i, j, k) = con(i, j, k) + ap0 * u0(i, j, k)
                ap(i, j, k) = aim(i, j, k) + aip(i, j, k) + ajm(i, j, k) + ajp(i, j, k) + akm(i, j, k) + akp(i, j, k) + ap0 - aps(i, j, k)
            end do
        end do
    end do
    
    do i = 2, nx1
	    do j = 2, ny1
	        do k = 2, nz1
                ap(i, j, k) = ap(i, j, k) / relax                      
                b(i, j, k) = b(i, j, k) + (1. - relax) * ap(i, j, k) * u(i, j, k)
            end do
        end do
    end do

end subroutine BuildMatrix


! Условие сходимости
subroutine ConvergenceCondition  
    rmax = 0.
    do i = 2, nx1
        do j = 2, ny1
            do k = 2, nz1
                rmax = dmax1(rmax, dabs(1. - u(i, j, k) / u1(i, j, k)))
            end do
        end do
    end do
    
    call MPI_AllReduce(rmax, grmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (id == 0) print *, grmax


end subroutine ConvergenceCondition  


subroutine ExchangeBoundary
integer l
    if (numProcesses == 1) then 
        return
    end if
    
    l = ny * nz
    
    if (mod(id, 2) /= 0) then
        if (id /= numProcesses - 1) then 
            call MPI_Send(u(nx-1, :, :), l, MPI_DOUBLE_PRECISION, id + 1, 200, MPI_COMM_WORLD, ierr)
        end if
            
        if (id /= 0) then 
            call MPI_Send(u(2, :, :), l, MPI_DOUBLE_PRECISION, id - 1, 201, MPI_COMM_WORLD, ierr)
        end if
        
        if (id /= 0) then 
            call MPI_Recv(u(1, :, :), l, MPI_DOUBLE_PRECISION, id - 1, 202, MPI_COMM_WORLD, status, ierr)
        end if
            
        if (id /= numProcesses - 1) then 
            call MPI_Recv(u(nx, :, :), l, MPI_DOUBLE_PRECISION, id + 1, 203, MPI_COMM_WORLD, status, ierr)
        end if
        
    else
        if (id /= 0) then 
            call MPI_Recv(u(1, :, :), l, MPI_DOUBLE_PRECISION, id - 1, 200, MPI_COMM_WORLD, status, ierr)
        end if
            
        if (id /= numProcesses - 1) then 
            call MPI_Recv(u(nx, :, :), l, MPI_DOUBLE_PRECISION, id + 1, 201, MPI_COMM_WORLD, status, ierr)
        end if
        
        if (id /= numProcesses - 1) then 
            call MPI_Send(u(nx-1, :, :), l, MPI_DOUBLE_PRECISION, id + 1, 202, MPI_COMM_WORLD, ierr)
        end if
            
        if (id /= 0) then 
            call MPI_Send(u(2, :, :), l, MPI_DOUBLE_PRECISION, id - 1, 203, MPI_COMM_WORLD, ierr)
        end if
        
    end if

end subroutine ExchangeBoundary



! Обработка граничных условий второго и третьего рода, расчет коэффициентов уравнения
subroutine EquationCoeffs
real(8) V
    V = (xu(3) - xu(2)) * (yu(3) - yu(2)) * (zu(3) - zu(2))
    
    ! Значения во внутренних точках линеаризованного источникового члена и удельной теплоемкости
    do i = 2, nx1
        do j = 2, ny1
            do k = 2, nz1
                con(i, j, k) = -x(i) * dexp(-t) * V
                !con(i, j, k) = ((x(i)**2 - y(j) + z(k)) * dcos(t) + 2. * dsin(t) * (x(i) - 1.)) * V
                aps(i, j, k) = 0.
                rho(i, j, k) = 1.
            end do
        end do
    end do  
    
    ! Коэффициенты температуропроводности на гранях
    do i = 2, nx
        do j = 2, ny
            do k = 2, nz
                gami(i, j, k) = 1.
                gamj(i, j, k) = 1.
                gamk(i, j, k) = 1.
            end do
        end do
    end do   
    
    
    if (id == 0) then 
        ! Граничные условия первого рода, не зависящие от времени на левой границе
        do j = 1, ny
            do k = 1, nz
                u(1, j, k) = AnalyticSolution(x(1), y(j), z(k), t)
            end do
        end do
    end if

    if (id == numProcesses - 1) then
        ! Граничные условия первого рода, не зависящие от времени на правой границе
        do j = 1, ny
            do k = 1, nz
                u(nx, j, k) = AnalyticSolution(x(nx), y(j), z(k), t)
            end do
        end do
    end if

    do i = 1, nx
        do k = 1, nz
            u(i, 1, k) = AnalyticSolution(x(i), y(1), z(k), t)
            u(i, ny, k) = AnalyticSolution(x(i), y(ny), z(k), t)
        end do
    end do
    
    do i = 1, nx
        do j = 1, ny
            u(i, j, 1) = AnalyticSolution(x(i), y(j), z(1), t)
            u(i, j, nz) = AnalyticSolution(x(i), y(j), z(nz), t)
        end do
    end do
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
!    if (id == 0) then 
!        do j = 1, ny
!            do k = 1, nz
!                u(1, j, k) = u(2, j, k)
!                
!                gami(2, j, k) = 0.
!            end do
!        end do
!    end if
!
!    if (id == numProcesses - 1) then
!        do j = 1, ny
!            do k = 1, nz
!                u(nx, j, k) = u(nx1, j, k)
!                
!                gami(nx, j, k) = 0.
!            end do
!        end do
!    end if
!
!    do i = 1, nx
!        do k = 1, nz
!            u(i, 1, k) = u(i, 2, k)
!            u(i, ny, k) = u(i, ny1, k)
!            
!            gamj(i, 2, k) = 0.
!            gamj(i, ny, k) = 0.
!        end do
!    end do
!    
!    do i = 1, nx
!        do j = 1, ny
!            u(i, j, 1) = u(i, j, 2)
!            u(i, j, nz) = u(i, j, nz1)
!            
!            gamk(i, j, 2) = 0.
!            gamk(i, j, nz) = 0.
!        end do
!    end do
    
    

end subroutine EquationCoeffs


subroutine Start

    nx1 = nx - 1
    ny1 = ny - 1
    nz1 = nz - 1

    t = 0.
    tmax = 5.
    dt = (xu(3) - xu(2))**2
    alliter = 0
    relax = 0.8
    
    call CalcNorm
    
    
    ! Поле скорости
    do i = 2, nx
        do j = 2, ny
            do k = 2, nz
                vu(i, j, k) = 1. !xu(i)**2
                vv(i, j, k) = 1. !x(i) * yu(j)
                vw(i, j, k) = 1. !-3. * x(i) * zu(k)
            end do
        end do
    end do
    
    
    
    ! Начальное условие
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                u(i, j, k) = x(i) + y(j) * z(k)
            end do
        end do
    end do
    
end subroutine Start


! Функция аналитического решения
real(8) function AnalyticSolution(x, y, z, t)
real(8) x, y, z, t
    AnalyticSolution = x*dexp(-t) + y*z
end function AnalyticSolution


end