module GridDecomposition
use GlobalDefinitions

implicit none
integer q
contains


subroutine GridDimension(n, left, right, nodes, edges)
integer n
real(8) left, right
real(8), allocatable :: nodes(:), edges(:)
integer i
real(8) h
    allocate(nodes(n), edges(n))

    h = (right - left) / dble(n - 2)

    edges(2) = left
    do i = 3, n
        edges(i) = edges(2) + (i - 2) * h
    end do

    nodes(1) = edges(2)
    do i = 2, n - 1
        nodes(i) = 0.5 * (edges(i) + edges(i + 1))
    end do
    nodes(n) = edges(n)
end subroutine GridDimension


subroutine CreateGridAndDecompose
    if (id == 0) then
        print *, "Creating grid..."
    end if

    ! Разбиение по y и z проводят все процессы
    call GridDimension(ny, yl, yr, y, yu)
    call GridDimension(nz, zl, zr, z, zu)
    
    if (numProcesses == 1) then
        call GridDimension(nx, xl, xr, x, xu)
        n = nx
    else 
        if (id == 0) then
            call SendDecompose
        else
            call RecvDecompose
        end if    
    end if    
end subroutine CreateGridAndDecompose


subroutine SendDecompose
integer pp, curEdge, gridSize, curPoint, i, k
real(8), allocatable, dimension(:) :: tx, txu, sx, sxu

    print *, "Decomposing..."
    
    ! Строим обычную сетку
    call GridDimension(nx, xl, xr, tx, txu)
    
    ! Далее распределяем сетки по процессам
    ! Кол-во расчетных точек на процесс
    !pp = (nx - 2) / numProcesses + 2
    pp = nx / numProcesses + 2
    
    allocate(sx(nx), sxu(nx))
    
    curPoint = 1
    curEdge = 2
    n = nx
    
    ! Цикл по всем процессам
    do i = 0, numProcesses - 1 
        if (i /= numProcesses - 1) then
            gridSize = pp
        
            do k = 1, pp
                sx(k) = tx(curPoint)
                curPoint = curPoint + 1    
            end do
        
            do k = 2, pp
                sxu(k) = txu(curEdge)
                curEdge = curEdge + 1    
            end do
    
            curPoint = curPoint - 2
            curEdge = curEdge - 1       
        else
            k = 1
            gridSize = 0
            do while (curPoint <= n)
                sx(k) = tx(curPoint)
                curPoint = curPoint + 1    
                k = k + 1
                gridSize = gridSize + 1
            end do   
            
            k = 2
            do while (curEdge <= n)
                sxu(k) = txu(curEdge)
                curEdge = curEdge + 1     
                k = k + 1
            end do     
        end if
        
        ! Каждому процессу, кроме самого себя, посылаем массив расчетных точек и граней
        if (i /= 0) then
            call MPI_Send(gridSize, 1, MPI_INTEGER, i, 100, MPI_COMM_WORLD, ierr)
            call MPI_Send(sx,  gridSize, MPI_DOUBLE_PRECISION, i, 101, MPI_COMM_WORLD, ierr)
            call MPI_Send(sxu, gridSize, MPI_DOUBLE_PRECISION, i, 102, MPI_COMM_WORLD, ierr)
        else
            ! Корневой процесс копирует сетку в расчетные массивы x, xu
            
            nx = gridSize   
            allocate(x(nx), xu(nx))  
            x = sx
            xu = sxu
            
            print *, id, ': Сетку получил. Размер: ', nx
!            print *, id, (xu(q), q = 2, nx)
!            print *, id, '---------------------------'
!            print *, id, (x(q), q = 1, nx)
!            print *, '==========================='
        end if
        
    end do

end subroutine SendDecompose


subroutine RecvDecompose
    ! Прием расчетных точек
    call MPI_Recv(nx, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, status, ierr)
    ! Выделяем память
    allocate(x(nx), xu(nx))
    call MPI_Recv(x, nx, MPI_DOUBLE_PRECISION, 0, 101, MPI_COMM_WORLD, status, ierr)
    ! Прием граней
    call MPI_Recv(xu, nx, MPI_DOUBLE_PRECISION, 0, 102, MPI_COMM_WORLD, status, ierr)
            
            print *, id, ': Сетку получил. Размер: ', nx
!            print *, id, (xu(q), q = 2, nx)
!            print *, id, '---------------------------'
!            print *, id, (x(q), q = 1, nx)
!            print *, '==========================='
end subroutine RecvDecompose


subroutine GatherData
integer i, j, k, curPoint, num, pr, w
real(8), allocatable :: tx(:)
    if (numProcesses == 1) then 
        return
    end if
    
    ! Формируем результат
    ! Некорневые процессы отправляют частное решение корневому
    if (id /= 0) then
        do k = 1, nz
            do j = 1, ny
                call MPI_Send(u(:, j, k), nx, MPI_DOUBLE_PRECISION, 0, 10 * id, MPI_COMM_WORLD, ierr)    
            end do
        end do
    
    ! Корневой процесс принимает решения от остальных процессов
    else
        u1 = u
        deallocate(u, u0)
        allocate(u(n, ny, nz), u0(n, ny, nz))
        allocate(tx(n))
    
        do j = 1, ny
            do i = 1, nx
                do k = 1, nz
                    u(i, j, k) = u1(i, j, k)
                end do
             end do
        end do
    
        curPoint = nx - 1
        
        do pr = 1, numProcesses - 1
        
            do k = 1, nz
                do j = 1, ny
                    call MPI_Recv(tx, n, MPI_DOUBLE_PRECISION, pr, 10 * pr, MPI_COMM_WORLD, status, ierr) 
                    call MPI_Get_Count(status, MPI_DOUBLE_PRECISION, num, ierr) 
                
                    do w = 0, num - 1
                        u(curPoint + w, j, k) = tx(w + 1)
                    end do
                end do
            end do        
            
            curPoint = curPoint + num - 2
                
        end do
        
        deallocate(x, xu)
        call GridDimension(n, xl, xr, x, xu)
        nx = n
        
    end if
end subroutine GatherData


end module