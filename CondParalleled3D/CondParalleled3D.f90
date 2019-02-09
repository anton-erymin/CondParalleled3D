use GlobalDefinitions
use GridDecomposition
use Solver

implicit none

real(8) startTime

    ! Инициализируем MPI и вводим данные
    call Initialization
    
    ! Строим сетку и производим декомпозицию
    xl = 0.
    xr = 6.
    
    yl = 0.
    yr = 1.
    
    zl = 0.
    zr = 1.
    
    call CreateGridAndDecompose
     
    ! Синхронизируем все процессы перед началом решения
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    if (id == 0) then
        startTime = MPI_Wtime()
    end if
    
    call Solve
    
    if (id == 0) then
        startTime = MPI_Wtime() - startTime
        print *, "Total time:", startTime, "s"
        write(1, *) "Total time:", startTime, "s"
        call AllFile
    end if
    
    
    ! Закрываем MPI
    call MPI_Finalize(ierr)

contains


subroutine Initialization

    ! Инициализируем MPI
    call MPI_Init(ierr)
    ! Получаем кол-во процессов в системе
    call MPI_Comm_Size(MPI_COMM_WORLD, numProcesses, ierr)
    ! Получаем ранг процесса
    call MPI_Comm_Rank(MPI_COMM_WORLD, id, ierr)   
    
    if (id == 0) then
        print *, "Initialization..."
        print *, "Number of processes:", numProcesses
        
        ! Вводим данные о разбиении
        write (*, *) "Divisions by X: "
        read (*, *) nx
        write (*, *) "Divisions by Y: "
        read (*, *) ny  
        write (*, *) "Divisions by Z: "
        read (*, *) nz    
        
        write (*, *) "Iterations: "
        read(*, *) numIterations
        
    end if
    
    ! Рассылаем остальным процессам введенные данные
    call MPI_Bcast(nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(numIterations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    if (id == 0) then 
        open(unit = 1, access = "append", file = 'stat.txt', status = 'unknown')
        open(unit = 3, file = 'delta.txt', status = 'unknown')
        open(unit = 4, file = 'info.txt', status = 'unknown')
    end if
    
end subroutine Initialization


subroutine AllFile
real(8) ret
integer i, j, k

    open(unit = 2, file = 'all.txt', status = 'unknown')

    write(2, *) 'VARIABLES = "X" "Y" "Z" "U" "U0" "Abs" "Rel"'
    write(2, *) "ZONE I=", nx, "J=", ny, "K=", nz

    do k = 1, nz
        do j = 1, ny
            do i = 1, n
                ret = AnalyticSolution(x(i), y(j), z(k), t)
                write(2, '(1P7E15.6)') x(i), y(j), z(k), u(i, j, k), ret, dabs(u(i, j, k) - ret), 100.0 * dabs((u(i, j, k) - ret) / norm)
                !write(2, '(1P7E15.6)') x(i), y(j), z(k), u(i, j, k), ret, dabs(gradu(i, j, k) - gradu0(i, j, k)), 100.0 * dabs((gradu(i, j, k) - gradu0(i, j, k)) / maxgrad)
                !write(2, '(1P6E15.6)') x(i), y(j), z(k), vu(i, j, k), vv(i, j, k), vw(i, j, k)
            end do
        end do
    end do

    endfile 2
    close(2)

end subroutine AllFile

end