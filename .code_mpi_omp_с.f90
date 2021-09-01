module types
    implicit none
    
    type vector
        double precision :: x, y
    end type
    
    type flow
        type(vector) :: speed
        double precision :: energy, density, pressure
    end type
    
    type direction
        double precision :: up, down, left, right
    end type
contains
end module types

module init
    use types
    implicit none
    
    double precision, parameter :: GasConstant = 8.3145, Mach = 2.0, AdiabaticIndex = 1.4, MolarMass = 28.98
    double precision, parameter :: deltaX = 0.1, deltaY = 0.1, deltaT = 0.01, width = 10., height = 10., time = 300.0
    integer, parameter :: countFrames = 100, countThreads = 4
    logical, parameter :: isMiddle = .true.
    
contains
    function initial(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: initial
        initial = table(i, j)
        
        initial%speed%x = 0.1e-15
        initial%speed%y = 0.1e-15
        initial%density = 1.0
        initial%pressure = 1.0
        initial%energy = initial%pressure / ((AdiabaticIndex - 1) * initial%density)
    end function initial
    
    function center(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: center
        center = table(i, j)
    end function center
    
    function left(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: left
        left = table(i, j)
        
        left%speed%x = -table(i + 1, j)%speed%x
        left%speed%y = table(i + 1, j)%speed%y
        left%density = table(i + 1, j)%density
        left%energy = table(i + 1, j)%energy
    end function left
    
    function right(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: right
        right = table(i, j)
        
        if (j < (int(0.5 * height / deltaY)) + 1) then
            right%speed%x = -table(i - 1,j)%speed%x
        else
            right%speed%x = table(i - 1,j)%speed%x
        end if
        right%speed%y = table(i - 1,j)%speed%y
        right%density = table(i - 1,j)%density
        right%energy = table(i - 1,j)%energy
    end function right
    
    function down(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: down
        down = table(i, j)
        
        down%speed%x = table(i, j + 1)%speed%x
        down%speed%y = -table(i, j + 1)%speed%y
        down%density = table(i, j + 1)%density
        down%energy = table(i, j + 1)%energy
        
        if (i == (int(0.5 * width / deltaX)) + 1 .and. j == 1) then
            down%energy = 30.
        end if
    end function down
        
    function up(table, i, j)
        integer :: i, j
        type(flow), dimension(:,:) :: table
        intent(in) :: i, j, table
        type(flow) :: up
        up = table(i, j)
        
        up%speed%x = table(i, j - 1)%speed%x
        up%speed%y = -table(i, j - 1)%speed%y
        up%density = table(i, j - 1)%density
        up%energy = table(i, j - 1)%energy
    end function up
end module

module transmission
    include 'mpif.h'
    
    integer :: ierr, numberThreads, currentThread, status(MPI_STATUS_SIZE)
    integer, dimension(1) :: rankD
    integer, dimension(2) :: rankDD
    integer, dimension(:), allocatable :: initIndexes, countIndexes
    double precision, dimension(:), allocatable :: buffer
    
contains
    subroutine transmission_init(x, y)
        integer :: i, x, y
        
        call MPI_INIT(ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, numberThreads, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, currentThread, ierr)
    
        allocate(initIndexes(numberThreads))
        allocate(countIndexes(numberThreads))
        allocate(buffer((x + 2) * (y + 2)))
        
        do i = 0, numberThreads
            initIndexes(i + 1) = 1 + i * (y + 2) / numberThreads
            countIndexes(i + 1) = (1 + i) * (y + 2) / numberThreads - initIndexes(i + 1) + 1
        end do
        
        rankD = (/(x + 2) * (countIndexes(currentThread + 1))/)
        rankDD = (/x + 2, y + 2/)
    end subroutine
    
    subroutine gather(input, destination, output)
        integer :: destination
        double precision, dimension(rankDD(1), countIndexes(currentThread + 1)) :: input
        double precision, dimension(rankDD(1), rankDD(2)) :: output
        
        call MPI_GATHERV(reshape(input, rankD), countIndexes(currentThread + 1) * rankDD(1), MPI_DOUBLE_PRECISION, buffer, countIndexes * rankDD(1), (initIndexes - 1) * rankDD(1), MPI_DOUBLE_PRECISION, destination, MPI_COMM_WORLD, ierr)
        if (currentThread == destination) then
            output = reshape(buffer, rankDD)
        end if
    end subroutine
    
    subroutine sendRecv(input, destinationOffset, output)
        integer :: destinationOffset
        double precision, dimension(rankDD(1)) :: input, output
        
        if (destinationOffset <= 0 .and. currentThread > -destinationOffset - 1 .or. destinationOffset > 0 .and. currentThread < numberThreads - destinationOffset) then
            call MPI_SEND(input, rankDD(1), MPI_DOUBLE_PRECISION, currentThread + destinationOffset, 0, MPI_COMM_WORLD, ierr)
        end if
        if (destinationOffset > 0 .and. currentThread > destinationOffset - 1 .or. destinationOffset <= 0 .and. currentThread < numberThreads + destinationOffset) then
            call MPI_RECV(output, rankDD(1), MPI_DOUBLE_PRECISION, currentThread - destinationOffset, 0, MPI_COMM_WORLD, status, ierr)
        end if
    end subroutine
end module transmission

program Main
    use init
    use transmission
    use omp_lib
    
    double precision :: start, finish
    integer :: x, y, t
    
    x = int(width/deltaX)
    y = int(height/deltaY)
    t = int(time/deltaT)
    
    call omp_set_num_threads(countThreads)
    
    call transmission_init(x, y)
    if (currentThread == 0) then
        print '("MPI_N =", 1x, i3)', numberThreads
        print '("OMP_N =", 1x, i3)', omp_get_max_threads()
    end if
    start = MPI_WTIME()
    call largeParticlesMethod(x, y, t)
    finish = MPI_WTIME() - start
    call MPI_REDUCE(finish, start, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (currentThread == 0) then
        print '("t =", 1x, f20.10, 1x, "s")', start
    end if
    
    call MPI_FINALIZE(ierr)
end program

subroutine largeParticlesMethod(x, y, t)
    use init
    use transmission
    use omp_lib
    
    integer :: x, y, t, i, j, k, percent
    double precision :: tmpX, tmpY
    type(direction) :: mass, route
    type(flow), dimension(x + 2, y + 2) :: table
    type(vector), dimension(x + 2, y + 2) :: tmpS
    double precision, dimension(x + 2, y + 2) :: tmpE, tmpD
    
    if (currentThread == 0) then
        print '("Starting large particle method:")'
    end if
    
    !Начальные условия
    do i = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
        if (i /= 1 .and. i /= size(table, 2)) then
            !$omp parallel do
            do j = 2, size(table, 1) - 1
                table(j, i) = initial(table, j, i)
            end do
            !$omp end parallel do
        end if
    end do
    
    !Передача записанных начальных условий в главный узел для вывода в файл
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%speed%x, 0, table%speed%x)
    call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%speed%y, 0, table%speed%y)
    call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%density, 0, table%density)
    call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%pressure, 0, table%pressure)
    call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%energy, 0, table%energy)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    !Вывод в файл в главном узле
    if (currentThread == 0 .and. isMiddle) then
        call appendTableToFile(table, size(table, 1), size(table, 2), "gnuPlot", .false., .true.)
    end if

    do k = 1, t
        !Граничные условия + Расчет давления v5.0
        if (currentThread == 0 .or. currentThread == numberThreads - 1) then
            !$omp parallel do
            do j = 1, size(table, 1)
                if (currentThread == 0) then
                    table(j, 1) = down(table, j, 1)
                end if
                if (currentThread == numberThreads - 1) then
                    table(j, size(table, 2)) = up(table, j, size(table, 2))
                end if
            end do
            !$omp end parallel do
        end if
        
        !$omp parallel do
        do i = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
            table(1, i) = left(table, 1, i)
            table(size(table, 1), i) = right(table, size(table, 1), i)
        end do
        !$omp end parallel do
        
        do i = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
            if (i /= 1 .and. i /= size(table, 2)) then
                !$omp parallel do
                do j = 2, size(table, 1) - 1
                    table(j, i) = center(table, j, i)
                end do
                !$omp end parallel do
            end if
        end do
        
        do i = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
            !$omp parallel do
            do j = 1, size(table, 1)
                table(j, i)%pressure = table(j, i)%density  * (AdiabaticIndex - 1) * (table(j, i)%energy - &
                     0.5 * (table(j, i)%speed%x ** 2 + table(j, i)%speed%y ** 2))
            end do
            !$omp end parallel do
        end do
        
        !Передача крайних столбцов необходимых параметров расчетной матрицы соседним вычислительным узлам
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call sendRecv(table(:, initIndexes(currentThread + 1))%pressure, -1, table(:, initIndexes(currentThread + 2))%pressure)
        call sendRecv(table(:, initIndexes(currentThread + 2) - 1)%pressure, 1, table(:, initIndexes(currentThread + 1) - 1)%pressure)
        call sendRecv(table(:, initIndexes(currentThread + 1))%speed%y, -1, table(:, initIndexes(currentThread + 2))%speed%y)
        call sendRecv(table(:, initIndexes(currentThread + 2) - 1)%speed%y, 1, table(:, initIndexes(currentThread + 1) - 1)%speed%y)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        
        !Эйлеров этап
        tmpX = deltaT/2
        tmpY = tmpX/deltaY
        tmpX = tmpX/deltaX
        tmpS = table%speed
        tmpE = table%energy
        do j = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
            if (j /= 1 .and. j /= size(table, 2)) then
                !$omp parallel do
                do i = 2, size(table, 1) - 1
                    tmpS(i, j)%x = table(i, j)%speed%x - &
                        tmpX * (table(i + 1, j)%pressure - table(i - 1, j)%pressure) / table(i, j)%density
                    tmpS(i, j)%y = table(i, j)%speed%y - &
                        tmpY * (table(i, j + 1)%pressure - table(i, j - 1)%pressure) / table(i, j)%density
                    tmpE(i, j) = table(i, j)%energy - (&
                        tmpX * ((table(i, j)%pressure + table(i + 1, j)%pressure) * &
                            (table(i, j)%speed%x + table(i + 1, j)%speed%x) - &
                            (table(i, j)%pressure + table(i - 1, j)%pressure) * &
                            (table(i, j)%speed%x + table(i - 1, j)%speed%x)) + &
                        tmpY * ((table(i, j)%pressure + table(i, j + 1)%pressure) * &
                            (table(i, j)%speed%y + table(i, j + 1)%speed%y) - &
                            (table(i, j)%pressure + table(i, j - 1)%pressure) * &
                            (table(i, j)%speed%y + table(i, j - 1)%speed%y))) / (2 * table(i, j)%density)
                end do
                !$omp end parallel do
            end if
        end do
        
        !Передача крайних столбцов необходимых параметров расчетной матрицы соседним вычислительным узлам
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call sendRecv(tmpS(:, initIndexes(currentThread + 1))%x, -1, tmpS(:, initIndexes(currentThread + 2))%x)
        call sendRecv(tmpS(:, initIndexes(currentThread + 2) - 1)%x, 1, tmpS(:, initIndexes(currentThread + 1) - 1)%x)
        call sendRecv(tmpS(:, initIndexes(currentThread + 1))%y, -1, tmpS(:, initIndexes(currentThread + 2))%y)
        call sendRecv(tmpS(:, initIndexes(currentThread + 2) - 1)%y, 1, tmpS(:, initIndexes(currentThread + 1) - 1)%y)
        call sendRecv(table(:, initIndexes(currentThread + 1))%density, -1, table(:, initIndexes(currentThread + 2))%density)
        call sendRecv(table(:, initIndexes(currentThread + 2) - 1)%density, 1, table(:, initIndexes(currentThread + 1) - 1)%density)
        call sendRecv(tmpE(:, initIndexes(currentThread + 1)), -1, tmpE(:, initIndexes(currentThread + 2)))
        call sendRecv(tmpE(:, initIndexes(currentThread + 2) - 1), 1, tmpE(:, initIndexes(currentThread + 1) - 1))
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        
        do j = initIndexes(currentThread + 1), initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1
            if (j /= 1 .and. j /= size(table, 2)) then
                !$omp parallel do private(tmpX, tmpY, mass, route)
                do i = 2, size(table, 1) - 1
                    !Лангранжев этап
                    tmpY = deltaT / 2
                    tmpX = tmpY * deltaY * (tmpS(i, j)%x + tmpS(i + 1, j)%x)
                    tmpY = tmpY * deltaX * (tmpS(i, j)%y + tmpS(i, j + 1)%y)
                    if (tmpX > 0) then
                        mass%right = abs(table(i, j)%density * tmpX)
                        route%right = 0.
                    else
                        mass%right = abs(table(i + 1, j)%density * tmpX)
                        route%right = 1.
                    end if
                    if (tmpY > 0) then
                        mass%up = abs(table(i, j)%density * tmpY)
                        route%up = 0.
                    else
                        mass%up = abs(table(i, j + 1)%density * tmpY)
                        route%up = 1.
                    end if
                    
                    tmpY = deltaT / 2
                    tmpX = tmpY * deltaY * (tmpS(i, j)%x + tmpS(i - 1, j)%x)
                    tmpY = tmpY * deltaX * (tmpS(i, j)%y + tmpS(i, j - 1)%y)
                    if (tmpX > 0) then
                        mass%left = abs(table(i - 1, j)%density * tmpX)
                        route%left = 1.
                    else
                        mass%left = abs(table(i, j)%density * tmpX)
                        route%left = 0.
                    end if
                    if (tmpY > 0) then
                        mass%down = abs(table(i, j - 1)%density * tmpY)
                        route%down = 1.
                    else
                        mass%down = abs(table(i, j)%density * tmpY)
                        route%down = 0.
                    end if
                    
                    !Заключительный этап
                    tmpX = deltaX * deltaY
                    tmpY = 1 / tmpX
                    tmpX = tmpX * table(i, j)%density
                    tmpD(i, j) = table(i, j)%density + 2 * tmpY * ((route%left - 0.5) * mass%left + &
                        (route%down - 0.5) * mass%down + (route%right - 0.5) * mass%right + (route%up - 0.5) * mass%up)
                    table(i, j)%speed%x = tmpY / tmpD(i, j) * &
                        (tmpS(i, j)%x * (tmpX - &
                        (1 - route%left) * mass%left - (1 - route%down) * mass%down - &
                        (1 - route%right) * mass%right - (1 - route%up) * mass%up) + &
                        tmpS(i - 1, j)%x * route%left * mass%left + tmpS(i, j - 1)%x * route%down * mass%down + &
                        tmpS(i + 1, j)%x * route%right * mass%right + tmpS(i, j + 1)%x * route%up * mass%up)
                    table(i, j)%speed%y = tmpY / tmpD(i, j) * &
                        (tmpS(i, j)%y * (tmpX - &
                        (1 - route%left) * mass%left - (1 - route%down) * mass%down - &
                        (1 - route%right) * mass%right - (1 - route%up) * mass%up) + &
                        tmpS(i - 1, j)%y * route%left * mass%left + tmpS(i, j - 1)%y * route%down * mass%down + &
                        tmpS(i + 1, j)%y * route%right * mass%right + tmpS(i, j + 1)%y * route%up * mass%up)
                    table(i, j)%energy = tmpY / tmpD(i, j) * &
                        (tmpE(i, j) * (tmpX - &
                        (1 - route%left) * mass%left - (1 - route%down) * mass%down - &
                        (1 - route%right) * mass%right - (1 - route%up) * mass%up) + &
                        tmpE(i - 1, j) * route%left * mass%left + tmpE(i, j - 1) * route%down * mass%down + &
                        tmpE(i + 1, j) * route%right * mass%right + tmpE(i, j + 1) * route%up * mass%up)
                end do
                !$omp end parallel do
            end if
        end do
        table%density = tmpD
        if (currentThread == 0 .and. 100 * k / t /= percent) then
            percent = 100 * k / t
            print '(3x,i3,"%")', percent
        end if

        if (.not. isMiddle .and. k == t .or. isMiddle .and. (t < countFrames .or. mod(k, t / countFrames) == 0)) then
            !Передача расчитанных промежуточных значений в главный узел для вывода в файл
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%speed%x, 0, table%speed%x)
            call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%speed%y, 0, table%speed%y)
            call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%density, 0, table%density)
            call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%pressure, 0, table%pressure)
            call gather(table(:, initIndexes(currentThread + 1) : initIndexes(currentThread + 1) + countIndexes(currentThread + 1) - 1 : 1)%energy, 0, table%energy)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            
            !Вывод промежуточных и конечных значений в файл в главном узле
            if (currentThread == 0) then
                call appendTableToFile(table, size(table, 1), size(table, 2), "gnuPlot", (k+t/countFrames) > t, .not. isMiddle)
            end if
        end if
    end do
    if (currentThread == 0) then
        print '("End calculating")'
    end if
end subroutine

subroutine appendTableToFile(table, x, y, path, isLast, isFirstSequence)
    use init
    implicit none
    
    logical, intent(in) :: isLast, isFirstSequence
    double precision :: tmp
    integer :: i, j, x, y, info
    character(len=*) :: path
    type(flow), dimension(x, y) :: table
    
    if (isFirstSequence) then
        open(1, file = trim(path)//"_values.dat", status = 'replace')
        open(2, file = trim(path)//"_consts.dat", status = 'replace')
    else
        open(1, file = trim(path)//"_values.dat", iostat=info, status = 'old', position = 'append')
        if (info == 2) then
            open(1, file = trim(path)//"_values.dat", iostat=info, status = 'new')
        end if
    end if
    if (isFirstSequence) then
        write (2, "(f20.10, 1x, f20.10, 1x, f20.10)") deltaX, deltaY, deltaT
    end if
    do i = 2, size(table, 1) - 1
        do j = 2, size(table, 2) - 1
            tmp = sqrt(table(i, j)%speed%x ** 2 + table(i, j)%speed%y ** 2)
            write(1, "(f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10, 1x, f20.10)") &
                (i - 1) * deltaX, (j - 1) * deltaY, table(i, j)%density, table(i, j)%energy, table(i, j)%pressure, &
                deltaX * table(i, j)%speed%x / tmp, deltaY * table(i, j)%speed%y / tmp, tmp, &
                (table(i, j)%pressure * MolarMass) / (table(i, j)%density * GasConstant)
        end do    
    end do
    if (.not. isLast) then
        write(1, "(a)")
        write(1, "(a)")
    end if
    close(1)
    close(2)
end subroutine