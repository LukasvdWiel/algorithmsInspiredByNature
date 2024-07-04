module parameters

implicit none

#ifdef debug
integer, parameter :: nCities = 10
integer, parameter :: nAnts = 2
#else
integer, parameter :: nCities = 100
integer, parameter :: nAnts = 10000
#endif


end module

program tsp
! Use a genetic algorithm to solves the travelling salesman problem.
! All cities are in the unit square.

use parameters

implicit none

double precision   :: pheromones(nCities, nCities)
double precision   :: distances(nCities, nCities)
double precision   :: x(nCities), y(nCities)
double precision   :: howNice(nCities)


integer            :: iCity, jCity, iAnt
integer            :: path(nCities)
integer            :: togo(nCities)

double precision   :: totalLength
double precision   :: addValue, max, r, selection
integer            :: togoEntry, weAreHere, destination
integer            :: from, to

integer, external  :: routeHasCity
double precision , external  :: length

double precision   :: progress(nAnts)

integer            :: shortestPath(nCities)
double precision   :: shortestLength


shortestLength = 999999d0

! initialize random cities
do iCity = 1, nCities
    call random_number(x(iCity))
    call random_number(y(iCity))
enddo


#ifdef debug
write(*,*) "Cities:"
do iCity=1,nCities
    write(*,*) iCity, x(iCity), y(iCity)
enddo
#endif

do iCity = 1, nCities
    do jCity = 1, nCities
        distances(iCity, jCity) = ((x(iCity) - x(jCity))**2 + (y(iCity) - y(jCity))**2  )
    enddo
enddo

#ifdef debug
write(*,*) "Distances:"
do iCity=1,nCities
    write(*,"(i3,10f6.2)") iCity, distances(iCity,:)
enddo
#endif

pheromones = 1d0



do iAnt = 1, nAnts
    ! The ant begint at city 1.
    ! It then chooses among the remaining cities based on density of pheromones en closeness of next city

#ifdef debug
    write(*,*) "-----------------------------------------------------"
    write(*,*) "----- start ant ", iAnt, "------------------"
    write(*,*) "-----------------------------------------------------"
#endif

    path = 0
    path(1) = 1

    ! create a path of the remaining cities.
    do iCity = 2, nCities

#ifdef debug
        write(*,*) "------ finding city", iCity, "-------"
#endif

        ! first select all the remaining cities that we have not yet visited
        ! (those do not yet appear in path)
        weAreHere = path(iCity-1)
        togo = 0
        togoEntry = 0
        do jCity = 1, nCities
            if (0 .eq. routeHasCity(path, jCity)) then
                ! yay, we found a city we have not visited yet
                togoEntry = togoEntry + 1
                togo(togoEntry) = jCity
            endif
        enddo

#ifdef debug
        write(*,*) "available cities: "
        write(*,"(10i3)") togo
#endif

        howNice = 0d0
        do jCity = 1, togoEntry
            ! more pheromones ->     nice!
            ! more distance   -> not nice!
            destination = togo(jCity)

            howNice(jCity+1) = pheromones(weAreHere, destination) / distances(weAreHere, destination)
        enddo

#ifdef debug
        write(*,*) "how nice to go to cities: "
        write(*,"(10f6.2)") howNice
#endif

        ! make a cumulative vector
        do jCity = 2, togoEntry+1
            howNice(jCity) = howNice(jCity-1) + howNice(jCity)
        enddo

#ifdef debug
        write(*,*) "cumulative how nice to go to cities: "
        write(*,"(10f6.2)") howNice
#endif

        ! select the next city based on this:
        call random_number(r)
        selection = r * howNice(togoEntry+1)

#ifdef debug
        write(*,*) "finding spot: ", selection
        write(*,*) "togoEntry: ", togoEntry
#endif


        do jCity = 1, togoEntry
            if ((selection .gt. howNice(jCity)) .and. &
                (selection .lt. howNice(jCity+1))) then
                ! found next destination
#ifdef debug
                write(*,"(a8,i3,a5,i3)") "between ", jCity, " and ", jCity+1
                write(*,*) "selected next destination: ", togo(jCity)
#endif
                path(iCity) = togo(jCity)
                cycle
            endif
        enddo

#ifdef debug
        write(*,*) "Path after this step: "
        write(*,"(10i3)") path
#endif

    enddo

#ifdef debug
    write(*,*) "The ant has walked the path: "
    write(*,"(10i3)") path
#endif


    ! write the first random walk to file
    if (iAnt .eq. 1) then
!        write(*,"(10i3)") path

        open(unit=41, file="initialPath.dat")
        do jCity = 1, nCities
            write(41,*) x(path(jCity)), y(path(jCity))
        enddo
        write(41,*) x(path(1)), y(path(1))
        close(41)
    endif




    ! update the pheromones based on the length of the path
    totalLength = length(x, y, path)
    if (totalLength .lt. shortestLength) then
        shortestLength = totalLength
        shortestPath = path
    endif

    progress(iAnt) = shortestLength
    addValue = 1d0 / totalLength

    

    do iCity = 2, nCities

        from = path(iCity-1)
        to   = path(iCity)

#ifdef debug
        write(*,*) "adapting pheronomes", from, to, addValue
#endif

        pheromones(from, to) = 0.999 * pheromones(from, to) + addValue
        pheromones(to, from) = 0.999 * pheromones(to, from) + addValue

    enddo


#ifdef debug
    write(*,*) "Updated pheromone matric:"
    do iCity=1,nCities
        write(*,"(i3,10f6.2)") iCity, pheromones(iCity,:)
    enddo
#endif


enddo


! all ants have marched, plot the best solution
open(unit=42, file="finalSolution.dat")
do iCity = 1, nCities
    write(42,*) x(shortestPath(iCity)), y(shortestPath(iCity))
enddo
write(42,*) x(path(1)), y(path(1))
close(42)

! print the progress of the solution
open(unit=44, file="convergence.dat")
do iAnt = 1, nAnts
    write(44,*) iAnt, progress(iAnt)
enddo
close(44)




end program

!----------------------------------------------------------

subroutine plotPheromones(iAnt, x, y, pheromones)

use parameters, only: nCities

implicit none

integer :: iAnt
double precision :: x(nCities), y(nCities)
double precision :: pheromones(nCities,nCities)

integer :: iCity, jCity

character(len=10) :: fileName





open(unit=43, file=fileName)
do iCity = 1, nCities-1
    do jCity = 1, nCities-1
        write(43,*) x(iCity), y(iCity)
        write(43,*) x(jCity), y(jCity)
        write(43,*) ">"
    enddo
enddo

close(43)


end subroutine

!----------------------------------------------------------

double precision function length(x, y, path)

use parameters, only: nCities

implicit none

double precision :: x(nCities), y(nCities)
integer          :: path(nCities)

double precision :: lentgh
integer          :: i

length = 0d0
do i=1,nCities-1
    length = length + sqrt( (x(path(i))-x(path(i+1)))**2 + &
                            (y(path(i))-y(path(i+1)))**2   )
enddo
! and at the final stretch top the home city
length = length + sqrt( (x(path(1))-x(path(nCities)))**2 + &
                        (y(path(1))-y(path(nCities)))**2   )

end function

!----------------------------------------------------------

integer function routeHasCity(route, city)

use parameters, only: nCities

implicit none

integer :: route(nCities)
integer :: city

integer :: i

routeHasCity = 0

do i=1, nCities
    if (route(i) .eq. city) then
        routeHasCity = city
        return
    endif
enddo

end function
