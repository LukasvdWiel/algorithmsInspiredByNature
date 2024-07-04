module parameters

implicit none

#ifdef debug
integer, parameter :: nCities = 10
integer, parameter :: populationSize = 10
integer, parameter :: nDeath = 5
integer, parameter :: nGenerations = 2
double precision, parameter :: mutationChance = 0.5
#else
integer, parameter :: nCities = 100
integer, parameter :: populationSize = 1000
integer, parameter :: nDeath = 500
integer, parameter :: nGenerations = 1000
double precision, parameter :: mutationChance = 0.02
#endif

integer, parameter :: nSurvivors = populationSize - nDeath

end module

program tsp
! Use a genetic algorithm to solves the travelling salesman problem.
! All cities are in the unit square.

use parameters
implicit none
double precision   :: x(nCities), y(nCities)
integer            :: paths(populationSize, nCities)
integer            :: iGeneration, iCity, iPop, i
double precision   :: lengths(populationSize)
integer            :: momID, dadID, childID
double precision   :: r

double precision, external :: length

integer            :: sortedIndex(populationSize)

integer            :: shortestPath(nCities)
double precision   :: shortestLength


shortestLength = 99999d0

! initialize random cities
do iCity = 1, nCities
    call random_number(x(iCity))
    call random_number(y(iCity))
enddo

#ifdef debug
write(*,*) "Cities:"
do i=1,nCities
    write(*,*) i, x(i), y(i)
enddo
#endif


! initialize random routes
do iPop = 1, populationSize
    call initialisePath(paths(iPop,:))
enddo

#ifdef debug
write(*,*) "initial routes:"
do i=1,populationSize
    write(*,"(10i3)") paths(i, :)
enddo
#endif


! print the data of an initial distribution for plotting (path 1)
open(unit=51, file="initialPath.dat")
do iCity = 1, nCities
    write(51,*) x(paths(1,iCity)), y(paths(1,iCity))
enddo
! wrap back to city 1
write(51,*) x(paths(1,1)), y(paths(1,1))
close(51)

open(unit=61, file="convergence.dat")

! run the evolution
do iGeneration = 1, nGenerations

#ifdef debug
    write(*,*) "-----------------------------------------------------"
    write(*,*) "----- start generation ", iGeneration, "------------------"
    write(*,*) "-----------------------------------------------------"
#endif


    ! determine the fitness of all the individuals
    do iPop = 1, populationSize
        lengths(iPop) = length(x, y, paths(iPop,:))
    enddo

#ifdef debug
do i=1,populationSize
    write(*,*) "lengths: ", i, lengths(i)
enddo
#endif

    ! determine the indices of the specimens that are dying,
    ! because they have the lowest fitness, and the indices
    ! of the survivors

    call sortLengths(lengths, sortedIndex)    

    if (lengths(1) .lt. shortestLength) then
        shortestLength = lengths(1)
        shortestPath = paths(sortedIndex(1),:)
    endif
    write(61,*) log10(dble(iGeneration)), shortestLength



    ! fill the vacated spots with offspring from the survivors

#ifdef debug
    write(*,*) "to be kept: "
    write(*,"(5i3)") sortedIndex(1:5)

    write(*,*) "to be replaced"
    write(*,"(5i3)") sortedIndex(6:10)
#endif


    do i = 1, nDeath

#ifdef debug
    write(*,*) "------------------- reproduction ", i, "------------"
#endif

        ! determine position of new child
        childID = sortedIndex(nSurvivors + i)

        ! determine positions of parents
        call random_number(r)
        momID = int(ceiling(r*nSurvivors-1))+1
        call random_number(r)
        dadID = int(ceiling(r*nSurvivors-1))+1
        dadID = mod(dadID + momID, nSurvivors)+1

        momID = sortedIndex(momID)
        dadID = sortedIndex(dadID)

#ifdef debug
        write(*,*) "replace path ", childID
        write(*,*) "use mom ", momID
        write(*,*) "use dad ", dadID  
        write(*,"(a5,10i3)") "mom: ", paths(momID,:)
        write(*,"(a5,10i3)") "dad: ", paths(dadID,:)
#endif

        ! and create the actual child
        call createOffspring(paths(momID,:), paths(dadID,:), paths(childID,:))

    enddo


    ! apply mutations by reversing a section of the path
    do iPop = 1, populationSize
        call random_number(r)
        if (r.lt.mutationChance) then
#ifdef debug
            write(*,*) "---- mutation of path", iPop, "------------"
#endif
            call reverseRandomSection(paths(iPop,:))
        endif
    enddo

enddo

close(61)



! order the fitness of all the individuals
do iPop = 1, populationSize
    lengths(iPop) = length(x, y, paths(iPop,:))
enddo
call sortLengths(lengths, sortedIndex)


! and print the data of the best one
open(unit=51, file="finalSolution.dat")
do iCity = 1, nCities

    write(51,*) x(shortestPath(iCity)), y(shortestPath(iCity))

!    write(51,*) x(paths(sortedIndex(1),iCity)), y(paths(sortedIndex(1),iCity))
enddo
! wrap back to city 1
write(51,*) x(shortestPath(1)), y(shortestPath(1))

close(51)




end program


!----------------------------------------------------------

subroutine createOffspring(mom, dad, child)

use parameters, only: nCities

implicit none

integer :: mom(nCities)
integer :: dad(nCities)
integer :: child(nCities)

integer :: checkMe, i, fillChildPos
integer :: nFromMom, momStart, momEnd
logical :: dadNumsUsed(nCities)

integer, external :: routeHasCity

double precision :: r

! we cannot randomly split and glue together bits of genes,
! because we must ensure that every city occurs exactly and only once.

! to this, we take a random bit from mom, and fill the
! missing bits up following the sequence from dad.

! keep in mind, though, the the journey is circular,
! so the bit inherited from mom can wrap around the integer edge

child = 0

nFromMom = int(floor(nCities * 0.5))

call random_number(r)
momStart = int(ceiling(r * nCities))
momEnd = momStart + nFromMom - 1
momEnd = mod(momEnd, nCities) ! possibly wrap

#ifdef debug
    write(*,*) "child gets mom bits", momStart, "to", momEnd
#endif


if (momEnd .gt. momStart) then
    ! there was no wrapping
    child(momStart:momEnd) = mom(momStart:momEnd)
else
    ! there is wrapping, copy mom's bit in two takes
    child(momStart:nCities) = mom(momStart:nCities)
    child(1:momEnd)         = mom(1:momEnd)
endif

#ifdef debug
write(*,"(10i3)") child
#endif

! it is possible that momEnd is the end of the array.
! prevent the new position from running out:
if (momEnd .eq. nCities) then
    fillChildPos = 1
else
    fillChildPos = momEnd + 1
endif

! Step through the father and copy unknown entries.
! TODO: speed this up with a looking table,
do i = 1,nCities
    checkMe = dad(i)
    if (0 .eq. routeHasCity(child, checkMe)) then
        ! yay, we have a city from dad the mother did not give.
        ! Add it!
#ifdef debug
        write(*,*) "dad pos", i, "has", checkMe, "which is unknown, insert in child pos", fillChildPos
#endif        
        child(fillChildPos) = checkMe
#ifdef debug
        write(*,"(10i3)") child
#endif


        if (fillChildPos .eq. nCities) then
            fillChildPos = 1
        else
            fillChildPos = fillChildPos + 1
        endif

    else
        ! This city is already give by the mother. Do nothing.
#ifdef debug
        write(*,*) "dad pos", i, "has", checkMe, "which is already in child. Skip it."
#endif

    endif
enddo

#ifdef debug
        write(*,*) "finished child"
        write(*,"(10i3)") child
#endif




end subroutine

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

subroutine initialisePath(path)

use parameters, only: nCities

implicit none

integer          :: path(nCities)
integer          :: n, m, temp, i
double precision :: r, s

do i=1,nCities
    path(i) = i
enddo

do i=1,nCities**2
    call switchTwoRandomPoints(path)
enddo

end subroutine

!----------------------------------------------------------

subroutine reverseRandomSection(path)

use parameters, only: nCities

implicit none

integer, parameter :: nMin = 2
integer, parameter :: nMax = int(floor(nCities * 0.5))

integer :: path(nCities)
integer :: section(nMax)
integer :: i
integer :: pos

integer :: reverseSize
integer :: sectionStart, sectionEnd
double precision :: r

! find start and end point of the section to reverse

#ifdef debug
write(*,*) "Reversing a section of path: "
write(*,"(10i3)") path
#endif

! determine the size. At least 2, at most half
call random_number(r)
reverseSize = nMin + int(ceiling(r * (nMax - nMin)))


! determine the boundaries of the section
call random_number(r)
sectionStart = int(ceiling(r * nCities))
sectionEnd = sectionStart + reverseSize - 1
sectionEnd = mod(sectionEnd, nCities) ! possibly wrap

#ifdef debug
write(*,*) "Reversing a section of size ", reverseSize, "from", sectionStart, "to", sectionEnd
#endif

! put the section in a buffer
do i=1,reverseSize
    pos = mod(sectionStart + i -2, nCities)+1
#ifdef debug
    write(*,*) "get path entry", pos, "in buffer pos", i
#endif
    section(i) = path(pos)
enddo

#ifdef debug
write(*,*) "Section buffered: ", section
#endif

! put the buffer back in the section in reverse order
do i=1,reverseSize
    pos = mod(sectionStart + i-2, nCities)+1
#ifdef debug
    write(*,*) "put section entry", reverseSize - i + 1, "in pos: ", pos
#endif
    path(pos) = section(reverseSize - i + 1)
enddo

#ifdef debug
write(*,*) "Reversal of section complere. New path:"
write(*,"(10i3)") path
#endif

end subroutine


!----------------------------------------------------------


subroutine switchTwoRandomPoints(path)
! used to randomize the initial solutions

use parameters, only: nCities

implicit none

integer          :: path(nCities)
integer          :: n, m, temp
double precision :: r, s

call random_number(r)
n = int(ceiling(r*nCities-1))+1
call random_number(s)
m = int(ceiling(s*nCities-1))+1
m = mod(m+n, nCities)+1

temp = path(n)
path(n) = path(m)
path(m) = temp

end subroutine

!----------------------------------------------------------

subroutine sortLengths(values, sortedIndex)

use parameters, only: populationSize

implicit none

double precision :: values(populationSize)
integer          :: sortedIndex(populationSize)

integer          :: i, j
integer          :: tempInt
double precision :: tempDoub

! use a shameless bubble sort

do i = 1, populationSize
    sortedIndex(i) = i
enddo

do i=1, populationSize-1
    do j=1, populationSize-1

        if (values(j) .gt. values(j+1)) then 
            tempDoub = values(j)
            values(j) = values(j+1)
            values(j+1) = tempDoub

            tempInt = sortedIndex(j)
            sortedIndex(j) = sortedIndex(j+1)
            sortedIndex(j+1) = tempInt
        endif

    enddo
enddo

end subroutine

!----------------------------------------------------------

