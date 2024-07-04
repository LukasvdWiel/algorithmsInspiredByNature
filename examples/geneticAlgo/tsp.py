import numpy
import random
import math
import sys

global nCities
global populationSize
global nDeath
global nGenerations
global mutationChance
global nSurvivors

############################################################
# Set these parameters to tune the model
nCities = 100
populationSize = 100
nDeath = 50
nGenerations = 1000
mutationChance = 0.02
############################################################

nSurvivors = populationSize - nDeath

#--------------------------------------------------------

def length(x, y, path):

	length = 0.0

#   print ("lx: ", x)
#   print ("ly: ", y)
#   print ("lp: ", path)



	for iCity in range (0, nCities-1):

		distance = math.sqrt(pow(x[path[iCity]] - x[path[iCity+1]], 2) + \
							 pow(y[path[iCity]] - y[path[iCity+1]], 2))

#	   print("distance", distance)
		length = length + distance
#   endfor

	distance = math.sqrt(pow(x[path[nCities-1]] - x[path[0]], 2) + \
						 pow(y[path[nCities-1]] - y[path[0]], 2))

	length = length + distance

#   print("length of route", length)

	return length

#--------------------------------------------------------

def routeHasCity(route, city):



	for iCity in range(0, nCities):
		if (route[iCity] == city):
			return iCity
#	  endif
#   endfor
	return -999

#--------------------------------------------------------

def initialisePath():

	path = numpy.zeros(nCities)

	for iCity in range(0, nCities):
		path[iCity] = iCity
#	endfor

#	randomize the path by switching a lot of numbers
	for iCity in range(0, nCities**2):
		r = random.uniform(0.0, 1.0)
		n = int(math.ceil(r * nCities-1))
		r = random.uniform(0.0, 1.0)
		m = int(math.ceil(r * nCities-1))
		m = (m+n)%nCities

		temp = path[n]
		path[n] = path[m]
		path[m] = temp
#	endfor

	return path

#--------------------------------------------------------

def reverseRandomSection(path):


#	print("reversing a section in", path)

	nMin = 2
	nMax = int(math.floor(nCities * 0.5))

	r = random.uniform(0.0, 1.0)
	reverseSize = nMin + int(math.ceil(r * (nMax - nMin)))

#	print("reverseSize: ", reverseSize) 

	r			 = random.uniform(0.0, 1.0)
	sectionStart = int(math.ceil(r * nCities))
#	print("start idx:", sectionStart)

	sectionEnd   = sectionStart + reverseSize - 1
#	print("a reverse from", sectionStart, "to", sectionEnd)

	sectionEnd   = sectionEnd%nCities # possibly wrap

#	print("b reverse from", sectionStart, "to", sectionEnd)


# put the to be reversed section in a buffer
	section = numpy.zeros(reverseSize)

	for i in range(0,reverseSize):
		pos = (sectionStart + i)%nCities
#		print("grab pos", pos)
		section[i] = path[pos]
#	endfor

#	print("reverse buffer", section)

# put the buffer back in the section in reverse order
	for i in range(0,reverseSize):
		pos = (sectionStart + i)%nCities
#		print("replace pos", pos, "by buffer entry", reverseSize, i, reverseSize - i-1)


		path[pos] = section[reverseSize - i-1]
#	endfor


#	print("new path with reverted section in it", path)


	return path


#--------------------------------------------------------

def sortLengths(values):
# use a shameless bubble sort

	sortedIndex = numpy.zeros(populationSize, dtype=numpy.int8)
	for iPop in range(0, populationSize,):
		sortedIndex[iPop] = iPop
#	endfor

	for i in range(0, populationSize-1):
		for j in range(0, populationSize-1):
			if (values[j] > values[j+1]):

				# use temp for floating points
				temp = values[j]
				values[j] = values[j+1]
				values[j+1] = temp

				# redefine temp to use for integers
				temp = sortedIndex[j]
				sortedIndex[j] = sortedIndex[j+1]
				sortedIndex[j+1] = temp


#			endif
#		endfor
#	endfor



	return [sortedIndex, values]

#--------------------------------------------------------

def createOffspring(mom, dad):

	# we cannot randomly split and glue together bits of genes,
	# because we must ensure that every city occurs exactly and only once.

	# to do this, we take a random bit from mom, and fill the
	# missing bits up following the sequence from dad.

	# keep in mind, though, the the journey is circular,
	# so the bit inherited from mom can wrap around the integer edge

	child = numpy.zeros(nCities, dtype=numpy.int8)
	child = child - 1

	r = random.uniform(0.0, 1.0)
	nFromMom = int(math.floor(nCities * 0.5))

#	print("# reproduction from mom: ", mom)
#	print("# reproduction from dad: ", dad)



#	print("# take", nFromMom, "numbers from mom")

	r = random.uniform(0.0, 1.0)
	momStart = int(math.ceil(r * nCities))


	momEnd = momStart + nFromMom - 1
	momEnd = momEnd%nCities # possibly wrap


#	print("# taking from mom:", momStart, momEnd)

	if (momEnd > momStart):
		# there was no wrapping
		child[momStart:momEnd+1] = mom[momStart:momEnd+1]
#		print("a copied bit", mom[momStart:momEnd+1])
	else:
		# there is wrapping, copy mom's bit in two takes
		child[momStart:nCities] = mom[momStart:nCities]
		child[0:momEnd+1]		= mom[0:momEnd+1]
#		print("b copied bit", mom[momStart:nCities])
#		print("c copied bit", mom[0:momEnd+1])

	# endif

#	print("# after copying mom bit: ", child)


	# it is possible that momEnd is the end of the array.
	# prevent the new position from running out:
	if (momEnd == nCities):
		fillChildPos = 1
	else:
		fillChildPos = momEnd+1
	# endif

	fillChildPos = fillChildPos%nCities

	# Step through the father and copy unknown entries.
	for i in range(0,nCities):
		checkMe = dad[i]
#		print("check dad pos", i, "which is", dad[i])
		if (-999 == routeHasCity(child, checkMe)):

#			print("city", dad[i], "not yet in solution")

			# yay, we have a city from dad the mother did not give.
			# Add it!
			child[fillChildPos] = checkMe
#			print("added ", checkMe, "to pos", fillChildPos, "so:", child)

			if (fillChildPos == nCities-1):
#				print("set fill pos to 0")
				fillChildPos = 0
			else:
				fillChildPos = fillChildPos+1
#				print("set fill pos to ", fillChildPos)

#			endif
#		else:
#			print("city", dad[i], "already in solution")

#		endif
#	endfor

#	print("# after copying dad bit: ", child)


	return child



#--------------------------------------------------------


shortestLength = 99999.0

# initialise array variables
paths = numpy.zeros([populationSize, nCities], dtype=numpy.int8)
x = numpy.zeros(nCities)
y = numpy.zeros(nCities)
lengths = numpy.zeros(populationSize)
sortedIndex = numpy.zeros(populationSize, dtype=numpy.int8)
shortestPath = numpy.zeros(nCities, dtype=numpy.int8)



# initialize random cities
for iCity in range(0, nCities):
	x[iCity] = random.uniform(0.0, 1.0)
	y[iCity] = random.uniform(0.0, 1.0)
#endfor

# initialize random routes
for iPop in range(0, populationSize):
	paths[iPop,:] = initialisePath()
#endfor

shortestPath = paths[0,:]


# print the data of an initial distribution for plotting (path 1)
fi = open("initialPath.dat", "w")
for jCity in range(0, nCities):
	fi.write(str(x[paths[1,jCity]]) + " " +  str(y[paths[1,jCity]]) + "\n")
# endfor
fi.write(str(x[paths[0]]) + " " + str(y[paths[0]]))
fi.close()

fc = open("convergence.dat", "w")

##########################################################################
# run the actual evolution
##########################################################################

for iGeneration in range(1, nGenerations):

#	print("----------------------------------------------------------------")
#	print("----  running generation ", iGeneration, "----------------------------")
#	print("----------------------------------------------------------------")


#	print("paths in this generation:")
#	for iPop in range(0, populationSize):
#		print("route ", iPop, "is:", paths[iPop,:])



	for iPop in range(0, populationSize):
		lengths[iPop] = length(x, y, paths[iPop,:])
#   endfor

#	for iPop in range(0, populationSize):
#		print("length: ", iPop, "is", lengths[iPop])



# Determine the indices of the specimens that are dying,
# because they have the lowest fitness, and the indices
# of the survivors
	sortedIndex, lengths = sortLengths(lengths)

# if we have a new record breaking route, save it
	if (lengths[0] < shortestLength):
		shortestLength = lengths[0]
		shortestPath = paths[sortedIndex[0],:]

#		print("shortest", shortestPath)

#   endif
	fc.write(str(math.log10(float(iGeneration)))  + " " + str(shortestLength) + "\n")


# fill the vacated spots with offspring from the survivors
	for i in range (0, nDeath):
		# determine position of new child
		childID = sortedIndex[nSurvivors + i]

		# determine positions of parents
		r = random.uniform(0.0, 1.0)
		momID = int(math.ceil(r*nSurvivors-1))
		r = random.uniform(0.0, 1.0)
		dadID = int(math.ceil(r*nSurvivors-1))
		dadID = (dadID + momID)%nSurvivors

		momID = sortedIndex[momID]
		dadID = sortedIndex[dadID]

#		print("make a child ",childID," with mom", momID, "is", paths[momID,:])
#		print("make a child ",childID," with dad", dadID, "is", paths[dadID,:])

		paths[childID,:] = createOffspring(paths[momID,:], paths[dadID,:])

#		print("make child", childID, "which is: ", paths[childID,:])

#	endfor

# apply mutations to all individuals
	for iPop in range(0, populationSize):
		r = random.uniform(0.0, 1.0)
#		print("determine mutation of individual", iPop, r)
		if (r < mutationChance):

			paths[iPop,:] = reverseRandomSection(paths[iPop,:])
#		endif
#	endfor

# endfor
##########################################################################


fc.close()


# and print the data of the best path to file
ff = open("finalSolution.dat", "w")
for iCity in range (0, nCities):
	ff.write(str(x[shortestPath[iCity]])  + " " + str(y[shortestPath[iCity]]) + "\n")
# endfor
# wrap back to city 1
ff.write(str(x[shortestPath[0]])  + " " + str(y[shortestPath[0]]) + "\n")
ff.close





