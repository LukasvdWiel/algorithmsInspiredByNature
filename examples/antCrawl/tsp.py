import numpy
import random
import math
import sys

global nCities
global nAnts

############################################################
# Set these parameters to tune the model
nCities = 100
nAnts = 10000
pheromoneEvaporation = 0.999
############################################################

def length(x, y, path):

	length = 0.0

	for iCity in range (0, nCities-1):

		distance = math.sqrt(pow(x[path[iCity]] - x[path[iCity+1]], 2) + \
							 pow(y[path[iCity]] - y[path[iCity+1]], 2))

		length = length + distance
#	endfor

	distance = math.sqrt(pow(x[path[nCities-1]] - x[path[0]], 2) + \
						 pow(y[path[nCities-1]] - y[path[0]], 2))

	length = length + distance

	return length


def routeHasCity(route, city):

	for iCity in range(0, nCities):
		if (route[iCity] == city):
			return iCity
#	   endif
#   endfor
	return 0



# declare variables

pheromones = numpy.zeros([nCities, nCities])
distance = 0.0
distances = numpy.zeros([nCities, nCities])
x = numpy.zeros(nCities)
y = numpy.zeros(nCities)
howNice = numpy.zeros(nCities)
progress = numpy.zeros(nAnts)

iCity = 0
jCity = 0
iAnt = 0

totalLength = 0.0
addValue = 0.0
max = 0.0
r = 0.0
selection = 0.0

togoEntry = 0
weAreHere = 0
destination = 0
ffrom = 0
to = 0

shortestLength = 99999.0

# initialize random cities
for iCity in range(0, nCities):
	x[iCity] = random.uniform(0.0, 1.0)
	y[iCity] = random.uniform(0.0, 1.0)
# endfor

for iCity in range(0, nCities):
	for jCity in range(0, nCities):
		distance = math.sqrt(pow(x[iCity] - x[jCity],2) + \
							 pow(y[iCity] - y[jCity],2))
		distances[iCity, jCity] = distance
		distances[jCity, iCity] = distance
#	endfor
# endfor

# set base pheronomes
for iCity in range(0, nCities):
	for jCity in range(0, nCities):
		pheromones[iCity, jCity] = 1.0
#	endfor
# endfor

# let the ants walk
for iAnt in range (0, nAnts):
	path = [0 for n in range(nCities)]
	# not needed to define the start city, as it is 0 here, in stead of 1

	for iCity in range(1, nCities):

		weAreHere = path[iCity-1]
		togo = [0 for n in range(nCities)]
		togoEntry = 0

		# determine the next city, by
		# 1) make a list of all the still unvisited cities
		for jCity in range(1, nCities):
			if (0 == routeHasCity(path, jCity)):
				# yay, we found a city we have not visited yet

				togoEntry = togoEntry + 1
				togo[togoEntry-1] = jCity
#			endif
#		endfor

		# 2) check how nice every city is, based on distance
		#	and deposited pheromones
		howNice = numpy.zeros(nCities)
		for jCity in range(0, togoEntry):
			destination = togo[jCity]

			howNice[jCity+1] = pheromones[weAreHere, destination] / \
						 		distances[weAreHere, destination]
#		endfor

		for jCity in range(1, togoEntry+1):
			howNice[jCity] = howNice[jCity-1] + \
							 howNice[jCity]
# 		endfor

		r = random.uniform(0.0, 1.0)
		selection = r * howNice[togoEntry]

		for jCity in range(0, togoEntry):
			if ((selection > howNice[jCity]) and \
				(selection < howNice[jCity+1])):

				path[iCity] = togo[jCity]

#			endif
#		endfor
#	endfor

# write initial solution to file, for plotting to compare with the final solution
	if (0 == iAnt):
		fi = open("initialPath.dat", "w")
		for jCity in range(0, nCities):
			fi.write(str(x[path[jCity]]) + " " +  str(y[path[jCity]]) + "\n")
#		endfor
		fi.write(str(x[path[0]]) + " " + str(y[path[0]]))
		fi.close()
#	endif

	totalLength = length(x, y, path)
	if (totalLength < shortestLength):
		shortestLength = totalLength
		shortestPath = path
#	endif

# update the pheromones
	progress[iAnt] = shortestLength
	addValue = 1.0 / totalLength

	for iCity in range(2, nCities):

		ffrom = path[iCity-1]
		to   = path[iCity]

		pheromones[ffrom, to] = pheromoneEvaporation * pheromones[ffrom, to] + addValue
		pheromones[to, ffrom] = pheromoneEvaporation * pheromones[to, ffrom] + addValue
#	endfor

# endfor


ff = open("finalSolution.dat", "w")
for iCity in range(0, nCities):
	ff.write(str(x[shortestPath[iCity]]) + " " + str(y[shortestPath[iCity]]) + "\n")
#endfor
ff.close

fc = open("convergence.dat", "w")
for iAnt in range(0, nAnts):
	fc.write(str(iAnt) + " " + str(progress[iAnt]) + "\n")
#endfor
fc.close
