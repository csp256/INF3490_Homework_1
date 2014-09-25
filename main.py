import os
import sys
import csv
import math
import time
import random
import matplotlib.pyplot as plt
from random import randrange
from generalTSP import *
from geneticTSP import *

"""
		ABANDON ALL HOPE YOU WHO ENTER HERE
"""

def createDirectory(n,timeLimit):
	directory = "Output_Images/"+str(n)+"_cities/"
	if timeLimit < 1.0:
		directory += str(int(timeLimit*1000))+"_milliseconds/"
	else:
		directory += str(int(timeLimit))+"_seconds/"
	if not os.path.exists(directory):
		os.makedirs(directory)
	return directory

 
def main():
	names, distances = readCSV('european_cities.csv')
	verbose = False
	varianceLengthDivisor = 5.0
	mutationRateInverse = 50.0
	preHeatCPU = True
	numberOfSearchTypes = 3
	
	numberOfRuns = 20
	numberOfCities = [10, 24]
	timeLimits = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
	populationMultipliers = [0.25, 0.5, 1.0, 2.0]
	# This commented decleration of populationSizes serves to remind the
	#   author that the parameter is actually decleared inside the loops below.
	# populationSizes = [math.trunc(n*1.0), math.trunc(n*3.0), math.trunc(n*10.0)]
	
	print "Run time will be at least",
	prediction = numberOfRuns*len(numberOfCities)*sum(timeLimits)*(2+len(populationMultipliers))
	print prediction,
	print "seconds."
	print " "
	t = time.time()
	
	# heat the processor up a bit first...
	# this can actually help achieve more consistent results!
	if preHeatCPU:
		exhaustiveTSP (names, distances, len(names), 60.0, False)
			
	for timeLimit in timeLimits:
		for n in numberOfCities:
			print "timeLimit =",timeLimit," ::  cities = ",n
			populationSizes = []
			for i in populationMultipliers:
				populationSizes.append(int(n*i))
			data          = [ [[],[],[]],   [[],[],[]] ]
			processedData = [ [],   [] ]
			for i in populationSizes:
				data.append( [[],[],[]] )
				processedData.append( [] )
			# each of the search types returns a third variable, depending upon which search.
			#   exhaustive  -->  percent completed
			#   hillClimb   -->  iterations
			#   genetic     -->  generations
			for runNumber in range(numberOfRuns):
				sys.stdout.write(str(runNumber+1)+" ")
				sys.stdout.flush()
				datum = exhaustiveTSP (names, distances, n, timeLimit, verbose) 
				for argumentNumber in range(len(datum)):
					data[0][argumentNumber].append(datum[argumentNumber])

				datum = hillClimbTSP  (names, distances, n, timeLimit, verbose, varianceLengthDivisor)
				for argumentNumber in range(len(datum)):
					data[1][argumentNumber].append(datum[argumentNumber])

				for populationSize in populationSizes:
					figureOffset = populationSizes.index(populationSize)
					datum = geneticTSP    (names, distances, n, timeLimit, verbose, varianceLengthDivisor, mutationRateInverse, populationSize, figureOffset)
					for argumentNumber in range(len(datum)):
						data[2+figureOffset][argumentNumber].append(datum[argumentNumber])
			print " "
			for searchType in range(2 + len(populationSizes)):
				processedData[searchType].append(min   (data[searchType][1])) # extracts minimum path length
				processedData[searchType].append(max   (data[searchType][1])) # extracts maximum path length
				processedData[searchType].append(avg   (data[searchType][1])) # calculates average path length
				processedData[searchType].append(stdDev(data[searchType][1])) # calculates standard deviation of path length
				processedData[searchType].append(avg   (data[searchType][2])) # average percent complete / iterations / generations
			for f in range(1,3+len(populationSizes)):
				plt.figure(f)
				plt.xlabel("Time Elapsed")
				plt.ylabel("Best known length (logarithmic)")
				plt.xlim(xmax=timeLimit)
				if 1 == f:
					title = "Exhaustive Search: "
				if 2 == f:
					title = "Hill Climbing: "
				if 3 <= f:
					title = "Genetic Search: "
				title += str(n)+" Cities, "+str(timeLimit)+" Seconds"
				if 1 == f:
					title += ", ~"+str(processedData[f-1][4])+"% Complete"
				if 2 == f:
					title += ", ~"+str(int(processedData[f-1][4]))+" Iterations"
				if 3 <= f:
					title += ", "+str(round(populationMultipliers[f-3],2))+" Pop. Mult."
					title += ", ~"+str(int(processedData[f-1][4]))+" Generations"
				title += ".\n"
				title += "Best: "   + str(int(round(processedData[f-1][0]))) + ". "
				title += "Worst: "  + str(int(round(processedData[f-1][1]))) + ". "
				title += "Avg: "    + str(int(round(processedData[f-1][2]))) + ". "
				title += "StdDev: " + str(int(round(processedData[f-1][3]))) + "."
				plt.title(title)
				if 1 == f:
					name = "exhaustiveSearch"
				if 2 == f:
					name = "hillClimbing"
				if 3 <= f:
					name = "geneticSearch"
					name += "_"+str(round(populationMultipliers[f-3],2))+"_populationMuliplier"
				directory = "Output_Images/"+str(n)+"_cities/"
				if timeLimit < 1.0:
					directory += str(int(timeLimit*1000))+"_milliseconds/"
				else:
					directory += str(int(timeLimit))+"_seconds/"
				if not os.path.exists(directory):
					os.makedirs(directory)
				plt.savefig(directory+name+".png")
				plt.clf()
			for f in range(len(populationSizes)):
				plt.figure(f+50)
				plt.xlabel("Time Elapsed")
				plt.ylabel("Fitness")
				plt.xlim(xmax=timeLimit)
				title = "Genetic Search: Population Diversity\n"
				title += str(n)+" Cities, "+str(timeLimit)+" Seconds"
				title += ", "+str(round(populationMultipliers[f],2))+" Population Multiplier"
				title += ".\n"
				plt.title(title)
				name = "populationDiversity"
				name += "_"+str(round(populationMultipliers[f],2))+"_populationMuliplier"
				directory = createDirectory(n, timeLimit)
				plt.savefig(directory+name+".png")
				plt.clf()
	print " "
	print "Run time was",round(time.time()-t,2),"seconds."
	print round(100.*(time.time()-t)/prediction,2),"% of prediction."
	return 0

if __name__ == '__main__':
	main()
