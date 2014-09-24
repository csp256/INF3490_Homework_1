import csv
import math
import time
import random
import matplotlib.pyplot as plt
from random import randrange
from generalTSP import *
from geneticTSP import *
 
def main():
	names, distances = readCSV('european_cities.csv')
	verbose = True
	varianceLengthDivisor = 5.0
	mutationRateInverse = 10.0
	
	numberOfRuns = 3
	numberOfCities = [10]#, 24]
	timeLimits = [0.1, 0.3]#, 1.0]#, 3.0, 10.0]
	# This commented decleration of populationSizes serves to remind the
	#   author that the parameter is actually decleared inside the loops below.
	# populationSizes = [math.trunc(n*1.0), math.trunc(n*3.0), math.trunc(n*10.0)]
	

	# heat the processor up a bit first
	# this actually helps achieve more consistent results!
	exhaustiveTSP (names, distances, len(names), 2.0, False)
	
	for timeLimit in timeLimits:
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		for n in numberOfCities:
			populationSizes = [math.trunc(n*1.0), math.trunc(n*3.0), math.trunc(n*10.0)]
			print "***********************************************************************"
			for populationSize in populationSizes:
				print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
				for i in range(numberOfRuns):
					print "-----------------------------------------------------------------------"
					exhaustiveTSP (names, distances, n, timeLimit, verbose) 
					hillClimbTSP  (names, distances, n, timeLimit, verbose, varianceLengthDivisor)
					geneticTSP    (names, distances, n, timeLimit, verbose, varianceLengthDivisor, mutationRateInverse, populationSize)
				#do analysis here
		
	return 0

if __name__ == '__main__':
	main()
