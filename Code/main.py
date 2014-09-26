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
		This program is for the fulfilment of Project 1 in
		INF3490, fall 2014, at the University of Oslo.

		Output is directed to the Output_Images directory.
		This probably won't work on Windows. Written for Linux.
"""

# it is much easier to make meaningful comparisons between similar plots if they are sorted into directories
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
	names, distances = readCSV('european_cities.csv') # 
	verbose = False # outputs data after each run. not needed if graphing is working (which it is)
	varianceLengthDivisor = 5.0 # Determines the distribution of widths chosen for stochastic 2-opt
	mutationRateInverse = 50.0 # Probability of mutation is 1/mutationRateInverse
	preHeatCPU = False # shorter runs happen first. I want the first run to be comparable to the last.
	numberOfSearchTypes = 3 # the necessity of how this variable is used is due to poor planning on my part
	
	numberOfRuns = 10
	numberOfCities = [17]#10, 24] 
	timeLimits = []#[0.2, 0.4]#, 1.0] # seconds
	populationMultipliers = [0.25, 0.5, 1.0]#, 2.0] # not static pop sizes, but this multiplies the number of cities to yield the pop size
	print " "
	print "Run time will be at least",
	prediction = numberOfRuns*len(numberOfCities)*sum(timeLimits)*(2+len(populationMultipliers))
	print prediction,
	print "seconds."
	print "(",round(prediction/60.,1),"minutes,",round(prediction/3600.,1),"hours)"
	print " "
	t = time.time()
	
	# heat the processor up a bit first...
	# this can actually help achieve more consistent results!
	# especially since shorter runs are done first. 
	if preHeatCPU:
		exhaustiveTSP (names, distances, len(names), 60.0, False) # run for one minute. excessive.

	exhaustiveTSP (names, distances, 10, 9999.9, True)

	# every piece of my code related to graphing is... suboptimal.
	# the actual search algorithms are better.
	for timeLimit in timeLimits:
		for n in numberOfCities:
			print "timeLimit =",timeLimit," ::  cities = ",n
			populationSizes = [] # we construct this explicitly in the loop...
			for i in populationMultipliers: # because it is a function of numberOfCities
				populationSizes.append(int(n*i))
			# the "data structure" i used for storing and processing the returned value of the is the ugliest part.
			# each search returns numberOfRuns-many copies of 3 variables: best permutation, best length, and...
			# each of the search types returns a third variable, depending upon which search.
			#   exhaustive  -->  percent completed
			#   hillClimb   -->  iterations
			#   genetic     -->  generations
			# the length of data is 2+populationSizes, where the 2 is from exhaustive & hillClimbing
			#   this explains the base case of two sets of lists, and the for loop appending to them below
			# indexing into data is thus done as:
			#   data[searchType][returnedValueType][runNumber]
			data          = [ [[],[],[]],   [[],[],[]] ] 
			# processed data collapses the innermost list (of length runNumber) into a flat list of 5 things
			#    0  -->  minimum path length of all runs
			#    1  -->  maximum path length of all runs 
			#    2  -->  avg path length
			#    3  -->  standard deviation of path (not really a useful metric with our sample size and data characteristics... but at least it is easy to compute..right??)
			#    4  -->  average of the third value returned (outlined in comments directly above; a function of search type)
			processedData = [ [],   [] ]
			for i in populationSizes:
				data.append( [[],[],[]] )
				processedData.append( [] )
			for runNumber in range(numberOfRuns):
				sys.stdout.write(str(runNumber+1)+" ") # cant use print because we want it all on one line
				sys.stdout.flush() # and to show up immediately
				datum = exhaustiveTSP (names, distances, n, timeLimit, verbose) # datum is the singular of data. stores the values returned by the search
				for argumentNumber in range(len(datum)): # there is probably a more python-esque way of doing this. send me an email if you know how?
					data[0][argumentNumber].append(datum[argumentNumber])

				datum = hillClimbTSP  (names, distances, n, timeLimit, verbose, varianceLengthDivisor)
				for argumentNumber in range(len(datum)):
					data[1][argumentNumber].append(datum[argumentNumber])

				for populationSize in populationSizes: # genetic search gets run len(populationSize) many times
					figureOffset = populationSizes.index(populationSize)
					datum = geneticTSP    (names, distances, n, timeLimit, verbose, varianceLengthDivisor, mutationRateInverse, populationSize, figureOffset)
					for argumentNumber in range(len(datum)):
						data[2+figureOffset][argumentNumber].append(datum[argumentNumber]) # magic number accounts for exhaustive and hillClimbing being stored in same data structure
			print " "
			for searchType in range(2 + len(populationSizes)):
				processedData[searchType].append(min   (data[searchType][1])) # extracts minimum path length
				processedData[searchType].append(max   (data[searchType][1])) # extracts maximum path length
				processedData[searchType].append(avg   (data[searchType][1])) # calculates average path length
				processedData[searchType].append(stdDev(data[searchType][1])) # calculates standard deviation of path length
				processedData[searchType].append(avg   (data[searchType][2])) # average percent complete / iterations / generations
			# this loop plots everything except genetic search population diversity
			for f in range(1,3+len(populationSizes)): # by convention figure numbers start with 1, instead of 0. the magic number justification is the same as above.
				plt.figure(f) # plotting to 
				plt.xlabel("Time Elapsed")
				plt.ylabel("Best known length (logarithmic)") # natural logarithm
				plt.xlim(xmax=timeLimit) # some data tends to go slightly past this, but it looks much better to have the axis terminate where the plot title says it should
				# create figure-number-dependent header
				# (i could not get plt.text() or similar functions to do anything! this was the best option.)
				if 1 == f:
					title = "Exhaustive Search: "
				if 2 == f:
					title = "Hill Climbing: "
				if 3 <= f:
					title = "Genetic Search: "
				title += str(n)+" Cities, "+str(timeLimit)+" Seconds"
				if 1 == f:
					title += ", ~"+str(100.*processedData[f-1][4])+"% Complete"
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
				# much easier to compare similar images if they are sorted into directories!
				directory = createDirectory(n, timeLimit)
				plt.savefig(directory+name+".png") # PNG seems to work best, although it is lossy, but GIF might have been better if it was supported!
				plt.clf() # clear the figure when we are done
			# population diversity (as a function of time) of each genetic search
			for f in range(len(populationSizes)):
				plt.figure(f+50) # we use a large figure offset spacing, 50, so we can create both figures in tandem, before they are saved to disk
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
