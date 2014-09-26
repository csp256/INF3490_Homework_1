import csv
import math
import time
import random
from random import randrange
from generalTSP import *

"""
	This file holds the geneticSearch algorithm. It also uses features found in generalTSP.py
	
	We call the routeLength "fitness" here. However, because less routeLength is better, 
	this genetic search algorithm actually selects routes with LOWER "fitness".
"""

def geneticTSP_randomInitialization(distances, n, populationSize, a, bestGenome, bestFitness):
	parents = []
	for i in range(0,populationSize):
		parents.append( perm(a, randrange(0, math.factorial(n)-1)) ) # selects a random permutation with uniform probability
	parentsFitness = []
	for i in range(0,populationSize):
		parentsFitness.append(routeLength(parents[i], distances)) 
		if (parentsFitness[i] < bestFitness):
			bestFitness = parentsFitness[i]
			bestGenome = parents[i][:] # slice notation used to perform deep copy
	return parents, parentsFitness, bestGenome, bestFitness

# because my implementation/analysis focus seems to favor small sizes, maybe we should	
def geneticTSP_tournamentSelection(parentsFitness): # draw 2 choose 1. simplest form of tournament selection.
	X = randrange(0, len(parentsFitness))
	Y = randrange(0, len(parentsFitness))
	if (parentsFitness[X] < parentsFitness[Y]):
		return X
	else:
		return Y

def geneticTSP_createChild(parents, parentsFitness, distances, mutationRateInverse, varianceLengthDivisor):
	X = geneticTSP_tournamentSelection(parentsFitness)
	Y = geneticTSP_tournamentSelection(parentsFitness)
	child1 = [] # holds the child's tour aka genome
	child2 = []
	# iterate through parent 1's tour, one city at a time. 50% of the time, add the city to the
	for i in range(0, len(parents[0])):
		if (randrange(0,2) == 0): # tweaking this probability (and this whole function) is probably the most interesting thing to try next!
			child1.append(parents[X][i])
		else: # i think it would be better to pull whole subsequences from the parents instead of having stochastic 50% chance... maybe a markov chain for if you will pull that subsequence? maybe some greedy reordering of selected subsequences, instead of just jamming them in whatever order...
			child2.append(parents[Y][i])
	# at this point, child1 & child2 have ~50% of their tour defined.
	child1 += [x for x in parents[Y] if x not in child1] # add all remaining cities in the order they appear in the OTHER parent
	child2 += [x for x in parents[X] if x not in child2]
	geneticTSP_mutate(child1, mutationRateInverse, varianceLengthDivisor) # custom form of 2-opt
	geneticTSP_mutate(child2, mutationRateInverse, varianceLengthDivisor)
	childFitness1 = routeLength(child1, distances)
	childFitness2 = routeLength(child2, distances)
	if (childFitness1 < childFitness2): 
		return child2, childFitness2 # this part of the function is really depressing if you think about it too much.
	return child1, childFitness1 # what effect would rotating / mirroring parents (or childrens) tour, perhaps only rarely, have?
	

def geneticTSP_mutate(a, mutationRateInverse, varianceLengthDivisor):
	if (0 == randrange(mutationRateInverse)):
		# i was going to inject some entropy at a relatively small (random) subsection of the tour
		# but it was slower than my custom 2-opt implementation
		# and gave worse slightly worse results for all test cases 
		"""
		r = randrange(len(a))
		a = rotate(a,r)
		perm(a,randrange(2*len(a)))
		a = rotate(a,-r)
		"""
		createCandidatePermutation(a, varianceLengthDivisor)
		
def geneticTSP_createNewGeneration(parents, parentsFitness, n, bestGenome, bestFitness, distances, mutationRateInverse, varianceLengthDivisor):
	children = []
	childrensFitness = []
	# create as many children as there are parents
	for i in range(0,len(parents)):
		child, childFitness = geneticTSP_createChild(parents, parentsFitness, distances, mutationRateInverse, varianceLengthDivisor)
		if (childFitness < bestFitness):
			bestGenome = child[:] # slice notation for deep copy
			bestFitness = childFitness
		children.append(child)
		childrensFitness.append(childFitness)
	return children, childrensFitness, bestGenome, bestFitness

	
def geneticTSP_mergeGenerations(parents,parentsFitness,children,childrensFitness):
	survivors = []
	survivorsFitness = []
	# every parent pairs up with a child.
	# only one of them is selected. 
	# the lists are not sorted because the children were selected, and thus ordered, in unbaised order
	for i in range(0,len(parents)):
		if (childrensFitness[i] < parentsFitness[i]):
			survivors.append(children[i])
			survivorsFitness.append(childrensFitness[i])
		else:
			survivors.append(parents[i])
			survivorsFitness.append(parentsFitness[i])
	return survivors, survivorsFitness
	
def geneticTSP_updateCensus(census, parentsFitness, dt):
	census[0].append(min(parentsFitness)) # ugly flat array styled storing
	census[1].append(max(parentsFitness)) # what which index corresponds to is straight-forward
	census[2].append(avg(parentsFitness))
	census[3].append(stdDev(parentsFitness))
	census[4].append(dt)
	return census


def geneticTSP(names, distances, n, timeAllowed, verbose, varianceLengthDivisor, mutationRateInverse, populationSize, figureOffset):
	t = time.time()
	a = range(0,n) # trivial permutation; that is, the worst starting path
	bestGenome = a[:] # slice notation lets us create a copy
	bestFitness = routeLength(a, distances) # just setting initial states
	parents, parentsFitness, bestGenome, bestFitness = geneticTSP_randomInitialization(distances, n, populationSize, a, bestGenome, bestFitness) # while very verbose, this does exactly what it says it does
	generations = 0
	dt = 0.0 
	bestFitnesses = [math.log(bestFitness)] # natural logarithm. this is used to create the plots of the form "geneticSearch_###_populationMuliplier.png"
	times = [dt] # i guess this gets stored twice... whoops. once here and in census.
	census = [[],[],[],[],[]] # this is used to create plots of the population diversity as a function of time. 
	census = geneticTSP_updateCensus(census, parentsFitness, dt) # see this function for a clear understanding of what is being stored, where
	while dt < timeAllowed: # looping off of generations / iterations / whatever does not allow a fair comparison... so we use a constant running time for each search type
		generations += 1
		children, childrensFitness, bestGenome, bestFitness = geneticTSP_createNewGeneration(parents, parentsFitness, n, bestGenome, bestFitness, distances, mutationRateInverse, varianceLengthDivisor)
		parents, parentsFitness = geneticTSP_mergeGenerations(parents, parentsFitness, children, childrensFitness)
		dt = time.time() - t
		if True: #generations % 180 == 0: # if the simulation runs too long you can use this to tune the resolution of the census data, to prevent terrible things from happening to matplotlib
			census = geneticTSP_updateCensus(census, parentsFitness, dt)
		logBestFitness = math.log(bestFitness)
		if logBestFitness < bestFitnesses[len(bestFitnesses)-1]: # used for plotting best-ever member, as a function of time. also note that we consider lower "fitness" as better (because it is really just path length). yes, this breaks the analogy, but the code is still easy to understand.
			bestFitnesses.append(logBestFitness)
			times.append(dt)
	bestGenome = standardizePermutation(bestGenome) # path is rotated until city 0 is at head of list. if next element is greater than the tail of the list, the list is reverted. this transformation does not change path length/fitness.
	census = geneticTSP_updateCensus(census, parentsFitness, dt) # add our last piece of data.
	bestFitnesses.append(math.log(bestFitness)) # particularly useful if we have not had a new best fitness in a while... elsewise the plotted curve terminates before the end time
	times.append(dt)
	plt.figure(3+figureOffset) # for magic number justification see main.py's plotting notes
	plt.plot(times,bestFitnesses) # remember, this just adds it to the figure, does not actually save the figure to an image file
	
	best  = [math.log(x) for x in census[0]] # natural logarithm
	worst = [math.log(x) for x in census[1]] # i like list comprehensions! i wonder if you can write a C++ macro to do the same thing?
	avg   = [math.log(x) for x in census[2]]
	# dont need log of sigma
	plusOneSigma  = [math.log(census[2][i]+census[3][i]) for i in range(len(census[0]))]
	minusOneSigma = [math.log(census[2][i]-census[3][i]) for i in range(len(census[0]))] 
	plt.figure(50+figureOffset) # we are going to plot the average of all individuals fitness
	p = plt.plot(census[4],avg) # we grab the handle p so we can reuse its color
	plt.plot(census[4],plusOneSigma,  '-.', color=p[0].get_color(), alpha=0.6) # faintly plot +/- 1 standard deviation
	plt.plot(census[4],minusOneSigma, '-.', color=p[0].get_color(), alpha=0.6)
	plt.plot(census[4],best, ':',  color=p[0].get_color(), alpha=0.3) # very faintly plot the best and worst
	plt.plot(census[4],worst,':',  color=p[0].get_color(), alpha=0.3) 
	
	if verbose: # not really useful anymore
		print "By genetic search of",generations,"generations"
		print "of",populationSize,"inidividuals the path"
		print bestGenome
		print "is psuedooptimal. It has a length of"
		print bestFitness
		print "and visits",n,"cities. This took"
		print dt, "seconds."
		print " "
	return bestGenome, bestFitness, generations
