import csv
import math
import time
import random
from random import randrange
from generalTSP import *

def geneticTSP_randomInitialization(distances, n, populationSize, a, bestGenome, bestFitness):
	parents = []
	for i in range(0,populationSize):
		parents.append( perm(a, randrange(0, math.factorial(n)-1)) )
	parentsFitness = []
	for i in range(0,populationSize):
		parentsFitness.append(routeLength(parents[i], distances))
		if (parentsFitness[i] < bestFitness):
			bestFitness = parentsFitness[i]
			bestGenome = parents[i][:]
	return parents, parentsFitness, bestGenome, bestFitness
	
def geneticTSP_tournamentSelection(parentsFitness): #draw 2 choose 1
	X = randrange(0, len(parentsFitness))
	Y = randrange(0, len(parentsFitness))
	if (parentsFitness[X] < parentsFitness[Y]):
		return X
	else:
		return Y

def geneticTSP_createChild(parents, parentsFitness, distances, mutationRateInverse, varianceLengthDivisor):
	X = geneticTSP_tournamentSelection(parentsFitness)
	Y = geneticTSP_tournamentSelection(parentsFitness)
	child1 = []
	child2 = []
	for i in range(0, len(parents[0])):
		if (randrange(0,2) == 0):
			child1.append(parents[X][i])
		else:
			child2.append(parents[Y][i])
	child1 += [x for x in parents[Y] if x not in child1]
	child2 += [x for x in parents[X] if x not in child2]
	geneticTSP_mutate(child1, mutationRateInverse, varianceLengthDivisor)
	geneticTSP_mutate(child2, mutationRateInverse, varianceLengthDivisor)
	childFitness1 = routeLength(child1, distances)
	childFitness2 = routeLength(child2, distances)
	if (childFitness1 < childFitness2):
		return child2, childFitness2
	return child1, childFitness1

def geneticTSP_mutate(a, mutationRateInverse, varianceLengthDivisor):
	if (0 == randrange(mutationRateInverse)):
		"""
		r = randrange(len(a))
		a = rotate(a,r)
		perm(a,randrange(len(a)))
		a = rotate(a,-r)
		"""
		createCandidatePermutation(a, varianceLengthDivisor)
		
def geneticTSP_createNewGeneration(parents, parentsFitness, n, bestGenome, bestFitness, distances, mutationRateInverse, varianceLengthDivisor):
	children = []
	childrensFitness = []
	for i in range(0,len(parents)):
		child, childFitness = geneticTSP_createChild(parents, parentsFitness, distances, mutationRateInverse, varianceLengthDivisor)
		if (childFitness < bestFitness):
			bestGenome = child[:]
			bestFitness = childFitness
		children.append(child)
		childrensFitness.append(childFitness)
	return children, childrensFitness, bestGenome, bestFitness

	
def geneticTSP_mergeGenerations(parents,parentsFitness,children,childrensFitness):
	survivors = []
	survivorsFitness = []
	for i in range(0,len(parents)):
		if (childrensFitness[i] < parentsFitness[i]):
			survivors.append(children[i])
			survivorsFitness.append(childrensFitness[i])
		else:
			survivors.append(parents[i])
			survivorsFitness.append(parentsFitness[i])
	return survivors, survivorsFitness
	
def geneticTSP_updateCensus(census, parentsFitness, dt):
	census[0].append(min(parentsFitness))
	census[1].append(max(parentsFitness))
	census[2].append(avg(parentsFitness))
	census[3].append(stdDev(parentsFitness))
	census[4].append(dt)
	return census


def geneticTSP(names, distances, n, timeAllowed, verbose, varianceLengthDivisor, mutationRateInverse, populationSize, figureOffset):
	t = time.time()
	a = range(0,n)
	bestGenome = a[:]
	bestFitness = routeLength(a, distances)
	parents, parentsFitness, bestGenome, bestFitness = geneticTSP_randomInitialization(distances, n, populationSize, a, bestGenome, bestFitness)
	generations = 0
	dt = 0.0 
	bestFitnesses = [math.log(bestFitness)]
	times = [dt]
	census = [[],[],[],[],[]]
	census = geneticTSP_updateCensus(census, parentsFitness, dt)
	while dt < timeAllowed:
		generations += 1
		children, childrensFitness, bestGenome, bestFitness = geneticTSP_createNewGeneration(parents, parentsFitness, n, bestGenome, bestFitness, distances, mutationRateInverse, varianceLengthDivisor)
		parents, parentsFitness = geneticTSP_mergeGenerations(parents, parentsFitness, children, childrensFitness)
		dt = time.time() - t
		census = geneticTSP_updateCensus(census, parentsFitness, dt)
		if math.log(bestFitness) < bestFitnesses[len(bestFitnesses)-1]:
			bestFitnesses.append(math.log(bestFitness))
			times.append(dt)
	bestGenome = standardizePermutation(bestGenome)
	census = geneticTSP_updateCensus(census, parentsFitness, dt)
	bestFitnesses.append(math.log(bestFitness))
	times.append(dt)
	plt.figure(3+figureOffset)
	plt.plot(times,bestFitnesses)
	
	best  = [math.log(x) for x in census[0]]
	worst = [math.log(x) for x in census[1]]
	avg   = [math.log(x) for x in census[2]]
	# dont need log of sigma
	plusOneSigma  = [math.log(census[2][i]+census[3][i]) for i in range(len(census[0]))]
	minusOneSigma = [math.log(census[2][i]-census[3][i]) for i in range(len(census[0]))]
	plt.figure(50+figureOffset)
	p = plt.plot(census[4],avg)
	plt.plot(census[4],plusOneSigma,  '-.', color=p[0].get_color(), alpha=0.6)
	plt.plot(census[4],minusOneSigma, '-.', color=p[0].get_color(), alpha=0.6)
	plt.plot(census[4],best, ':',  color=p[0].get_color(), alpha=0.3)
	plt.plot(census[4],worst,':',  color=p[0].get_color(), alpha=0.3)
	
#	plt.plot(census[4],census[1])
#	plt.plot(census[4],census[0])
	
	if verbose:
		print "By genetic search of",generations,"generations"
		print "of",populationSize,"inidividuals the path"
		print bestGenome
		print "is psuedooptimal. It has a length of"
		print bestFitness
		print "and visits",n,"cities. This took"
		print dt, "seconds."
		print " "
	return bestGenome, bestFitness, generations
