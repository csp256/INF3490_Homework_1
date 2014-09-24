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
	

def geneticTSP(names, distances, n, timeAllowed, verbose, varianceLengthDivisor, mutationRateInverse, populationSize):
	t = time.time()
	a = range(0,n)
	bestGenome = a[:]
	bestFitness = routeLength(a, distances)
	parents, parentsFitness, bestGenome, bestFitness = geneticTSP_randomInitialization(distances, n, populationSize, a, bestGenome, bestFitness)
	generations = 0
	dt = 0.0 
	while dt < timeAllowed:
		generations += 1
		children, childrensFitness, bestGenome, bestFitness = geneticTSP_createNewGeneration(parents, parentsFitness, n, bestGenome, bestFitness, distances, mutationRateInverse, varianceLengthDivisor)
		parents, parentsFitness    = geneticTSP_mergeGenerations   (parents, parentsFitness, children, childrensFitness)
		dt = time.time() - t
	bestGenome = standardizePermutation(bestGenome)
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
