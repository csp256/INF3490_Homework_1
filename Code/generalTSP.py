import csv
import math
import time
import random
import matplotlib.pyplot as plt
from random import randrange

"""
		This file holds general functions, as well as the 
		exhaustiveSearh & hillClimbing functions.
"""

def avg(l): # average
	return sum(l)/len(l)

def stdDev(l): # standard deviation
	mean = avg(l)
	s = 0
	for i in l: # calculates standard deviation of path length
		s += (i - mean)**2
	s /= 1.0*len(l)
	return math.sqrt(s)

def rotate(l,n): # cyclic list rotation
    return l[n:] + l[:n]
    
def standardizePermutation(a): # path is rotated until city 0 is at head of list. if next element is greater than the tail of the list, the list is reverted. this transformation does not change path length/fitness.
	a = rotate(a,a.index(0))
	if a[len(a)-1] < a[1]:
		a.reverse()
		a = rotate(a, -1)
	return a

def routeLength(a, distances):
	l = len(a)
	dist = distances[a[l-1]][a[0]] # let our base case be our edge case: returning to initial city
	for i in range(1, l):
		dist += distances[a[i-1]][a[i]]
	return dist

# this function was found on StackOverflow
def perm(sequence, index): # for a sequence (list) of length N, and integer 0<= index < factorial(N) this returns a unique permutation
    sequence = list(sequence)
    result = []
    for x in xrange(len(sequence)):
        idx = index % len(sequence)
        index /= len(sequence)
        result.append( sequence[idx] )
        # constant time non-order preserving removal
        sequence[idx] = sequence[-1]
        del sequence[-1]
    return result

# reverses a subsection of the list. does not include end.
def subsequenceReversal(a, start, end):
	a[start:end] = a[start:end][::-1]

# provided code wrapped in a function
def readCSV(filename):
	distance_table = []
	with open(filename) as file:
		reader = csv.reader(file,delimiter=';')
		city_names = reader.next()
		for row in reader:
			distance_table.append([float(cell) for cell in row])
	return (city_names, distance_table)

# integer truncation (round towards 0) of (absolute value of (normal distribution, 0 mean, with given variance.) )
# most commonly returned value is 0 for all values of variance
# this tends towards a fairly greedy search
def myNormal(variance):
	return math.trunc(abs(random.gauss(0,variance)))
	
# reverse a subsection of the tour. the start of the subsection is uniformly chosen.
# the length of the subsection is stochastic and grows according to the tour length. see myNormal()
# also handles all edge cases
def createCandidatePermutation(a, varianceLengthDivisor):
	l = len(a)
	start = randrange(l)
	end = start + 2 + myNormal(l / varianceLengthDivisor) # magic number insures that no trivial reversals occur
	if l <= end:
		end %= l
	if end < start:
		temp = start
		start = end
		end = temp
	subsequenceReversal(a,start,end)

# memory complexity is O(n), not O(n!) as in the python combinatorics function enumerations
# but is ultimately a bit slower
def exhaustiveTSP(names, distances, i, timeLimit, verbose): 
	t = time.time()
	a = range(0,i) # trivial (uninformed) tour
	bestPermutation = a[:]
	bestDistance = routeLength(a,distances)
	dt = 0.0
	factorial = math.factorial(i)
	x = 0
	bestDistances = [math.log(bestDistance)]
	times = [dt]
	while x < factorial:
		x += 1
		b = perm(a, x) # next tour in the sequence. see perm()
		currDistance = routeLength(b,distances)
		if (currDistance < bestDistance):
			bestDistance = currDistance
			bestPermutation = b[:]
		dt = (time.time() - t)
		bestDistances.append(math.log(bestDistance))
		times.append(dt)
		if timeLimit <= dt:# and False: # "and False" turns off premature exiting for testing
			bestPermutation = standardizePermutation(bestPermutation)
			plt.figure(1)
			plt.plot(times,bestDistances)
			if verbose:
				print "Exhaustive search exceeded the time limit of",timeLimit,"seconds."
				print "The best permutation (so far) is"
				print bestPermutation
				print "with a path length of"
				print bestDistance
				print "while visiting",i,"cities."
				print x,"paths were considered. (",x/(1.0*math.factorial(i)),"% of total)"
				print " "
			return bestPermutation, bestDistance, x/(1.0*math.factorial(i))
	bestPermutation = standardizePermutation(bestPermutation)
	plt.figure(1) # see main.py discussion of plotting conventions for magic number justification
	plt.plot(times,bestDistances)
	if verbose:
		print "By exhaustive search the permutation "
		print bestPermutation
		print "is optimal. It has a path length of "
		print bestDistance
		print "and visits ", i, " cities. This took "
		print dt, "seconds."
		print "(",dt/timeLimit,"% of timeLimit,",timeLimit,")"
		print " "
	return bestPermutation, bestDistance, 1.0
	

def hillClimbTSP(names, distances, i, timeLimit, verbose, varianceLengthDivisor):
	t = time.time()
	a = range(0,i)
	bestPermutation = a[:]
	bestDistance = routeLength(a,distances)
	dt = 0.0
	iterations = 0
	bestDistances = [math.log(bestDistance)]
	times = [dt]
	while dt < timeLimit:
		iterations += 1
		createCandidatePermutation(a, varianceLengthDivisor)
		currDistance = routeLength(a, distances)
		dt = (time.time() - t)
		if currDistance < bestDistance:
			bestPermutation = a[:]
			bestDistance = currDistance
			bestDistances.append(math.log(bestDistance))
			times.append(dt)
	bestPermutation = standardizePermutation(bestPermutation)
	bestDistances.append(math.log(bestDistance))
	times.append(dt)
	plt.figure(2) # see main.py discussion of plotting conventions for magic number justification
	plt.plot(times, bestDistances)
	if verbose:
		print "By hillclimbing the permutation "
		print bestPermutation
		print "is psuedooptimal. It has a path length of "
		print bestDistance
		print "and visits ", i, " cities. This took "
		print dt, "seconds."
		print " "
	return bestPermutation, bestDistance, iterations
