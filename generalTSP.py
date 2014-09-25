import csv
import math
import time
import random
import matplotlib.pyplot as plt
from random import randrange

def avg(l):
	return sum(l)/len(l)

def stdDev(l):
	mean = avg(l)
	s = 0
	for i in l: # calculates standard deviation of path length
		s += (i - mean)**2
	s /= 1.0*len(l)
	return math.sqrt(s)

def rotate(l,n):
    return l[n:] + l[:n]
    
def standardizePermutation(a):
	a = rotate(a,a.index(0))
	if a[len(a)-1] < a[1]:
		a.reverse()
		a = rotate(a, -1)
	return a

def routeLength(a, distances):
	l = len(a)
	dist = distances[a[l-1]][a[0]]
	for i in range(1, l):
		dist += distances[a[i-1]][a[i]]
	return dist

def perm(sequence, index):
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

def subsequenceReversal(a, start, end):
	a[start:end] = a[start:end][::-1]

def readCSV(filename):
	distance_table = []
	with open(filename) as file:
		reader = csv.reader(file,delimiter=';')
		city_names = reader.next()
		for row in reader:
			distance_table.append([float(cell) for cell in row])
	return (city_names, distance_table)

def myNormal(variance):
	return math.trunc(abs(random.gauss(0,variance)))
	
def createCandidatePermutation(a, varianceLengthDivisor):
	l = len(a)
	start = randrange(l)
	end = start + 2 + myNormal(l / varianceLengthDivisor)
	if l <= end:
		end %= l
	if end < start:
		temp = start
		start = end
		end = temp
	subsequenceReversal(a,start,end)

def exhaustiveTSP(names, distances, i, timeLimit, verbose): 
	t = time.time()
	a = range(0,i)
	bestPermutation = a[:]
	bestDistance = routeLength(a,distances)
	dt = 0.0
	factorial = math.factorial(i)
	x = 0
	bestDistances = [math.log(bestDistance)]
	times = [dt]
	while x < factorial:
		x += 1
		b = perm(a, x)
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
	plt.figure(1)
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
	plt.figure(2)
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
