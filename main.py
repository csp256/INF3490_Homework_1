import csv
import math
import time
import random
from random import randrange

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
	
def createCandidatePermutation(a,l):
	start = randrange(l)
	end = start + 2 + myNormal(l/5)
	if l <= end:
		end %= l
	if end < start:
		temp = start
		start = end
		end = temp
	subsequenceReversal(a,start,end)

def exhaustiveTSP(names, distances, i): 
	t = time.time()
	a = range(0,i)
	bestPermutation = a[:]
	bestDistance = routeLength(a,distances)
	for x in range(1, math.factorial(i)):
		b = perm(a, x)
		currDistance = routeLength(b,distances)
		if (currDistance < bestDistance):
			bestDistance = currDistance
			bestPermutation = b[:]
	print "By exhaustive search the permutation "
	print bestPermutation
	print "is optimal. It has a path length of "
	print bestDistance
	print "and visits ", i, " cities. This took "
	print (time.time() - t), "seconds."
	print " "

def hillClimb(names, distances, i, droughtLimit):
	t = time.time()
	a = range(0,i)
	bestPermutation = a[:]
	bestDistance = routeLength(a,distances)
	drought = 0 # define the "drought" as the number of consecutive failed attempts to generate a better solution
	while drought < droughtLimit:
		createCandidatePermutation(a, len(a))
		currDistance = routeLength(a, distances)
		if currDistance < bestDistance:
			bestPermutation = a[:]
			bestDistance = currDistance
			drought = 0
		else:
			drought += 1
	print "By hillclimbing the permutation "
	print bestPermutation
	print "is psuedooptimal. It has a path length of "
	print bestDistance
	print "and visits ", i, " cities. This took "
	print (time.time() - t), "seconds."
	print " "
 
def main():
	names, distances = readCSV('european_cities.csv')
	n = 8
	exhaustiveTSP(names, distances, n)
	hillClimb(names,distances, n, 1000)
	
	return 0

if __name__ == '__main__':
	main()

