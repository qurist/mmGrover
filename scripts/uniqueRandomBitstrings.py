# coding: utf-8

from random import *
import sys
def randbitstr(n): return [randint(0,1) for b in range(1,n+1)]
def distinct(a, b, n): return bool(sum([a[i]^b[i] for i in range(n)]))
def isnew(a,strlist): 
	if not strlist: return 1
	else: return reduce(lambda x,y: x and y, [distinct(a,i,len(a)) for i in strlist])

n=int(sys.argv[1])
m=int(sys.argv[2])

count=0
strings = []

while count!=m:
    a = randbitstr(n)
    if isnew(a,strings): 
        strings.append(a)
        count+=1

# print output as a Grover instance
print("c Marked items generated uniformly at random.")
print("p %d %d 1"%(n,m))
for i in range(m):
 	s = ''.join(str(j) for j in strings[i])
	print("%s %d"%(s,-1))

