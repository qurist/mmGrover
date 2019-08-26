# In this script we will generate a small cyclic grover schedule for a given input size n

from numpy import pi
import sys

n = int(sys.argv[1])

beta   = pi/n
gamma  = pi
T = int(sys.argv[2])
t = int(sys.argv[3])
cyclic = int(sys.argv[4])


print("c Small cyclic QAOA schedule for Grover on %d bits"%n)
print("p %d %d %d"%(T, t, cyclic))
for i in range(T):
    print("%f %f"%(beta, gamma))
