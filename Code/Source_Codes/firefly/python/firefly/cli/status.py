import argparse
import matplotlib.pyplot as plt

from ..pso import PSO



# parse command line arguments
parser = argparse.ArgumentParser(prog="firefly status", description="plot status of optimization")
parser.add_argument("filename")
parser.add_argument("-L", "--logscale", action="store_true")
args = parser.parse_args()



# load data
particles, gbest, xstd = PSO.loadnc(args.filename)
i = [G.iterations for G in gbest]
deltas = [G.delta for G in gbest]
values = [G.value for G in gbest]



# 1. Normalized swarm radius: standard deviation of particle positions
fig = plt.figure()
plt.plot(xstd)
if args.logscale:
    plt.xscale('log')
    plt.yscale('log')
plt.title("Normalized standard deviation of particle positions")



# 2. Relative change of best 
fig = plt.figure()
plt.plot(i, deltas)
if args.logscale:
    plt.xscale('log')
    plt.yscale('log')
plt.title("Relative change of best known position")



# 3. Best known value
fig = plt.figure()
plt.plot(i, values)
if args.logscale:
    plt.xscale('log')
    plt.yscale('log')
plt.title("Best known value")

plt.show()
