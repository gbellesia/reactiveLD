#!/usr/bin/env python
from multiprocessing import Pool
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import math
import numpy
import random
import matplotlib.pyplot
import scipy.stats

rho = 0.1
total = 2000
X = 100.0
volFracts = numpy.linspace(0.01, 0.3, 25)

repeats = 1

def work(args):
    T, inData = args

    return output

pool = Pool(4)

out = []

dt = 0.001
# For all the possible volume fractions, compute an effective diffusion coefficient
# Volume fraction = area of particle * number of particles / volume_of_surface
for v in volFracts:
    total = int(X * X * v / numpy.pi)
    #print total
    #X = math.sqrt(numpy.pi * total / v)
    #distsSq = []

    printSteps = 10000

    # Run a simulation on an X*X surface
    # Timestep = 0.01
    # Diffusion constant = 0.5
    # Radius = 1
    # 200 timesteps, so run to t = 2.0s
    config = { "X" : X,
               "Y" : X,
               #"Z" : X,
	       "dt" : dt,
               "steps" : 400000,
               "printSteps" : printSteps,
               "monatomic" : [],
               "binary" : [],
               "atomsDist" : [{"type" : 1, "number" : total}],
               "types" : [{ "type" : 1, "radius" : 1.0, "D" : 0.5 }]}

    [infdid, inData] = tempfile.mkstemp(dir = ".")
    
    h = os.fdopen(infdid, "w")
    h.write(json.dumps(config))
    h.close()

    [outfdid, outData] = tempfile.mkstemp(dir = ".")
    os.close(outfdid)
    #print outData
    
    h = subprocess.Popen(shlex.split("./bd_run {0} {1} {2} pizza".format(inData, outData, int(random.random() * 1000000000))), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    h.wait()

    print 'hi'
    
    #print stdout
        #print stdout, stderr
    
    #map of atom id -> trajectories
    pos = {}
    
    h = open(outData, "r")

    # Read in all the output data here. It's a fun format
    # The output of this is 'pos', which is a map of atomid -> trajectories (like stated above)
    # So that's like:
    #
    # { atom1 : [(x0, y0, z0), (x1, y1, z1), ... (xN, yN, zN)],
    #   atom2 : [(x0, y0, z0), (x1, y1, z1), ... (xN, yN, zN)],
    #         .
    #         .
    #         .
    #   atomN : [(x0, y0, z0), (x1, y1, z1), ... (xN, yN, zN)] }
    #
    readTime = False
    readNext = False
    lines = list(h.readlines())
    for line in lines:
        if readTime:
            time = int(line)

            readTime = False
        if readNext:
            ind = line.split()

            if re.match("[0-9]+", ind[0]):
                if time == 0:
                    pos[int(ind[0])] = []

                pos[int(ind[0])].append([float(ind[2]), float(ind[3]), float(ind[4])])
            else:
                readNext = False
        if re.match("ITEM: ATOMS", line):
            readNext = True
        if re.match("ITEM: TIMESTEP", line):
            readTime = True

    h.close()

    dists = {}

    os.remove(outData)

    # Compute a graph of time separation vs. MSD
    #
    #  So in a time series, x0, x1, .... xN
    #
    #  There are N choose 2 unique pairs of positions that can be chosen. All of these are separated by their own (possibly not unique) time separations.
    #
    #  So there are N - 1 pairs of distances separated by a single dt, N - 2 pairs separated by 2 * dt, etc.
    #
    #  ----
    #
    #  In each simulation there are also a bunch of atoms. So if there are M atoms, then we have (N - 1) * M samples of how far an atom diffuses in one timestep. (N - 2) * M samples of how far an atom diffuses in two timesteps, etc.
    #
    #  ----
    #
    #  So for each time step separation (in the range 1 * dt to 0.25 * steps * dt), we compute the squared displacement of *all* atoms and lump this in one distribution, and take the mean (so it's the MSD of every molecular movement with the given time separation).
    #
    #  From the plot of MSD vs. time step separation we interpolate MSD = 4 * D * t and get values for D

    for atomid in pos:
        steps = len(pos[atomid])
        for s in range(1, steps / 4):
            traj = pos[atomid]

            distsSq = []
            for i in range(0, len(traj) - s):
                dx = min(abs(traj[i + s][0] - traj[i][0]), X - abs((traj[i + s][0] - traj[i][0])))
                dy = min(abs(traj[i + s][1] - traj[i][1]), X - abs((traj[i + s][1] - traj[i][1])))
                dz = min(abs(traj[i + s][2] - traj[i][2]), X - abs((traj[i + s][2] - traj[i][2])))

                distSquared = dx * dx + dy * dy + dz * dz
                distsSq.append(distSquared)

            if s not in dists:
                dists[s] = []

            dists[s].extend(distsSq)

    output = {}
    for d in dists:
        output[d] = [numpy.mean(dists[d]), scipy.stats.sem(dists[d])]

    distsSq = output#work((0, inData))

    times, valuesStds = zip(*sorted(distsSq.items(), key = lambda x : x[0]))
    times = numpy.array(times) * printSteps * dt;

    values, stds = zip(*valuesStds)

    values = numpy.array(values)
    stds = numpy.array(stds)

    length = len(times)

    # Compute the slope of MSD = 4 * D * t and divide by 4!

    a, b = numpy.polyfit(times[length / 2 :], values[length / 2 :], 1) / 4
    a2, b = numpy.polyfit(times[length / 2 :], values[length / 2 :] - stds[length / 2 :], 1) / 4
    a3, b = numpy.polyfit(times[length / 2 :], values[length / 2 :] + stds[length / 2 :], 1) / 4

    # Print the results! a is the diffuion constant!

    print "total particles {0}, diffusion constant {1}, volume fraction {2}, side length {3}, diff min {4}, diff max {5}".format(total, a, v, X, a2, a3)


    out.append((total, (values, stds)))
        
    os.remove(inData)
    ciccio = numpy.vstack((values,stds)).T
    numpy.savetxt('diffusion0009.out', ciccio)
for total, dataSet in out:
    matplotlib.pyplot.plot(dataSet[0])
    matplotlib.pyplot.plot(dataSet[0] - dataSet[1])
    matplotlib.pyplot.plot(dataSet[0] + dataSet[1])

matplotlib.pyplot.show()
matplotlib.pyplot.legend(totals)
