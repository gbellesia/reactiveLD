#!/usr/bin/env python
from multiprocessing import Pool
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import numpy
import random
import matplotlib.pyplot

X = 20

totals = [10, 50, 200, 300, 400, 500]
repeats = 1

def work(args):
    T, inData = args

    return output

pool = Pool(4)

out = []

dt = 0.1
for total in totals:
    #distsSq = []
    config = { "X" : X,
               "Y" : X,
               "Z" : X,
               "dt" : dt,
               "steps" : 201,
               "printSteps" : 1,
               "monatomic" : [],
               "binary" : [],
               "atomsDist" : [{"type" : 1, "number" : total},
                              {"type" : 2, "number" : 0}],
               "types" : [{ "type" : 1, "radius" : 1.0, "D" : 0.5 },
                          { "type" : 2, "radius" : 1.0, "D" : 0.5 }]}

    [infdid, inData] = tempfile.mkstemp(dir = ".")
    
    h = os.fdopen(infdid, "w")
    h.write(json.dumps(config))
    h.close()

    [outfdid, outData] = tempfile.mkstemp(dir = ".")
    os.close(outfdid)
    #print outData
    
    h = subprocess.Popen(shlex.split("./bd_run {0} {1} {2}".format(inData, outData, int(random.random() * 1000000000))), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = h.communicate()
    
        #print stdout, stderr
    
    #map of atom id -> trajectories
    pos = {}
    
    h = open(outData, "r")
    
    readTime = False
    readNext = False
    for line in h.readlines():
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

    for atomid in pos:
        steps = len(pos[atomid])
        for s in range(1, steps / 4):
            traj = pos[atomid]

            distsSq = []
            for i in range(0, steps / 4):
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
        output[d] = numpy.mean(dists[d])

    distsSq = work((0, inData))#map(work, zip(range(0, repeats), repeats * [inData]))

    #for key, val in sorted(distsSq.items(), key = lambda x : x[0]):
    #    print key, val

    times, values = zip(*sorted(distsSq.items(), key = lambda x : x[0]))
    times = numpy.array(times) * dt;

    length = len(times)

    a, b = numpy.polyfit(times[length / 2 :], values[length / 2 :], 1) / 6

    print "total particles {0}, diffusion constant {1}".format(total, a)

    out.append((total, values))

    #meansq = numpy.mean(distsSq)

    #print numpy.sqrt(meansq), numpy.mean(distsSq), numpy.sqrt(numpy.var(distsSq))#, distsSq, 6 * 1e-4 * 1e8 * 0.05e-9 * 100
    #matplotlib.pyplot.hist(distsSq)
    #matplotlib.pyplot.show()
        
    os.remove(inData)

for total, dataSet in out:
    matplotlib.pyplot.plot(dataSet)

matplotlib.pyplot.show()
matplotlib.pyplot.legend(totals)
