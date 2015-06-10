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
import pickle
import random
import matplotlib.pyplot

totals = [512]

rho = 0.1
X = (8.0 * totals[0] / rho)**(1.0/3.0) 
# from normalized density = (2r)^3 * N / V 
repeats = 1

def work(args):
    T, inData = args

    return output

pool = Pool(4)

out = []

dt = 0.001
realizations = 100
allDistances = []
for r in range(0, realizations):
    distsSq = []
    config = { "X" : 1,
               "Y" : 1,
               "Z" : 1,
               "steps" : 100001,
               "printSteps" : 1000,
               "dt" : 0.0001,
               "simType" : "boxWithWalls",
               "atoms" : [{ "type" : 1, "x" : 0.5, "y" : 0.5, "z" : 0.5 }],
               "atomsDist" : [{ "type" : 2, "number" : 500 },
                              { "type" : 3, "number" : 500 },
                              { "type" : 4, "number" : 1000 }],
               "types" : [{ "type" : 1, "radius" : 0.001, "D" : 1e-3 },
                          { "type" : 2, "radius" : 0.001, "D" : 1e-3 },
                          { "type" : 3, "radius" : 0.0025, "D" : 1e-3 },
                          { "type" : 4, "radius" : 0.0050, "D" : 1e-3 }]}

    [infdid, inData] = tempfile.mkstemp(dir = ".")
    
    print r
    h = os.fdopen(infdid, "w")
    h.write(json.dumps(config))
    h.close()

    [outfdid, outData] = tempfile.mkstemp(dir = ".")
    os.close(outfdid)
    
    h = subprocess.Popen(shlex.split("./bd_run {0} {1} {2} pizza".format(inData, outData, int(random.random() * 1000000000))), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = h.communicate()
    
    #map of atom id -> trajectories
    pos = []
    
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
                if int(ind[1]) == 1:
                    pos.append([float(ind[2]), float(ind[3]), float(ind[4])])
            else:
                readNext = False
        if re.match("ITEM: ATOMS", line):
            readNext = True
        if re.match("ITEM: TIMESTEP", line):
            readTime = True

    h.close()

    distances = []
    
    for f in pos:
        d = numpy.linalg.norm(numpy.array([0.5, 0.5, 0.5]) - numpy.array(f))**2.0
        distances.append(d)

    allDistances.append(distances)

    os.remove(inData)
    os.remove(outData)

print numpy.array(allDistances)
ftmp = open('tempfile', 'w')
ftmp.write(pickle.dumps(numpy.array(allDistances)))
ftmp.close()
print numpy.mean(numpy.array(allDistances), axis = 0)

#for total, dataSet in out:
#    matplotlib.pyplot.plot(dataSet)

#matplotlib.pyplot.show()
#matplotlib.pyplot.legend(totals)
