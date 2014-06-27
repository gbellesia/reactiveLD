#!/usr/bin/env python

import json
import math
import os
import random
import shlex
import subprocess
import sys
import tempfile
import random
import matplotlib.pyplot

Zr = 0

if len(sys.argv) < 2:
    print "Include an output file -- ./rev.py outfilename"
    exit(-1)

outputfile = sys.argv[1]

distances = []

for i in range(0, 30000):
    if os.path.isdir("./"):
        tmpdir = "./"
    else:
        tmpdir = "/tmp"

    infile = tempfile.NamedTemporaryFile(dir = tmpdir, delete = False)
    outfile = tempfile.mktemp(dir = tmpdir)

    theta = math.pi * random.random()
    phi = 2 * math.pi * random.random()

    config = { "X" : 10000,
               "Y" : 10000,
               "Z" : 10000,
               "dt" : 1.0,
               "steps" : 2,
               "printSteps" : 1,
               "monatomic" : [{ "C" : 3, "A" : 1, "B" : 2, "k" : ??????? }],
	       "binary" : [{ "A" : 1, "B" : 2, "C" : 3,"k" : 2.010619298 }],
               "atoms" : [{"x" : 5000,
                           "y" : 5000,
                           "z" : 5000,
                           "type" : 1}],
               "types" : [{ "type" : 1, "radius" : 0.5, "D" : 0.5 },
                          { "type" : 2, "radius" : 0.5, "D" : 0.5 },
                          { "type" : 3, "radius" : 0.5, "D" : 0.5 }],}
    
    atom2 = { "x" : 5000 + 1.00000001 * math.sin(theta) * math.cos(phi),
              "y" : 5000 + 1.00000001 * math.sin(theta) * math.sin(phi),
              "z" : 5000 + 1.00000001 * math.cos(theta),
              "type" : 1 }

    config["atoms"].append(atom2)

    infile.write(json.dumps(config))
    infile.close()

    #print './bd_run {0} {1}'.format(infile.name, outfile)
    handle = subprocess.Popen(shlex.split('./bd_run {0} {1} {2}'.format(infile.name, outfile, random.randint(0, 1000000))), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    handle.communicate()

    infile.close()

    f = open(outfile, 'r')
    for j in range(0, 14):
        f.readline()
    count = int(f.readline().strip())
    for j in range(0, 5):
        f.readline()

    if count == 0:
        Zr += 1
    else:
        i1, t, x, y, z = f.readline().strip().split()
        i1, t, x2, y2, z2 = f.readline().strip().split()

        x = float(x)
        y = float(y)
        z = float(z)

        x2 = float(x2)
        y2 = float(y2)
        z2 = float(z2)

        distances.append(math.sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2) + (z - z2) * (z - z2)))

        #print x, y, z, x2, y2, z2, distances[-1]

    print i
    f.close()
    os.remove(outfile)
    os.remove(infile.name)

f = open(outputfile, 'w')
f.write("\n".join(map(repr, distances)))
f.close()

matplotlib.pyplot.xlabel('Distances')
matplotlib.pyplot.hist(distances, 100)
matplotlib.pyplot.show()
