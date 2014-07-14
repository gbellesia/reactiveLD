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
    print "Include an output file -- ./irr.py outfilename"
    exit(-1)

outputfile = sys.argv[1]

distances = []

for i in range(0, 100000):
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
	       "simType" : "periodicBox",
               "dt" : 0.001,
               "steps" : 18000,
               "printSteps" : 1,
               "monatomic" : [],
               "binary" : [{ "A" : 1, "B" : 1, "k" : 2.010619298 }],
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
    lll = f.readlines()
    f.close()
    if (len(lll[len(lll)-1].split()) == 4 and len(lll[len(lll)-2].split()) == 4):
        i1, x, y, z = lll[len(lll)-1].split()
        i1, x2, y2, z2 = lll[len(lll)-2].split()
        x = float(x)
        y = float(y)
        z = float(z)
        x2 = float(x2)
        y2 = float(y2)
        z2 = float(z2)
        dx = (x - x2)
        dy = (y - y2)
        dz = (z - z2)
	#print math.sqrt(dx*dx + dy*dy + dz*dz)
        distances.append(math.sqrt(dx*dx + dy*dy + dz*dz))

        #print x, y, z, x2, y2, z2, distances[-1]

    print i
    os.remove(outfile)
    os.remove(infile.name)

f = open(outputfile, 'w')
f.write("\n".join(map(repr, distances)))

#matplotlib.pyplot.xlabel('Distances')
#matplotlib.pyplot.hist(distances, 100)
#matplotlib.pyplot.show()
