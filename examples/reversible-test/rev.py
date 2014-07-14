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
	       "simType" : "periodicBox",
               "dt" : 0.01,
               "steps" : 2001,
               "printSteps" : 2000,
               "monatomic" : [{"A" : 3, "B" : 1, "C" : 2, "k" : 0.1}],
               "binary" : [{ "A" : 1, "B" : 2, "C" : 3, "k" : 2.010619298 }],
               "atoms" : [{"x" : 5000,
                           "y" : 5000,
                           "z" : 5000,
                           "type" : 1}],
               "types" : [{ "type" : 1, "radius" : 0.5, "D" : 0.0 },
                          { "type" : 2, "radius" : 0.5, "D" : 1.0 },
                          { "type" : 3, "radius" : 0.5, "D" : 0.0 }],}
    
    atom2 = { "x" : 5000 + 1.00000001 * math.sin(theta) * math.cos(phi),
              "y" : 5000 + 1.00000001 * math.sin(theta) * math.sin(phi),
              "z" : 5000 + 1.00000001 * math.cos(theta),
              "type" : 2 }

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

    # We set up the input to only print out two timesteps: 0 and 2000
    # The intput to get this is steps = 2001, printSteps = 2000
    # The code prints at time 0 by default, and then every time the stepCounter % printSteps == 0
    # For the standard c loop: for(int stepCounter = 0; stepCounter < steps; stepCounter++)
    # stepCounter goes to 2000 inside the loop *only* if steps is 2001
    # Fun nuances of code. Yeeeeeeah.
    #
    # We're looking for an output file that looks like:
    # 2
    # go spurs
    # 1 5000 5000 5000
    # 1 rndx rndy rndz
    # (either 1 or 2 here)
    # go spurs
    # (either 1 or two lines of stuff here)

    # THe best way to do this is look for the second 'go spurs'
    # When we find it, read the next two lines

    gospurs = 0
    readNext = False
    dataLines = []
    for line in lll:
        if readNext:
            dataLines.append(line)

        if line.strip() == "go spurs":
            gospurs += 1

        if gospurs == 2:
            readNext = True

    #print dataLines

    #print len(dataLines)

    if len(dataLines) > 1:
        i1, x, y, z = dataLines[0].split()
        i2, x2, y2, z2 = dataLines[1].split()
        x = float(x)
        y = float(y)
        z = float(z)
        x2 = float(x2)
        y2 = float(y2)
        z2 = float(z2)
        dx = (x - x2)
        dy = (y - y2)
        dz = (z - z2)
        distances.append(math.sqrt(dx*dx + dy*dy + dz*dz))
    
    print i
    os.remove(outfile)
    os.remove(infile.name)

f = open(outputfile, 'w')
f.write("\n".join(map(repr, distances)))
