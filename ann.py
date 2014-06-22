#!/usr/bin/env python

import matplotlib.pyplot
import re
import sys

if len(sys.argv) < 2:
    print "Correct usage: ./ann.py [outputfileofsimulation]"
    print "Where the simulation is: ./bd_run ann.in [outputfileofsimulation] 50003409"
    exit(-1)

h = open(sys.argv[1], 'r')

readNext = False
counts = []
for line in h.readlines():
    if readNext:
        s = line.split()
        if re.match("[0-9]+", s[0]):
            #print int(line)
            if int(s[1]) == 1:
                count += 1
        else:
            counts.append(count)
            readNext = False
    if re.match("ITEM: ATOMS", line):
        readNext = True
        count = 0
        
h.close()

matplotlib.pyplot.plot(counts)
matplotlib.pyplot.show()

