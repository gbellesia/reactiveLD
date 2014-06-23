# Pizza.py script to viz E Coli particles

import sys

if len(sys.argv) < 4:
    print "Error, input file required pizza.py -f viz.py [infile]"

d = dump(sys.argv[3])
d.map(1, "id", 2, "type", 3, "x", 4, "y", 5, "z")
g = gl(d)
g.arad([1, 2], [1.0, 1.0])
v = vcr(g)
