#!/usr/bin/env python

import re

h = open('ann.out.txt', 'r')

readNext = False
for line in h.readlines():
    if readNext:
        s = line.split()
        if re.match("[0-9]+", s[0]):
            #print int(line)
            if int(s[1]) == 1:
                count += 1
        else:
            print count
            readNext = False
    if re.match("ITEM: ATOMS", line):
        readNext = True
        count = 0
        
h.close()
