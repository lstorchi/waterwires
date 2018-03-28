import argparse
import pybel
import math
import glob 
import sys
import os

from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("-z","--minz", help="minimum difference along z", type=float,\
                default=10.0)
parser.add_argument("-x","--minx", help="minimum difference along x", type=float,\
                default=5.0)
parser.add_argument("-y","--miny", help="minimum difference along y", type=float,\
                default=5.0)


args = parser.parse_args()

minx = args.minx
miny = args.miny
minz = args.minz

os.chdir("./")

for file in glob.glob("*.pdb"):
    print(file)
    mol = list(pybel.readfile("pdb", file))

    if len(mol) != 1:
        print "Each file should contain a single molecule set"
        eit(1)

    acoordslist = []
    for a in mol[0].atoms:
        if a.atomicnum == 8:
            acoordslist.append(a.coords)

    sort_acoordslist = sorted(acoordslist, key=itemgetter(2))

    if (len(sort_acoordslist) > 1):
        lastidx = 0
        lastatom = sort_acoordslist[lastidx]

        wires = []
        wire = []
        wire.append(lastidx)
        for i in range(1,len(sort_acoordslist)):
            diffz = sort_acoordslist[i][2] - \
                    lastatom[2]
            wirecompleted = True
            if (diffz <= minz):
                diffx = math.fabs(sort_acoordslist[i][0] - \
                        lastatom[0])
                diffy = math.fabs(sort_acoordslist[i][1] - \
                        lastatom[1])

                if (diffx <= minx) and (diffy <= miny):
                    wire.append(i)
                    wirecompleted = False

            if wirecompleted:
                wires.append(wire)
                wire = []
                lastidx = i
                lastatom = sort_acoordslist[lastidx]

        if len(wires) == 0:
            wires.append(wire)
        else:
            lastwire = wires[-1]
            areequal = False
            if len(wire) == len(lastwire):
                valmatch = [i for i, j in zip(wire, lastwire) if i == j]
                if valmatch == len(wire):
                    areequal = True

            if not areequal:
                wires.append(wire)

        print "Num. of wires: ", len(wires)
        for wire in wires:
            print "   ", len(wire)
