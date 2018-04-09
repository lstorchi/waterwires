import openbabel
import argparse
import pybel
import math
import glob 
import sys
import os

from operator import itemgetter

###############################################################################

def if_exist_rm (fname):
    if os.path.exists(fname):
        os.remove(fname)

###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-z","--minz", help="minimum difference along z", type=float,\
                default=5.0)
parser.add_argument("-x","--minx", help="minimum difference along x", type=float,\
                default=5.0)
parser.add_argument("-y","--miny", help="minimum difference along y", type=float,\
                default=5.0)
parser.add_argument("-d","--mind", help="minimum distance", type=float,\
                default=5.0)


args = parser.parse_args()

minx = args.minx
miny = args.miny
minz = args.minz
mind = args.mind

os.chdir("./")

if_exist_rm ("final.txt")
fp = open("final.txt", "w")
wires_centroids = {}

for file in glob.glob("*.pdb"):
    print(file)
    mol = list(pybel.readfile("pdb", file))
    basename = file[:-4]

    if len(mol) != 1:
        print "Each file should contain a single molecule set"
        exit(1)

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
            diffx = math.fabs(sort_acoordslist[i][0] - \
                    lastatom[0])
            diffy = math.fabs(sort_acoordslist[i][1] - \
                    lastatom[1])
            diffz = math.fabs(sort_acoordslist[i][2] - \
                    lastatom[2])
 
            dist = math.sqrt(diffx**2 + diffy**2 + diffz**2)
            wirecompleted = True
            if (diffz <= minz) and (diffx <= minx) and (diffy <= miny) \
                    and (dist <= mind):
                    wire.append(i)
                    wirecompleted = False

            if wirecompleted:
                wires.append(wire)
                wire = []
                lastidx = i
                lastatom = sort_acoordslist[lastidx]
                wire.append(lastidx)
            else:
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
        tot = 0
        molidx = 1
        for wire in wires:
            wirelen = len(wire)
            print "   ", wirelen
            fp.write(str(wirelen) + "\n")
            tot += wirelen

            if not (wirelen in wires_centroids):
                wires_centroids[wirelen] = []

            obmol = openbabel.OBMol()
            lidx = 0
            xcentr = 0.0
            ycentr = 0.0
            zcentr = 0.0
            for idx in wire:
                x = sort_acoordslist[idx][0]
                y = sort_acoordslist[idx][1]
                z = sort_acoordslist[idx][2]

                xcentr += x
                ycentr += y
                zcentr += z

                newobatom = openbabel.OBAtom()
                newobatom.SetAtomicNum(8)
                newobatom.SetId(lidx)
                newobatom.SetIdx(lidx)
                newobatom.SetVector(x, y, z)
                obmol.AddAtom(newobatom)
                lidx = lidx + 1

            xcentr = xcentr / float(wirelen)
            ycentr = ycentr / float(wirelen)
            zcentr = zcentr / float(wirelen)

            wires_centroids[wirelen].append((xcentr, ycentr, zcentr))

            for idx in range(1,len(wire)):
                obmol.AddBond(idx-1, idx, 1)

            pybelmol = pybel.Molecule(obmol)
            pybelmol.write("mol2", basename + "_" + str(molidx)+".mol2", True)
            molidx = molidx + 1

        if tot != len(sort_acoordslist):
            print "error in total number"

fp.close()

for num in wires_centroids:
    fname = str(num) + ".xyz"
    if_exist_rm (fname)

    fp = open(fname, "w")
    fp.write(str(len(wires_centroids[num])) + "\n")
    fp.write(str(num) + "\n")

    for c in wires_centroids[num]:
        fp.write("H %12.5f %12.5f %12.5f\n"%(c[0], c[1], c[2]))
    fp.close()
