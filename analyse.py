import argparse
import pybel
import glob 
import sys
import os

from operator import itemgetter

parser = argparse.ArgumentParser()
parser.add_argument("-x","--axis", help="axis to look for", type=str,\
                default="z")

args = parser.parse_args()

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
    for a in sort_acoordslist:
        print a
