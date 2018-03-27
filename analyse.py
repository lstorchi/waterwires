import argparse
import pybel
import glob 
import sys
import os

import os.path

os.chdir("./")

for file in glob.glob("*.pdb"):
    print(file)
