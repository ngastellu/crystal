#!/usr/bin/env python

import numpy as np
from math import *
from atom import Atom
from crystal import Crystal
from scipy.interpolate import interp1d, interp2d


input_file = raw_input("Input file: ")

crystem = Crystal(input_file)
os.chdir('/'.join(input_file.split('/')[:-1]))


