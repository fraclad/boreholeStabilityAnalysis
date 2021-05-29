import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Helvetica"
from boreholeStress import *

# define angle domain and steps
inclination = np.linspace(0,90,10)
azimuth = np.linspace(0,360,37)

# create tile arrays
azimTiled = np.array([])
for i in inclination:
    add = np.tile(i, len(azimuth))
    azimTiled = np.concatenate((azimTiled, add))
incTiled = np.tile(azimuth, len(inclination))

# need to loop this
S = boreholeStress(3.12,3.12,1,0,0,65,0,45,30,0.25,0, viz = False)
S1 = max(S["Stmax"])
S3 = max(S["Stmin"])
C = getMaxStrength(S1, S3, MohrCoulombSlope(1))

