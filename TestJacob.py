# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 21:26:25 2022

@author: Jacob Mehr'
"""

import boreholeStress as zoback
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri
s1=20
s2=18
s3=17
pp=12
# alpha,beta,gamma = 0,90,0 Normal Faulting
# alpha,beta,gamma = 90,0,90 Strike Slip Faulting
# alpha,beta,gamma = 90,0,0 Reverse Faulting
alpha=0
beta=90
gamma=0
# Delta = angle away from North
# phi = deviation from vertical 
delta=90
phi=90
v= .4 #Poisson Ratio (Usually .25)
to=1.923
thetacoh = 30
TH=0
viz = False
#Pmmin,Pmax = zoback.boreholeStress(s1, s2, s3, pp, alpha, beta, gamma, delta, phi, v,to,thetacoh,TH, viz)
#print(Pmmin,Pmax)

window = []
ygraph = []
xgraph = []
for i in range(0,95,5):
    for j in range (0,365,5):
        print(i,j)
        Pmmin,Pmax = zoback.boreholeStress(s1, s2, s3, pp, alpha, beta, gamma, j,i , v,to,thetacoh,TH,viz)
        delta2 = np.abs(Pmmin-Pmax)
        jpi = math.radians(j)
        ipi = math.radians(i)
        x = np.cos(jpi)
        y =np.sin(jpi)
        r = np.sin(ipi)
        xplot = x*r
        yplot = y*r
       
        window.append(delta2)
        xgraph.append(xplot)
        ygraph.append(yplot)

        #plt.show()
triang = tri.Triangulation(xgraph, ygraph)
fig1, ax1 = plt.subplots()
tpc = ax1.tripcolor(triang, window, shading='flat')
fig1.colorbar(tpc)
plt.grid()