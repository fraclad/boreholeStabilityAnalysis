import numpy as np
from numpy import sin, cos, tan
import matplotlib.pyplot as plt
import math
plt.rcParams["font.family"] = "Helvetica"

def degToRad(ang, reverse = False):
    res = (np.pi/180)*ang
    if reverse == True:
        res = ang * 180 / np.pi
    return res

def boreholeStress(s1, s2, s3, pp, alpha, beta, gamma, delta, phi, v,to,thetacoh,TS, viz = True):
    thetacoh = thetacoh
    Ss = np.eye(3)
    
    Ss[0][0] = s1
    Ss[1][1] = s2
    Ss[2][2] = s3
    Pp = pp
    
    Rs = np.zeros((3,3))
    
    alp = degToRad(alpha)
    bet = degToRad(beta)
    gam = degToRad(gamma)
    
    Rs[0][0] = cos(alp)*cos(bet)
    Rs[0][1] = sin(alp)*cos(bet)
    Rs[0][2] = -1 * sin(bet)
    Rs[1][0] = cos(alp)*sin(bet)*sin(gam) - sin(alp)*cos(gam)
    Rs[1][1] = sin(alp)*sin(bet)*sin(gam) + cos(alp)*cos(gam)
    Rs[1][2] = cos(bet)*sin(gam)
    Rs[2][0] = cos(alp)*sin(bet)*cos(gam) + sin(alp)*sin(gam)
    Rs[2][1] = sin(alp)*sin(bet)*cos(gam) - cos(alp)*sin(gam)
    Rs[2][2] = cos(bet)*cos(gam)
    
    Rb = np.zeros((3,3))
    
    de = degToRad(delta)
    phi = degToRad(phi)
    
    Rb[0][0] = -cos(de)*cos(phi)
    Rb[0][1] = -sin(de)*cos(phi)
    Rb[0][2] = sin(phi)
    Rb[1][0] = sin(de)
    Rb[1][1] = -cos(de)
    Rb[1][2] = 0
    Rb[2][0] = cos(de)*sin(phi)
    Rb[2][1] = sin(de)*sin(phi)
    Rb[2][2] = cos(phi) 
    
    S = Rb @ Rs.transpose() @ Ss @ Rs @ Rb.transpose()
    Seff = S.copy()
    Seff[0][0] = Seff[0][0] - Pp
    Seff[1][1] = Seff[1][1] - Pp
    Seff[2][2] = Seff[2][2] - Pp
    ###
    v = v
    theta = np.linspace(0,2*np.pi)
    i = pp - (pp*.5)
    thetacoh = degToRad(thetacoh)
    while i < pp*3:
        Sigmas=[0]
        dp = i - pp
        Szz = Seff[2][2] - 2*v*(Seff[0][0] - Seff[1][1])*cos(2*theta) - 4*v*Seff[0][1]*sin(2*theta)
        Stt = Seff[0][0] + Seff[1][1] - 2*(Seff[0][0] - Seff[1][1])*cos(2*theta)  - 4*Seff[0][1]*sin(2*theta) - dp
        Stz = 2*(Seff[1][2]*cos(theta) - Seff[0][2]*sin(theta))
        Srr = np.abs(dp)
        Stmax = 0.5*(Szz + Stt + ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Stmin = 0.5*(Szz + Stt - ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Srr = np.repeat(dp, len(Stmin))
        StmaxNorm = Stmax 
        StminNorm = Stmin
        SrrNorm = Srr
        Stmaxmax = np.amax(StmaxNorm)
        Stminmin = np.amax(StminNorm)
        Srrmin = np.amin(Srr)
        Sigmas.append(Srrmin)
        Sigmas.append(Stmaxmax)
        Sigmas.append(Stminmin)
        Sigmas.pop(0)
        Sigmassort = np.sort(Sigmas,axis=None)
        sider1 = to*np.cos(thetacoh)
        sider2 = .5*np.sin(thetacoh)*(Sigmassort[2]+Sigmassort[0])
        sider = sider1+sider2
        sideleft = .5*(Sigmassort[2]-Sigmassort[0])
        sider = np.round(sider,2)
        sideleft = np.round(sideleft,2)
        diff = np.abs(sider-sideleft)
        if diff<.1 :
            Pmin = i
            i = pp*4
            print(Pmin)
            break
        else:
            i+=.1
    j = Pmin + (Pmin*.25)
    while j < Pmin*4:
        Sigmas=[0]
        dp = j - pp
        Szz = Seff[2][2] - 2*v*(Seff[0][0] - Seff[1][1])*cos(2*theta) - 4*v*Seff[0][1]*sin(2*theta)
        Stt = Seff[0][0] + Seff[1][1] - 2*(Seff[0][0] - Seff[1][1])*cos(2*theta)  - 4*Seff[0][1]*sin(2*theta) - dp
        Stz = 2*(Seff[1][2]*cos(theta) - Seff[0][2]*sin(theta))
        Srr = dp
        Stmax = 0.5*(Szz + Stt + ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Stmin = 0.5*(Szz + Stt - ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Srr = np.repeat(dp, len(Stmin))
        StmaxNorm = Stmax 
        StminNorm = Stmin
        SrrNorm = np.abs(Srr)
        Stmaxmax = np.amin(StmaxNorm) 
        Stminmin = np.amin(StminNorm) 
        Srrmin = np.amin(Srr) 
        Sigmas.append(np.abs(Srrmin))
        Sigmas.append(np.abs(Stmaxmax))
        Sigmas.append(np.abs(Stminmin))
        Sigmas.pop(0)
        Sigmassort = np.sort(Sigmas,axis=None)
        sider1 = 3*Sigmassort[0]
        sider2 = -1*Sigmassort[2]  + TS - pp 
        sider = sider1 + sider2
        sider = np.abs(np.round(sider,4))
        j = np.round(j,4)
        diff = np.abs(j-sider)
        if diff<=.3:
            Pmax = j
            print(Pmax)
            break
        else:
            j+=.1

  ###

    if viz == True:
       plt.plot(degToRad(theta[:51], reverse = True), StmaxNorm[:51], 
                 label = "$\sigma_{t \ max}$", color = "red")
       plt.plot(degToRad(theta[:51], reverse = True), StminNorm[:51], 
                 label = "$\sigma_{t \ min}$", color = "blue")
       plt.plot(degToRad(theta[:51], reverse = True), SrrNorm[:51], 
                 label = "$\sigma_{rr}$", color = "green")
       plt.xlabel("borehole angle ($\degree$)")
       plt.ylabel("Effective stress (normalized with $\sigma_1$)")
       plt.legend()
       plt.tight_layout()
       plt.show()
    
    resultDict = {"Szz":Szz,
                  "stt":Stt,
                  "Stz":Stz,
                  "Srr":dp,
                  "Stmax":Stmax,
                  "Stmin":Stmin
        }

    return Pmin,Pmax

def MohrCoulombSlope(mu):
    """
    

    Parameters
    ----------
    mu : float
        frcition coeff.

    Returns
    -------
    res : float
        the slope in MC failure criterion.

    """
    res = ((mu**2 + 1)**0.5 + mu)**2
    return res

def getMaxStrength(s1,s3,slope):
    """
    

    Parameters
    ----------
    s1 : float
        maximum effective stress.
    s3 : float
        minimum effective stress.
    slope : object
        function to calculate the slope in MC criterion 

    Returns
    -------
    res : float
        max strength allowable with given stress configs.

    """
    res = s1 - slope*s3
    return res














''''
  while j < pp*3:
        Sigmas=[0]
        dp = j - pp
        Szz = Seff[2][2] - 2*v*(Seff[0][0] - Seff[1][1])*cos(2*theta) - 4*v*Seff[0][1]*sin(2*theta)
        Stt = Seff[0][0] + Seff[1][1] - 2*(Seff[0][0] - Seff[1][1])*cos(2*theta)  - 4*Seff[0][1]*sin(2*theta) - dp
        Stz = 2*(Seff[1][2]*cos(theta) - Seff[0][2]*sin(theta))
        Srr = dp
        Stmax = 0.5*(Szz + Stt + ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Stmin = 0.5*(Szz + Stt - ((Szz - Stt)**2 + 4*(Stz)**2)**0.5)
        Srr = np.repeat(dp, len(Stmin))
        StmaxNorm = Stmax 
        StminNorm = Stmin
        SrrNorm = np.abs(Srr)
        Stmaxmax = np.amin(StmaxNorm) 
        Stminmin = np.amin(StminNorm) 
        Srrmin = np.amin(Srr) 
        Sigmas.append(np.abs(Srrmin))
        Sigmas.append(np.abs(Stmaxmax))
        Sigmas.append(np.abs(Stminmin))
        Sigmas.pop(0)
        Sigmassort = np.sort(Sigmas,axis=None)
        sider1 = 3*Sigmassort[0]
        sider2 = -1*Sigmassort[2]  + TS - pp 
        sider = sider1 + sider2
        sider = np.abs(np.round(sider,4))
        j = np.round(j,4)
        diff = np.abs(j-sider)
        if diff<.5:
            Pmax = j
            print(Pmax)
            break
        else:
            j+=.1
'''


























