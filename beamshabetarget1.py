# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 11:19:44 2021

@author: raraJ
"""
import numpy as np
from EasyBeam import Beam2D
from scipy import optimize
import matplotlib.pyplot as plt

b = 1        # mm
h = 2          # mm
# F = -100        # N
l = 1000        # mm
E = 210000      # MPa
rho = 7.85e-9   # t/mm^3
I = b*h**3/12   # mm^4
A = b*h         # mm^2
nu = 0.3
nEl=100

xBeam=np.linspace(0,l,nEl+1)
#box with corner located at
xx = 500                           ## desired x position
yy = 70                            ##pivot point
dd = 300                           ##desired y postion 
lx = int((xx/l)*nEl)               ## x elements
ly = int(((l-xx)/l)*nEl)           ## y Elements

# Positive slope line
y = np.zeros(ly+1)
y[:] = np.linspace(yy-5,dd,ly+1)

# # y cordinates
ytarget = np.zeros((nEl+1,))
ytarget [:lx] = xBeam[:lx] *(yy-5)/xx
ytarget [lx:] = y

# x cordianates
xtarget = np.linspace(0,l,nEl+1)
xc = (( 0.5 * ( dd - yy-5 )**2 - (y-((dd + yy-5)/2))**2)**0.5
                                                  + (xx-((dd-yy-5)/2)))
xtarget[lx:]=xc

ShapeTarget = np.zeros(((nEl + 1) * 2,))
ShapeTarget[1::2] = ytarget
ShapeTarget[0::2] = xtarget
J=np.vstack([(xtarget),(ytarget)])



plt.plot(J[0,:],J[1,:])
plt.show()

def SetupFE():
    Cantilever = Beam2D()
    Cantilever.SizingVariables = [["h", "b"]]
    Cantilever.stiffMatType = "Euler-Bernoulli"

    Cantilever.Nodes = [[]] * (nEl + 1)
    for i in range(nEl + 1):
        Cantilever.Nodes[i] = [l * i / nEl, 0.0]

    Cantilever.El = [[]] * (nEl)
    for i in range(nEl):
        Cantilever.El[i] = [i + 1, i + 2]

    Cantilever.Properties = [["Prop1", rho, E, nu, 1, h, b]]
    Cantilever.PropID = ["Prop1"] * nEl
    Cantilever.Disp = [[1, [0.0, 0.0, 0.0]],]
    return Cantilever


def CalcStroke(nodesAct):
    nAct = len(nodesAct)
    ii = 0
    uMat = np.zeros(((nEl + 1) * 3, nAct * 2))
    for i in range(nAct):
        for j in range(2):
            if j == 0:
                loading = [1.0, 0, ""]
            elif j == 1:
                loading = [0, 1.0, ""]
            x = [0] * (len(nodesAct) * 2)
            x[ii] = 1.0
            Cantilever = SetupFE()
            for jj in range(nAct):
                Cantilever.Disp.append([nodesAct[jj], [x[jj * 2], x[jj * 2 + 1], "f"]])
            Cantilever.StaticAnalysis()
            uMat[:, ii] = Cantilever.u
            # Cantilever.PlotDisplacement(component='y', scale=20)
            ii += 1

    dofControl = []
    nodesControl = range(Cantilever.nN)
    for i in range(Cantilever.nN):
        dofControl.append(nodesControl[i] * 3)
        dofControl.append(nodesControl[i] * 3 + 1)
    H = uMat[dofControl, :]

    res = optimize.lsq_linear(H, ShapeTarget)
    uStroke = res.x

    # Validate uStroke and calculate RMSE
    Cantilever = SetupFE()
    for jj in range(nAct):
        Cantilever.Disp.append(
            [nodesAct[jj], [uStroke[jj * 2], uStroke[jj * 2 + 1], "f"]]
        )
    Cantilever.StaticAnalysis()
    # only x and y (no rotation!)
    #Cantilever.PlotDisplacement(component="y", scale=20)
    Cantilever.PlotStress(stress='max', scale=20)
    eRMS = np.sqrt(
        np.sum((Cantilever.u[dofControl] - ShapeTarget) ** 2) / len(ShapeTarget)
    )
    return eRMS, uStroke, Cantilever.F

# definition of actuators
nodesAct = [[]] * 4
eRMS = [[]] * 4
uStroke = [[]] * 4
F = [[]] * 4
nodesAct[0] = np.array(range(101, 102, 1)).tolist()
nodesAct[1] = np.array(range(60, 102, 41)).tolist()
nodesAct[2] = np.array(range(51, 102, 25)).tolist()
nodesAct[3] = np.array(range(1, 102, 1)).tolist()
# nodesAct[3] = np.array(range(41, 102, 20)).tolist()
# nodesAct[6] = np.array(range(2, 52, 2)).tolist()


# Parameter study
for i in range(len(nodesAct)):
    eRMS[i], uStroke[i], F[i] = CalcStroke(nodesAct[i])
    print(eRMS[i])