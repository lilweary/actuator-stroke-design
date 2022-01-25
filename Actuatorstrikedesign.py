import numpy as np
from EasyBeam import Beam2D
from scipy import optimize
import matplotlib.pyplot as plt

b = 10          # mm
h = 2          # mm
# p = -10     # N
l = 1000        # mm
E = 210000      # MPa
rho = 7.85e-9   # t/mm^3
I = b*h**3/12   # mm^4
A = b*h         # mm^2
nu = 0.3
nEl=100

xBeam=np.linspace(0,l,nEl+1)

#obstacle @(x0,y0)
x0=600
y0=50
m= 50 #margin 
y1=250 # target Coordinates (x1,y1)>(x0,y0)
# pivot point max Deformation at yy
xx = (x0+m)
yy = y0

linearly increasing
p=np.zeros(nEl+1)
p=(yy*6*E*I)/(2*xx**3)
p1=(y1*6*E*I)/((xx**2)*(3*l-xx))
y = np.zeros(nEl+1)
for i in range(nEl+1):
    # if xBeam[i] <=xx:    
        y[i]=((p*(xBeam[i])**2) / (6*E*I))*((3*xx - xBeam[i]))
    # else:
    #   y[i]=((p1*(xx)**2) / (6*E*I))*((3*xBeam[i]-xx))

lx=int((xx/l)*nEl)
# y cordinate
# linearly increasing
ytarget = np.zeros((nEl+1,))
ytarget [:lx+1] = y[:lx+1]
R=np.linspace(yy,y1, (nEl+1-lx))
ytarget [lx+1:] = R[1:]

# x cordianates
xtarget = np.linspace(0,l,nEl+1)

#target Shape
ShapeTarget = np.zeros(((nEl + 1) * 2,))
ShapeTarget[1::2] = ytarget

J=np.vstack([(xBeam),(ytarget)])
plt.plot(xtarget,J[1,:])
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
    Cantilever.Disp = [[1, [0.0, 0.0, 0.0]]]
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
            ii += 1

   dofControl = []
    nodesControl = range(Cantilever.nN)
    for i in range(Cantilever.nN):
        dofControl.append(nodesControl[i] * 3)
        dofControl.append(nodesControl[i] * 3 + 1)
    H = uMat[dofControl, :]
    uStroke = np.linalg.lstsq(H, ShapeTarget, rcond=None)[0]

    # Validate uStroke and calculate RMSE
    Cantilever = SetupFE()
    for jj in range(nAct):
        Cantilever.Disp.append(
            [nodesAct[jj], [uStroke[jj * 2], uStroke[jj * 2 + 1], 'f']]
        )
    Cantilever.StaticAnalysis()
    # Cantilever.PlotMesh(
    #     FontMag=1, NodeNumber=False, ElementNumber=False, Loads=False, BC=False
    # )
    Cantilever.PlotStress(stress='max', scale=20)
    eRMS = np.sqrt(
        np.sum((Cantilever.u[dofControl] - ShapeTarget) ** 2)
        / len(ShapeTarget)
    )
    FAct_dof = Cantilever.F[dofControl]
    FAct = np.zeros(len(nodesAct),)
    for i in range(len(nodesAct)):
        FAct[i] = np.sqrt(FAct_dof[i*2]**2+FAct_dof[i*2+1]**2)
    FActMax = np.max(FAct)
    return eRMS, uStroke, FActMax, np.max(Cantilever.sigmaMax)


# definition of actuators
nodesAct = [[]] * 4
eRMS = [[]] * 4
uStroke = [[]] * 4
F = [[]] * 4
sigmaMax = [[]]*4
nodesAct[0] = np.array(range(101, 102, 1)).tolist()
#nodesAct[3] = np.array(range(1, nEl, 1)).tolist()
nodesAct[1] = np.array(range(65, 102, 36)).tolist()
nodesAct[2] = np.array(range(65,  102, 18)).tolist()
nodesAct[3] = np.array(range(56, 102, 15)).tolist()

# Parameter study
nActList =[]
for i in range(len(nodesAct)):
    nActList.append(len(nodesAct[i]))
    eRMS[i], uStroke[i], F[i], sigmaMax[i] = CalcStroke(nodesAct[i])
    print(eRMS[i])

def adjust_spines(ax, spines, color):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(
                ('outward', 24)
            )
            if loc == 'bottom':
                spine.set_color('black')
            else:
                spine.set_color(color)

fig = plt.figure()
ax = fig.add_subplot()
ax2 = ax.twinx()
ax3 = ax.twinx()
ax.plot(nActList, eRMS, ".-", color="tab:blue")
ax2.plot(nActList, sigmaMax, ".-", color="tab:orange")
ax3.plot(nActList, F, ".-", color="tab:green")

plt.xticks(nActList)

ax.yaxis.set_label_coords(-0.1, 1.035)
ax2.yaxis.set_label_coords(1.1, 1.035)
ax3.yaxis.set_label_coords(1.5, 1.035)

ax.set_xlabel("number of actuators")
ax.set_ylabel(
    "root mean square\nerror $\\varepsilon_{\\mathrm{RMS}}$",
    rotation='horizontal',
    horizontalalignment='right',
    verticalalignment='baseline',
)
ax2.set_ylabel(
    "maximum stress\n$\\sigma_{\\max}$ [MPa]",
    rotation='horizontal',
    horizontalalignment='left',
    verticalalignment='baseline',
)
ax3.set_ylabel(
    "maximum actuator\nforce $F_{\mathrm{act},\\max}$ [N]",
    rotation='horizontal',
    horizontalalignment='left',
    verticalalignment='baseline',
)
ax.set_xlim([np.min(nActList), np.max(nActList)])
ax.set_ylim([0, np.max(eRMS)])
ax2.set_ylim([0, np.max(sigmaMax)])
ax3.set_ylim([0, np.max(F)])
ax.tick_params(direction='in')
ax2.tick_params(direction='in')
ax3.tick_params(direction='in')
adjust_spines(ax, ['bottom', 'left'], 'tab:blue')
adjust_spines(ax2, ['right'], 'tab:orange')
adjust_spines(ax3, ['right'], 'tab:green')
ax.tick_params(axis='y', colors='tab:blue')
ax2.tick_params(axis='y', colors='tab:orange')
ax3.tick_params(axis='y', colors='tab:green')
ax.yaxis.label.set_color('tab:blue')
ax2.yaxis.label.set_color('tab:orange')
ax3.yaxis.label.set_color('tab:green')
ax3.spines["right"].set_position(("axes", 1.5))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['left'].set_visible(False)
