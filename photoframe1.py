# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 23:45:30 2021

@author: raraJ
"""
from EasyBeam import Beam2D
import numpy as np
import matplotlib.pyplot as plt
b = 10          # mm
h = 20          # mm
F = -100        # N
l = 1000        # mm
E = 210000      # MPa
rho = 7.85e-9   # t/mm^3
I = b*h**3/12   # mm^4
A = b*h         # mm^2
nu = 0.3

# Initialisiern des Problems
Frame = Beam2D()

# Knoten [mm]
Frame.Nodes = [[0,0],[600,0],[600,190],[600,400]]
# Frame.Nodes = [[0,0],[600,0],[800,0],[1000,0]]

# Elemente: verbindet die Knoten
Frame.El = [[]]*(3)
for i in range(3):
    Frame.El[i] = [i+1, i+2]


Frame . Properties = [[ 'Prop1 ', rho , E, nu , 1, h, b]]
Frame . PropID = ['Prop1 '] * (5)

Frame.Disp = [[1, [ 0 , 0, 0]],[2, [ 'f' , 1, 0]],[4, [ 0, 0 , 0]]]              
Frame.Load = [[2, [  0, 0, 0]],[3, [6000, 0 , 0]],[4, [ 0, 0 , 0]]]

# Frame.Disp = [[1, [  0, 0, 0]],[2, [  0, 1, 0]],[3, [  -200,190, 0]],
#             [4, [  -400, 400, 0]]]
              
# Frame.Load = [[3, [1000, 0, 0]]]

Frame.nSeg = 20

# Lösen
Frame.StaticAnalysis()
# Frame.EigenvalueAnalysis(nEig=len(Frame.DoF))

# Grafische Darstellung
Frame.PlotMesh(FontMag=2)
Frame.PlotStress(stress="all", scale=100)
Frame.PlotDisplacement(component="all", scale=100)
# Frame.PlotMode(scale=30)


