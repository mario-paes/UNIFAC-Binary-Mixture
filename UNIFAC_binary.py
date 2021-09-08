# -*- coding: utf-8 -*-


from thermo.unifac import UFIP, UFSG, UNIFAC

#Calculating a sample example

GE = UNIFAC.from_subgroups(chemgroups=[{1:4, 3:2}, {50:1}], T=56.1+273.15, xs=[0.4, 0.6], version=0, interaction_data=UFIP, subgroups=UFSG)
print ( GE.gammas())

# this means that the mixture is 
# 0.4 mole fraction of species 1
# species 1 consists of two functional groups
# {1:4, 3:2} means subgroup 1 (CH3) with a count of 4
# and subgroup 3 (CH) with a count of 2
# so this is 2,2-Dimethylbutane

# the 2nd compound has only one functional group
# 50 (CHCl3) at count 1
# this is Chloroform

T_K = 273.15 + 56.1
import numpy as np
xv = np.linspace(0,1,30)
g1v = np.empty_like(xv)
g2v = np.empty_like(xv)
cc=0
for x in xv:
    
    GE = UNIFAC.from_subgroups(chemgroups=[{1:4, 3:2}, {50:1}], T=56.1+273.15, xs=[x, 1-x], version=0, interaction_data=UFIP, subgroups=UFSG)
    g1v[cc]=GE.gammas()[0]
    g2v[cc]=GE.gammas()[1]
    cc+=1

import pylab as pl

fig1=pl.figure()
pl.plot(xv,g1v,color='blue')
pl.plot(xv,g2v,color='green')
pl.xlabel('x1')
pl.ylabel('activity coef')
pl.title('Activity Coef of 2,2-Dimethylbutane (1) Chloroform (2) System from UNIFAC via Thermo.py')
pl.show(fig1)