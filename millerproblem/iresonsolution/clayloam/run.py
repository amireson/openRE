import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from Richards import run_RE

# Soil properties:
#import sys
#sys.path.insert(0,'..')
from setup import *

z=np.arange(dz,zN,dz)
n=len(z)

# Initial condition:
psi0=z-zN

t=np.arange(0,tN+dt,dt)
nt=len(t)
BC_T=psiT+np.zeros(nt)
BC_B=psiB+np.zeros(nt)
dum,WB,runtime=run_RE(dt,t[:2],dz,zN,n,psi0,BC_T,BC_B,pars)
psi,WB,runtime=run_RE(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)

z=np.hstack([0,z,zN])
#for i in range(nt):
#    pl.plot(np.hstack([psiT,psi[i,:],psiB]),z)
#    
#pl.ylim(10,0)
#pl.show()
psi=np.hstack([psiT,psi[-1,:],psiB])

# Save data:
fname='output/psi.csv'
f=open(fname,'w')
f.write('z, psi\n')


for i,j in zip (z,psi): f.write('%.4f, %.4f\n'%(i,j))
f.close()

f=open('output/runtime.csv','w')
f.write('%.3f\n'%runtime)
f.close()

fname='output/mb.csv'
WB.to_csv(fname,float_format='%12.8f')
