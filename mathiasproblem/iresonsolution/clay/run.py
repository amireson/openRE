import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

import sys
sys.path.insert(0,'..')
from Richards import run_RE

def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se


# Soil properties:
from setup import *

z=np.arange(dz,zN,dz)
n=len(z)

# Initial condition:
psi0=np.zeros(n)+psiB

t=np.arange(0,tN+dt,dt)
nt=len(t)
BC_T=psiT+np.zeros(nt)

to=[5,10,20,50,100]
tint=np.rint(t*24*60).astype('int')
ti=[np.where(tint==i)[0][0] for i in to]

BC_B=psiB+np.zeros(nt)
psi,WB,runtime=run_RE(dt,t[:2],dz,zN,n,psi0,BC_T,BC_B,pars)
psi,WB,runtime=run_RE(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)

theta=thetaFun(psi[ti,:],pars)
#pl.plot(z,theta.T)
#pl.show()

np.savetxt('theta.csv', np.vstack([z,theta]).T, delimiter=',',fmt='%10.4f')

#for i in range(nt):
#    pl.plot(np.hstack([psiT,psi[i,:],psiB]),z)
#    
#pl.ylim(10,0)
#pl.show()
#psi=np.hstack([psiT,psi[-1,:],psiB])
#
# Save data:
#fname='theta.csv'
#
#
#for i,j in zip (z,psi): f.write('%.4f, %.4f\n'%(i,j))
#f.close()
#
f=open('runtime.csv','w')
f.write('%.3f\n'%runtime)
f.close()
#
fname='mb.csv'
WB.to_csv(fname,float_format='%12.8f')
