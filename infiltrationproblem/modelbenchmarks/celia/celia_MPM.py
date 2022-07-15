import numpy as np
from numba import jit 
import pandas as pd
import matplotlib.pyplot as pl
from time import time

# MyNumba.py
from numba import types
from numba.typed import Dict

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

@jit(nopython=True)
def thetafun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

@jit(nopython=True)
def Cfun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

@jit(nopython=True)
def Kfun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

@jit(nopython=True)
def solverfun(R,C,Kmid,dt,dz,n):
    # Initialize arrays
    a=np.zeros(n)
    b=np.zeros(n)
    c=np.zeros(n)
    y=np.zeros(n)

    # Construct matrix
    a=Kmid[:-1]/dz
    b=-(Kmid[:-1]+Kmid[1:])/dz-C*dz/dt
    c=Kmid[1:]/dz

    # Lower boundary:
    b[0]=-Kmid[1]/dz-C[0]*dz/dt

    # Upper boundary:
    b[-1]=-Kmid[-2]/dz-C[-1]*dz/dt

    A=np.diag(a[1:],-1)+np.diag(b,0)+np.diag(c[:-1],1)

    # Construct RHS
    y[:]=R[:]

    # Boundary conditions - nothing to do
    #y[0]=y[0]
    #y[-1]=y[-1]

    # Solve:
    dell = np.linalg.solve(A, y)

    return dell

@jit(nopython=True)
def Rfun(psiiter,psiin,qT,KB,C,Kmid,dtheta,dt,dz,n):
    # This solves the Picard residual term:
    zero=np.array([0.])
    psigrid=np.hstack((zero,psiiter,zero))

    x1=dtheta/dt*dz
    x2=-(Kmid[1:]-Kmid[:-1])
    x3=-Kmid[1:]*(psigrid[2:]-psigrid[1:-1])
    x4=Kmid[:-1]*(psigrid[1:-1]-psigrid[:-2])

    # Lower boundary:
    x2[0]=-Kmid[1]
    x4[0]=KB

    # Upper boundary:
    x2[-1]=Kmid[-2]
    x3[-1]=qT

    R=x1+x2+x3+x4

    return R

@jit(nopython=True)
def iterfun(psiin,pars,qT,psiT,psiB,dt,dz,n):
    # psiin = psi^n
    # psiiter = psi^n+1,m
    # psiout = psi^n+1,m+1

    tolerance=1e-10
    maxcount=1000
    Rmax=1.
    zero=np.array([0.])

    # Initialize arrays
    psiiter=np.zeros(len(psiin))
    psiout=np.zeros(len(psiin))

    # Initial guess: psi_n+1^1 = psi_n
    psiiter[:]=psiin[:]

    count=0.
    while count <= 1 or (Rmax >= tolerance and count<= maxcount):
        # Get C,K:
        C=Cfun(psiiter,pars)
        K=Kfun(np.hstack((zero, psiiter, zero)),pars)
        Kmid=(K[1:]+K[:-1])/2.
        dtheta=thetafun(psiiter,pars)-thetafun(psiin,pars)
        # Get R
        R=Rfun(psiiter,psiin,qT,K[1],C,Kmid,dtheta,dt,dz,n)
        # Solve for del
        dell=solverfun(R,C,Kmid,dt,dz,n)
        # Update psi estimates at different iteration levels
        psiout[:]=psiiter[:]+dell[:]

#        # Plot to check convergence:
#        pl.plot(z,psiiter)
#        pl.plot(z,psiout)
#        pl.show()

        #err=psiout-psiiter
        psiiter[:]=psiout[:]
        Rmax=np.abs(np.max(R))
        count+=1

    #print('Iteration count = %d'%(count-1))

    return psiout

@jit(nopython=True)
def massbal(psi,qT,psiT,psiB,pars,n,dt,dz):

    # Initial storage:
    theta=thetafun(psi.reshape(-1),pars)
    theta=np.reshape(theta,psi.shape)
    S=np.sum(theta*dz,1)
    S0=S[0]
    SN=S[-1]

    # Inflow:
    Kin=Kfun(psi[:,0],pars)
    QIN=-Kin
    QIN[0]=0.
    QINsum=np.sum(QIN)*dt

    # Outflow:
    QOUT=qT
#     Kout=(Kfun(psi[:,-1],pars)+Kfun(psiT,pars))/2.
#     QOUT=-Kout*((psiT-psi[:,-1])/dz+1.)
    QOUT[0]=0.
    QOUTsum=np.sum(QOUT)*dt

    # Balance:
    dS=SN-S0
    dQ=QINsum-QOUTsum
    err=dS/dQ

#    print('Delta storage = %.3f'%dS)
#    print('Flow at base (+ve upwards) = %.3f'%QINsum)
#    print('Flow at surface (+ve upwards) = %.3f'%QOUTsum)
#    print('Delta flow = %.3f'%dQ)
#    print('Error metric = %.3f'%err)
    return QIN,QOUT,S,err

@jit(nopython=True)
def ModelRun(dt,dz,n,nt,psi,qT,psiB,psiT,pars):
    # Solve:
    for j in range(1,nt):
        psi[j,:]=iterfun(psi[j-1,:],pars,qT[j],psiT,psiB,dt,dz,n)

    QIN,QOUT,S,err=massbal(psi,qT,psiT,psiB,pars,n,dt,dz)

    return psi,QIN,QOUT,S,err
