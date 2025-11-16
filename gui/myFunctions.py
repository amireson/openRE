import pandas as pd
import numpy as np

from Richards_FD import run_RE as run_RE_FD
from Richards_FP import run_RE as run_RE_FP
from Richards_ZF import run_RE as run_RE_ZF

def HydProps(psi_min,psi_max,thetaR, thetaS,alpha,n,KS):
    global pars
    psi=np.linspace(psi_min,psi_max,301)
    pars={}
    pars['thetaR']=thetaR
    pars['thetaS']=thetaS
    pars['alpha']=alpha
    pars['n']=n
    pars['m']=1-1/n
    pars['Ks']=KS
    pars['neta']=0.5
    pars['Ss']=1e-6

    theta=thetaFun(psi,pars)
    C=CFun(psi,pars)
    K=KFun(psi,pars)

    return psi,theta,C,K

def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

def CFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

def KFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

def runRE(RunTime,TimeStep,SoilDepth,SpaceStep,tI,Ipulses,psi_ini,lowerBC,IC):
    
    global pars
    
    # Time grid:
    tN=float(RunTime)

    # Spatial grid:
    dz=float(SpaceStep)
    zN=float(SoilDepth)
    z=np.arange(dz/2,zN,dz)
    n=len(z)
    # z=np.hstack([0,z,zN])
    #z=z[-1]-z

    # Initial condition:
    if IC=='HS':
        psi0=z-zN
    elif IC=='FP':
        psi0=np.zeros(n)+float(psi_ini)
    else:
        print('error')

    
    
    psiB=float(psi_ini)
    
    dt=float(TimeStep)
    t=np.arange(0,tN+dt,dt)
    nt=len(t)

    # Boundary conditions:
    tI=[float(i) for i in tI.split(',')]
    Ipulses=[float(i) for i in Ipulses.split(',')]

    I=np.zeros(len(t))
    c=0

    for ti in range(len(t)):
        if t[ti]>=tI[c+1]:
            c+=1
        I[ti]=Ipulses[c]
    
    
    # I=(0.5-np.cos(t*2*np.pi/365)/2)*(maxInf-minInf)+minInf
    
    BC_T=I+np.zeros(nt)
    BC_B=psiB+np.zeros(nt)
    
    if lowerBC=='FD':
        psi,WB,runtime=run_RE_FD(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)
    elif lowerBC=='FP':
        psi,WB,runtime=run_RE_FP(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)
    elif lowerBC=='ZF':    
        psi,WB,runtime=run_RE_ZF(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)

    theta=thetaFun(psi,pars)    
    return t,z,psi,WB,theta
