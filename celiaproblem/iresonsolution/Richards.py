# 01_header.py
import numpy as np
import pandas as pd
import time
from scipy.integrate import ode, solve_ivp
from numba import jit 

# 02_MyNumba.py
from numba import types
from numba.typed import Dict

def MakeDictArray():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64[:],)
    return d

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

# 03_celiaprops.py
@jit(nopython=True)
def CFun(psi,pars):
    x3=(pars['alpha']+np.abs(psi)**pars['beta'])**2.
    x2=pars['alpha']/x3
    x1=pars['beta']*np.abs(psi)**(pars['beta']-1)*x2
    C=(pars['thetaS']-pars['thetaR'])*x1
    return C

@jit(nopython=True)
def KFun(psi,pars):
    x2=pars['A']+np.abs(psi)**pars['gamma']
    x1=pars['A']/x2
    K=pars['Ks']*x1
    return K

@jit(nopython=True)
def thetaFun(psi,pars):
    x3=pars['alpha']+np.abs(psi)**pars['beta']
    x2=pars['alpha']/x3
    x1=(pars['thetaS']-pars['thetaR'])*x2
    theta=pars['thetaR']+x1
    return theta


# 04_Cinv_AN.py
@jit(nopython=True)
def CinvFun(psi,psi_n,pars):
    Cinv=1/CFun(psi,pars)
    return Cinv

# 05_BC_t1t1.py
@jit(nopython=True)
def BoundaryFluxes(BC_T,BC_B,pars,dz,psiTn,psiBn):
    # Inputs:
    #  BC_T = specified flux at surface or specified pressure head at surface;
    #  BC_B = specified flux at base or specified pressure head at base;
    #  pars = soil hydraulic properties
    # psiTn = pressure head at node 0 (uppermost node)
    # psiBn = pressure head at node -1 (lowermost node)

    # Upper BC: Type 1 specified pressure head
    psiT=BC_T
    Kin=(KFun(np.array([psiT]),pars)+KFun(np.array([psiTn]),pars))/2.
    qT=-Kin[0]*((psiTn-psiT)/dz-1.)

    # Lower BC: Type 1 specified pressure head
    psiB=BC_B
    Kout=(KFun(np.array([psiBn]),pars)+KFun(np.array([psiB]),pars))/2.
    qB=-Kout[0]*((psiB-psiBn)/dz-1.)

    return qT,qB

# 06_richardsFlux.py
# Functions called by the ODE solver:
def odefun_blockcentered(t,DV,pars,n,BC_T,BC_B,dz,psi_n):
    return odefuncall(t,DV,pars,n,BC_T,BC_B,dz,psi_n)

@jit(nopython=True)
def odefuncall(t,DV,pars,n,BC_T,BC_B,dz,psi_n):

    # In this function, we use a block centered grid approch, where the finite difference
    # solution is defined in terms of differences in fluxes. 

    # Unpack the dependent variable:
    QT=DV[0]
    QB=DV[-1]
    psi=DV[1:-1]
    psi_n=psi_n[1:-1]

    #qT=np.interp(t,tT,qT)
    q=np.zeros(n+1)
    K=np.zeros(n+1)
    
    K=KFun(psi,pars)
    Kmid=(K[1:]+K[:-1])/2.
    
    # Boundary fluxes:
    qT,qB=BoundaryFluxes(BC_T,BC_B,pars,dz,psi[0],psi[-1])
    q[0]=qT
    q[-1]=qB

    # Internal nodes
    q[1:-1]=-Kmid*((psi[1:]-psi[:-1])/dz-1)

    # Continuity
    Cinv=CinvFun(psi,psi_n,pars)
    dpsidt=-Cinv*(q[1:]-q[:-1])/dz

#    # Change in cumulative fluxes:
#    dQTdt=qT
#    dQBdt=qB

    # Pack up dependent variable:
    dDVdt=np.hstack((np.array([qT]),dpsidt,np.array([qB])))

    return dDVdt


# 08_solve_ode_RF_BDF.py
def run_RE(dt,t,dz,zN,n,psi0,BC_T,BC_B,parsIN):
    # 4. scipy function "ode", with the jacobian, solving one step at a time:
    
    pars=MakeDictFloat()
    for k in parsIN: pars[k]=parsIN[k]

    DV=np.zeros((len(t),n+2))
    DV[0,0]=0.       # Cumulative inflow
    DV[0,-1]=0.      # Cumulative outflow
    DV[0,1:-1]=psi0  # Matric potential

    r = ode(odefun_blockcentered)
    r.set_integrator('vode',method='BDF',uband=1,lband=1)
    
    tic=time.time()
    for i,ti in enumerate(t[:-1]):
        r.set_initial_value(DV[i,:], 0)

        params=(pars,n,BC_T[i],BC_B[i],dz,DV[i,:])
        #r.set_jac_params(*params)
        r.set_f_params(*params)
        r.integrate(dt)
        DV[i+1,:]=r.y

    runtime=time.time()-tic
    print('ode, with jac runtime = %.2f seconds'%(runtime))

    # Unpack output:
    QT=DV[:,0]
    QB=DV[:,-1]
    psi=DV[:,1:-1]
    qT=np.hstack([0,np.diff(QT)])/dt
    qB=np.hstack([0,np.diff(QB)])/dt

    # Water balance terms
    theta=thetaFun(psi.reshape(-1),pars)
    theta=np.reshape(theta,psi.shape)
    S=np.sum(theta*dz,1)

    # Pack output into a dataframe:
    WB=pd.DataFrame(index=t)
    WB['S']=S
    WB['QIN']=qT
    WB['QOUT']=qB

    return psi,WB,runtime

