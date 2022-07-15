import numpy as np
from numba import jit 

@jit(nopython=True)
def Cfun(psi,pars):
    x3=(pars['alpha']+np.abs(psi)**pars['beta'])**2.
    x2=pars['alpha']/x3
    x1=pars['beta']*np.abs(psi)**(pars['beta']-1)*x2
    C=(pars['thetaS']-pars['thetaR'])*x1
    return C

@jit(nopython=True)
def Kfun(psi,pars):
    x2=pars['A']+np.abs(psi)**pars['gamma']
    x1=pars['A']/x2
    K=pars['Ks']*x1
    return K

@jit(nopython=True)
def thetafun(psi,pars):
    x3=pars['alpha']+np.abs(psi)**pars['beta']
    x2=pars['alpha']/x3
    x1=(pars['thetaS']-pars['thetaR'])*x2
    theta=pars['thetaR']+x1
    return theta

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
    A=np.diag(a[1:],-1)+np.diag(b,0)+np.diag(c[:-1],1)

    # Construct RHS
    y[:]=R[:]

    # Solve:
    dell = np.linalg.solve(A, y)

    return dell

@jit(nopython=True)
def Rfun(psiiter,psiin,psiT,psiB,C,Kmid,dtheta,dt,dz,n):
    # This solves the Picard residual term:
    psigrid=np.hstack((psiB,psiiter,psiT))
    
    x1=dtheta/dt*dz
    x2=-(Kmid[1:]-Kmid[:-1])
    x3=-Kmid[1:]*(psigrid[2:]-psigrid[1:-1])/dz
    x4=Kmid[:-1]*(psigrid[1:-1]-psigrid[:-2])/dz

    R=x1+x2+x3+x4

    return R

@jit(nopython=True)
def iterfun(psiin,pars,psiT,psiB,dt,dz,n):
    # psiin = psi^n
    # psiiter = psi^n+1,m
    # psiout = psi^n+1,m+1

    tolerance=1e-10
    maxcount=1000
    Rmax=1.

    # Initialize arrays
    psiiter=np.zeros(len(psiin))
    psiout=np.zeros(len(psiin))

    # Initial guess: psi_n+1^1 = psi_n
    psiiter[:]=psiin[:]

    count=0.
    while count <= 1 or (Rmax >= tolerance and count<= maxcount):
        # Get C,K:
        C=Cfun(psiiter,pars)
        K=Kfun(np.hstack((psiB, psiiter, psiT)),pars)
        Kmid=(K[1:]+K[:-1])/2.
        dtheta=thetafun(psiiter,pars)-thetafun(psiin,pars)
        # Get R
        R=Rfun(psiiter,psiin,psiT,psiB,C,Kmid,dtheta,dt,dz,n)
        # Solve for del
        dell=solverfun(R,C,Kmid,dt,dz,n)
        # Update psi estimates at different iteration levels
        psiout[:]=psiiter[:]+dell[:]
        
        psiiter[:]=psiout[:]
        Rmax=np.abs(np.max(R))
        count+=1

    #print('Iteration count = %d'%(count-1))

    return psiout

@jit(nopython=True)
def massbal(psi,psiT,psiB,pars,n,dt,dz):

    # Initial storage:
    theta=thetafun(psi.reshape(-1),pars)
    theta=np.reshape(theta,psi.shape)
    S=np.sum(theta*dz,1)
    S0=S[0]
    SN=S[-1]

    # Inflow:
    Kin=(Kfun(psiB,pars)+Kfun(psi[:,0],pars))/2.
    QIN=-Kin*((psi[:,0]-psiB)/dz+1.)
    QIN[0]=0.
    QINsum=np.sum(QIN)*dt

    # Outflow:
    Kout=(Kfun(psi[:,-1],pars)+Kfun(psiT,pars))/2.
    QOUT=-Kout*((psiT-psi[:,-1])/dz+1.)
    QOUT[0]=0.
    QOUTsum=np.sum(QOUT)*dt

    # Balance:
    dS=SN-S0
    dQ=QINsum-QOUTsum
    err=dS/dQ
    
    return QIN,QOUT,S,err

@jit(nopython=True)
def ModelRun(dt,dz,n,nt,psi,psiB,psiT,pars):
    # Solve:
    for j in range(1,nt):
        psi[j,:]=iterfun(psi[j-1,:],pars,psiT,psiB,dt,dz,n)

    QIN,QOUT,S,err=massbal(psi,psiT,psiB,pars,n,dt,dz)

    return psi,QIN,QOUT,S,err
