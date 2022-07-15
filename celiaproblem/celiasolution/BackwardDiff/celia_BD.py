import numpy as np

def Cfun(psi,pars):
    x3=(pars['alpha']+np.abs(psi)**pars['beta'])**2.
    x2=pars['alpha']/x3
    x1=pars['beta']*np.abs(psi)**(pars['beta']-1)*x2
    C=(pars['thetaS']-pars['thetaR'])*x1
    return C

def Kfun(psi,pars):
    x2=pars['A']+np.abs(psi)**pars['gamma']
    x1=pars['A']/x2
    K=pars['Ks']*x1
    return K

def thetafun(psi,pars):
    x3=pars['alpha']+np.abs(psi)**pars['beta']
    x2=pars['alpha']/x3
    x1=(pars['thetaS']-pars['thetaR'])*x2
    theta=pars['thetaR']+x1
    return theta

def solverfun(psi,C,Kmid,dt,dz,n,psiT,psiB):
    # Initialize arrays
    a=np.zeros(n)
    b=np.zeros(n)
    c=np.zeros(n)
    y=np.zeros(n)

    # Construct matrix
    a=-Kmid[:-1]*dt/C/dz**2.
    b=1.+(Kmid[:-1]+Kmid[1:])*dt/C/dz**2.
    c=-Kmid[1:]*dt/C/dz**2.
    A=np.diag(a[1:],-1)+np.diag(b,0)+np.diag(c[:-1],1)

    # Construct RHS
    y[:]=psi+(Kmid[1:]-Kmid[:-1])*dt/C/dz

    # Boundary conditions
    y[0]=y[0]+psiB*Kmid[0]*dt/C[0]/dz**2.
    y[-1]=y[-1]+psiT*Kmid[-1]*dt/C[-1]/dz**2.

    # Solve:
    psinext = np.linalg.solve(A, y)

    return psinext

def iterfun(psiin,pars,psiT,psiB,dt,dz,n):
    # Here, psiin, psiiter and psiout are 1D array in z
    tolerance=1e-3
    R=1.
    psiiter=np.zeros(len(psiin))
    psiout=np.zeros(len(psiin))
    psiiter[:]=psiin[:]
    while R >= tolerance:
        C=Cfun(psiiter,pars)
        K=Kfun(np.hstack([psiB, psiiter, psiT]),pars)
        Kmid=(K[1:]+K[:-1])/2.
        psiout=solverfun(psiiter,C,Kmid,dt,dz,n,psiT,psiB)
        R=1e-10

    return psiout

def massbal(psi,psiT,psiB,pars,n,dt,dz):

    # Initial storage:
    theta=thetafun(psi,pars)
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

def ModelRun(dt,dz,n,nt,psi,psiB,psiT,pars):
    # Solve:
    for j in range(1,nt):
        psi[j,:]=iterfun(psi[j-1,:],pars,psiT,psiB,dt,dz,n)

    QIN,QOUT,S,err=massbal(psi,psiT,psiB,pars,n,dt,dz)

    return psi,QIN,QOUT,S,err
