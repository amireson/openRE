import numpy as np
from numba import jit 

@jit(nopython=True)
def thetafun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

@jit(nopython=True)
def Cfun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    theta=pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return theta/pars['thetaS']*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

@jit(nopython=True)
def Kfun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

@jit(nopython=True)
def ThomasAlg(a,b,c,n,f):

    beta=np.zeros(n)
    y=np.zeros(n)
    
    beta[0]=c[0]/b[0]
    y[0]=f[0]/b[0]
    for i in range(1,n):
        beta[i]=c[i]/(b[i]-a[i]*beta[i-1])
        y[i]=(f[i]-a[i]*y[i-1])/(b[i]-a[i]*beta[i-1])
    
    x=np.zeros(n)
    x[-1]=y[-1]
    for i in range(n-2,-1,-1):
        x[i]=y[i]-beta[i]*x[i+1]
    
    return x

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

    dell=ThomasAlg(a,b,c,n,R)
    #A=np.diag(a[1:],-1)+np.diag(b,0)+np.diag(c[:-1],1)

    # Construct RHS
    #y[:]=R[:]

    # Solve:
    #dell = np.linalg.solve(A, y)

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

    #print(count)
    #print(Rmax)
    #print('')
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
    
#    print('Delta storage = %.3f'%dS)
#    print('Flow at base (+ve upwards) = %.3f'%QINsum)
#    print('Flow at surface (+ve upwards) = %.3f'%QOUTsum)
#    print('Delta flow = %.3f'%dQ)
#    print('Error metric = %.3f'%err)
    return QIN,QOUT,S,err

@jit(nopython=True)
def ModelRun(dt,dz,n,nt,psi,psiB,psiT,pars):
    # Solve:
    for j in range(1,nt):
        #print(j)
        psi[j,:]=iterfun(psi[j-1,:],pars,psiT,psiB,dt,dz,n)

    QIN,QOUT,S,err=massbal(psi,psiT,psiB,pars,n,dt,dz)
    #QIN=0.
    #QOUT=0.
    #S=0.
    #err=0.

    return psi,QIN,QOUT,S,err
