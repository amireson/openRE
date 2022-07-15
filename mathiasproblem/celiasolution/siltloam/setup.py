pars={}
pars['thetaR']=0.131
pars['thetaS']=0.396
pars['alpha']=0.00423
pars['n']=2.06 
pars['m']=1-1/pars['n']
pars['Ks']=4.96
pars['neta']=0.5
pars['Ss']=0.

# Time grid:
dt=0.01/24./60.
tN=100./24./60.

# Spatial grid:
dz=0.05
zN=20.

# Boundary conditions:
Se0=0.01
SeI=0.99
psiT=-(SeI**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
psiB=-(Se0**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
