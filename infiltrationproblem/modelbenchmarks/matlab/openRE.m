function openRE()

% Run the model

% Load some driving data (flux across upper boundary):
qT=csvread('../../input/infiltration.dat',1,1);
qT=qT(:,1)/1000;
qB=zeros(numel(qT),1);

t=1:numel(qT);
t=t(:);
dt=t(2)-t(1);

size(t)
size(qT)

% Soil properties:;
pars.thetaR=0.131;
pars.thetaS=0.396;
pars.alpha=0.423;
pars.n=2.06;
pars.m=1-1/pars.n;
pars.Ks=0.0496;
pars.neta=0.5;
pars.Ss=0.000001;


% Spatial grid:;
dz=0.1;
zN=1.5;
z=dz/2:dz:zN;
n=numel(z);

% Initial condition:;
psi0=zeros(n,1)-3.59;

% Solve
[psi,WB,runtime]=run_RE(dt,t,dz,zN,n,psi0,qT,qB,pars);

figure()
hold on
plot(WB.t,1000*(WB.S-WB.S(1)),'-b')
plot(WB.t,1000*cumsum(WB.QIN-WB.QOUT),'.r')
legend('Cumulative change in storage','Cumulative net flux')
grid on
xlabel('Time (days)')
ylabel('Storage (mm)')

end

function [psi,WB,runtime]=run_RE(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars);

    DV=zeros(numel(t),n+2);
    DV(1,1)=0.;          % Cumulative inflow
    DV(1,end)=0.;        % Cumulative outflow
    DV(1,2:end-1)=psi0;  % Matric potential
    
    JPat=spdiags(ones(n+2,3),[-1 0 1],n+2,n+2);
    options=odeset('JPattern',JPat);
    tic;
    for i=1:numel(t)-1;
        [~,DVdum]=ode15s(@odefun,[0,dt],DV(i,:),options,pars,n,BC_T(i),BC_B(i),dz,DV(i,:));
        DV(i+1,:)=DVdum(end,:);
    end
    
    runtime=toc
    

    % Unpack output:
    QT=DV(:,1);
    QB=DV(:,end);
    psi=DV(:,2:end-1);
    qT=[0;diff(QT)];
    qB=[0;diff(QB)];

    % Water balance terms
    theta=thetaFun(psi,pars);
    S=sum(theta*dz,2);

    % Pack output into a dataframe:
    WB.t=t;
    WB.S=S;
    WB.QIN=qT;
    WB.QOUT=qB;
 
end


% vgprops.py
function [theta]=thetaFun(psi,pars)
    Se=(1+(psi*-pars.alpha).^pars.n).^(-pars.m);
    Se(psi>0.)=1.0;
    theta=pars.thetaR+(pars.thetaS-pars.thetaR).*Se;
end

function [C]=CFun(psi,pars)
    Se=(1+(psi*-pars.alpha).^pars.n).^(-pars.m);
    Se(psi>0.)=1.0;
    dSedh=pars.alpha*pars.m/(1-pars.m).*Se.^(1/pars.m).*(1-Se.^(1/pars.m)).^pars.m;
    C=Se*pars.Ss+(pars.thetaS-pars.thetaR)*dSedh;
end

function [K]=KFun(psi,pars)
    Se=(1+(psi*-pars.alpha).^pars.n).^(-pars.m);
    Se(psi>0.)=1.0;
    K=pars.Ks*Se.^pars.neta.*(1-(1-Se.^(1/pars.m)).^pars.m).^2;
end

% Cinv_AN.py
function [Cinv]=CinvFun(psi,psi_n,pars)
    Cinv=1./CFun(psi,pars);
end

% BC_t2FD.py
function [qT,qB]=BoundaryFluxes(BC_T,BC_B,pars,dz,psiTn,psiBn)
    % Inputs:
    %  BC_T = specified flux at surface or specified pressure head at surface;
    %  BC_B = specified flux at base or specified pressure head at base;
    % For free drainage BC_B must be an arbitrary array, that is not used.
    %  pars = soil hydraulic properties
    % psiTn = pressure head at node 0 (uppermost node)
    % psiBn = pressure head at node -1 (lowermost node)

    % Upper BC: Type 2 specified flux
    qT=BC_T;

    % Lower BC: Free drainage
    qB=KFun(psiBn,pars);
end

% richardsFlux.py
function [dDVdt]=odefun(t,DV,pars,n,BC_T,BC_B,dz,psi_n)

    % In this function, we use a block centered grid approch, where the finite difference
    % solution is functionined in terms of differences in fluxes. 

    % Unpack the dependent variable:
    QT=DV(1);
    QB=DV(end);
    psi=DV(2:end-1);

    %qT=interp(t,tT,qT)
    q=zeros(n+1,1);
    K=zeros(n+1,1);
    
    K=KFun(psi,pars);
    Kmid=(K(1:end-1)+K(2:end))/2.;
    
    % Boundary fluxes:
    [qT,qB]=BoundaryFluxes(BC_T,BC_B,pars,dz,psi(1),psi(end));
    q(1)=qT;
    q(end)=qB;

    % Internal nodes
    q(2:end-1)=-Kmid.*((psi(2:end)-psi(1:end-1))/dz-1);

    % Continuity
    Cinv=CinvFun(psi,psi_n,pars);
    dpsidt=-Cinv.*(q(2:end)-q(1:end-1))/dz;

%    % Change in cumulative fluxes:
%    dQTdt=qT
%    dQBdt=qB

    % Pack up dependent variable:
    dDVdt=[qT;dpsidt;qB];

end
