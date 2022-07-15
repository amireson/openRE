function HorizontalInfiltrationModel
%Set color to previous release color order
co = [0 0 1;    
    0 0.5 0;
    1 0 0;
    0 0.75 0.75;
    0.75 0 0.75;
    0.75 0.75 0;
    0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co)

SeI=0.01; %Initial effective saturation
Se0=0.99; %Boundary effective saturation

%Van genuchten paramaters
SoilName{1}='a) Hygiene sandstone';
SoilName{2}='b) Touchet Silt Loam G.E. 3';
SoilName{3}='c) Silt Loam G.E. 3';
SoilName{4}='d) Guelph Loam (drying)';
SoilName{5}='e) Beit Netofa Clay';
fname{1}='sand';
fname{3}='siltloam';
fname{5}='clay';
thS=[0.250 0.469 0.396 0.520 0.446];
thR=[0.153 0.190 0.131 0.218 0];
Ks=[108 303 4.96 31.6 0.082]; %(cm/day)
alpha=[0.0079 0.005 0.00423 0.0115 0.00152]; %(1/cm)
n=[10.4 7.09 2.06 2.03 1.17];
eta=0.5;
%Dimensionlise results
m=1-1./n;
for k=1:numel(m)
    [phi(:,k),Se,sig2(k)]=ModelFun(SeI,Se0,m(k),eta);
end
z=phi.*sqrt(Ks./alpha./(thS-thR));
theta=(thS-thR).*Se+thR;
thetaI=(thS-thR).*SeI+thR;
S=sqrt((thS-thR).*Ks./alpha.*sig2);
%Define times of interest
% t=[10 20 30 40 50]/60/24;
t=[5 10 20 50 100]/60/24;
xMAX=[100 90 12 18 1.0];

%figure(1)
%clf
for k=[1,3,5]%Choose which soil you want to look at
    %Determine horizontal distance (cm)
    xPLOT=[z(:,k).*sqrt(t); repmat(xMAX(k),1,5)];
    thetaPLOT=[theta(:,k); thetaI(k)];
    %subplot(3,2,k)
    mat=[thetaPLOT xPLOT];
    csvwrite([fname{k},'.csv'],mat);
    %plot(xPLOT,thetaPLOT,'.-','linewidth',1,'markersize',5)
    %hold on
    [zMOL,thetaMOL]=MOLFun(t,xMAX(k),thR(k),thS(k),alpha(k),n(k),Ks(k),eta,SeI,Se0);
    %plot(zMOL,thetaMOL,'g-o')
    %xlabel('Distance (cm)')
    %ylabel('Moisture content (-)')
    %legend('10 min','20 min','30 min','40 min','50 min','MOL','location','southwest')
    %title([SoilName{k} ' (S = ' num2str(S(k),4) ' cm day^1^/^2)'])
    %set(gca,'xscale','log')
    %xlim([1 xMAX(k)])
    
    %xlswrite('PsuedospectralData.xlsx',[thetaPLOT xPLOT],SoilName{k})
end

%**************************************************************************

function [z,theta]=MOLFun(t,xMAX,thR,thS,alpha,n,Ks,eta,SeI,Se0)

%Define grid points
N=200;
zB=linspace(0,xMAX,N+1)';
z=(zB(2:end,1)+zB(1:end-1,1))/2;

%Define boundary and initial conditions
m=1-1/n;
psi0=-(Se0.^(-1/m)-1).^(1/n)/alpha;
psiI=-(SeI.^(-1/m)-1).^(1/n)/alpha;
psiI=repmat(psiI,size(z,1),1);

%Solve Richards Equation
JPat=spdiags(ones(N,3),[-1 0 1],N,N);
options=odeset('JPattern',JPat,'MaxStep',1);
%wait=waitbar(0,'Please wait...');
wait=0;
[t,psi]=ode15s(@MyODEFUN,[0 t],psiI,options,z,zB,thR,thS,alpha,n,Ks,eta,psi0,t(1),t(end),wait);
%close(wait)

%De     ter     mine moisture contents
Se=(1+abs(alpha*psi).^n).^-m;
theta=Se.*(thS-thR)+thR;
%Delete initial conditions
theta(1,:)=[];

%**************************************************************************

function [dpsidt,STORAGE,theta,q]=MyODEFUN(t,psi,z,zB,thR,thS,alpha,n,Ks,eta,psi0,t0,tMax,wait)
%Add psi0 into psi vector
psi=[psi0;psi];
z=[zB(1);z];
%Apply soil moisture characteristics
Ss=0; %The value of this makes a difference with the clay so I set it to zero
m=1-1/n;
Se=(1+abs(alpha*psi).^n).^-m;
dSedpsi=m./(1-m).*Se.*(Se.^(1/m)-1)./psi;
ind=psi>=0;
Se(ind)=1;
dSedpsi(ind)=0;
K=Ks*Se.^eta.*(1-(1-Se.^(1/m)).^m).^2;
theta=Se.*(thS-thR)+thR;
C=theta/thS*Ss+(thS-thR).*dSedpsi;
%Calculate hydraulic head noting that z is depth
h=psi-z*0; %Assume no gravity for horizontal infiltration
%Calculate head gradient
dhdz=diff(h,1,1)./diff(z,1,1);
%Calculate Darcy flux
K0=K(1);
KB=(K(2:end,:)+K(1:end-1,:))/2;
KB(1)=K0;
q=-KB.*dhdz;
%Apply boundary conditions
q=[q;0];
%Remove boundary node
psi(1)=[];
C(1)=[];
theta(1)=[];
%Apply continuity equation
dthetadt=-diff(q,1,1)./diff(zB,1,1);
%Apply chain rule to get change in psi
dpsidt=dthetadt./C;
%Calculate total storage
STORAGE=sum(theta.*diff(zB,1,1),1);
%Monitor progress
%if ~isempty(wait)
%    waitbar((t-t0)/(tMax-t0),wait)
%end

%**************************************************************************

function [phi,th,sig2]=ModelFun(thI,th0,m,eta)
N=100; %Number of Chebyshev nodes
[z,D]=chebdif(N,2); %Get differentitation matrices
dzdth=2/(th0-thI); %Chebyshev node scaling factor
E1=dzdth*D(:,:,1); %First-order
E2=dzdth^2*D(:,:,2); %Second-order
%Determine coefficients for integration
IntCoefs=pi/(N-1)/dzdth*sqrt(1-z.^2)';
I=eye(N); %Identity matrix
%Determine theta values for each z value
th=(th0+thI)/2+(th0-thI)/2*z;
%Determine diffusivity for each z value
L=(1-th.^(1/m)).^m;
Dbar=(1-m)/m*th.^(eta-1/m).*(1-L).^2./L;
OF=1; %Initialise objective function
i=2:N-1; %Inner node index
F=ones(N,1); %Initial guess
while OF>1e-6 %Newton iteration
    %Determine square of sorptivity
    sig2=IntCoefs*[2*(th-thI).*Dbar./F];
    Q=2*Dbar/sig2./F;
    R=[E2(i,:)*F+I(i,:)*Q;F(N)-0;F(1)-1];
    dR=[E2(i,:)+I(i,:)*diag(-Q./F);I(N,:);I(1,:)];
    Fold=F; %Store previous iteration
    F=max(eps,F-dR\R); %Update F and ensure > 0
    OF=max(abs(F-Fold)); %Define objective function
end
%Determine phi for each theta value
phi=sqrt(sig2)*E1*F;

%**************************************************************************

function [x, DM] = chebdif(N, M)

%  The function [x, DM] =  chebdif(N,M) computes the differentiation
%  matrices D1, D2, ..., DM on Chebyshev nodes.
%
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced
%  accuracy suggested by W. Don and S. Solomonoff in
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric
%  identities to avoid the computation of differences
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:   "Spectral Differencing with a Twist", by
%  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp.

%  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by
%  JACW, May 2003.

I = eye(N);                          % Identity matrix.
L = logical(I);                      % Logical identity matrix.

n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

k = [0:N-1]';                        % Compute theta vector.
th = k*pi/(N-1);

x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

T = repmat(th/2,1,N);
DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity.
DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick.
DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.

C = toeplitz((-1).^k);               % C is the matrix with
C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))
Z(L) = zeros(N,1);                      % with zeros on the diagonal.

D = eye(N);                          % D contains diff. matrices.

for ell = 1:M
    D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
    D(L) = -sum(D');                            % Correct main diagonal of D
    DM(:,:,ell) = D;                                   % Store current D in DM
end

%**************************************************************************

function p = chebint(fk, x)

%  The function p = chebint(fk, x) computes the polynomial interpolant
%  of the data (xk, fk), where xk are the Chebyshev nodes.
%  Two or more data points are assumed.
%
%  Input:
%  fk:  Vector of y-coordinates of data, at Chebyshev points
%       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  x:   Vector of x-values where polynomial interpolant is to be evaluated.
%
%  Output:
%  p:    Vector of interpolated values.
%
%  The code implements the barycentric formula; see page 252 in
%  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
%  (Note that if some fk > 1/eps, with eps the machine epsilon,
%  the value of eps in the code may have to be reduced.)

%  J.A.C. Weideman, S.C. Reddy 1998

fk = fk(:); x = x(:);                    % Make sure data are column vectors.

N = length(fk);
M = length(x);

xk = sin(pi*[N-1:-2:1-N]'/(2*(N-1)));    % Compute Chebyshev points.

w = ones(N,1).*(-1).^[0:N-1]';          % w = weights for Chebyshev formula
w(1) = w(1)/2; w(N) = w(N)/2;

D = x(:,ones(1,N)) - xk(:,ones(1,M))';  % Compute quantities x-x(k)
D = 1./(D+eps*(D==0));                  % and their reciprocals.

p = D*(w.*fk)./(D*w);                   % Evaluate interpolant as
% matrix-vector products.
