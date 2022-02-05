%Integrating multi-fidelity blood flow data with reduced-order data assimilation
%by Milad Habibi, Roshan M D'Souza, Scott Dawson, Amirhossein Arzani
%Paper: https://arxiv.org/abs/2104.01971
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code used in Section 3.2 of the paper: 
% Test case 2: Idealized 2D cerebral aneurysm model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code Courtesy of Mr. Milad Habibi (from Arzani lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need the optimal_SVHT_coef.m and KalMapRE.m code to run





clc
clear all;
close all;

data1=load ('Gr');
Velocity=data1.velocity;
VelocityP=Velocity;
VelocityP(3:3:end,:)=[];
VelocityPx=VelocityP(1:2:end,:);       % CReating VX
VelocityPy=VelocityP(2:2:end,:);
VelocityP=[VelocityPx;VelocityPy];


Velocity1=load('Un');
Velocity1=Velocity1.velocity;
Velocity1(3:3:end,:)=[];
Velocity1x=Velocity1(1:2:end,:);
Velocity1y=Velocity1(2:2:end,:);
Velocity1=[Velocity1x;Velocity1y];

Velocity2=load('Ex');
Velocity2=Velocity2.velocity;
Velocity2(3:3:end,:)=[];
Velocity2x=Velocity2(1:2:end,:);
Velocity2y=Velocity2(2:2:end,:);
Velocity2=[Velocity2x;Velocity2y];
%Mapping Matrix
H=load('H');
H=H.H_matrix;
Hz=zeros(size(H));
H=[H,Hz;Hz,H];


% fbDMD
X1   = Velocity1(:,1:end-1);
X2  = Velocity1(:,2:end);
[U, S, V] = svd(X1, 'econ');

SIG1=diag(S);
dt=0.95/1500;
beta=(size(S,2)/size(S,1));
thresh1=optimal_SVHT_coef(beta, 0)*median(SIG1);
r =  length(find(diag(S)>thresh1));
%r=98;
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde1 = U_r' * X2 * V_r / S_r;
%backward DMD
[U2, S2, V2] = svd(X2, 'econ');

SIG2=diag(S2);
beta2=(size(S2,2)/size(S2,1));
thresh2=optimal_SVHT_coef(beta2, 0)*median(SIG2);
% r2 =  length(find(diag(S2)>thresh2))+1;
r2=r;
U_r2 = U2(:, 1:r2); % truncate to rank-r
S_r2 = S2(1:r2, 1:r2);
V_r2 = V2(:, 1:r2);
Atilde2 = U_r2' * X1 * V_r2 / S_r2;

Atilde=(Atilde1*inv(Atilde2))^0.5;

[W_r, D] = eig(Atilde);
Phi = X2 * V_r / S_r * W_r; % DMD modes
INVPhi=pinv(Phi);
% aa=vecnorm(Phi);
% Phi=aa.*GramSchmidt(Phi);
lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
x1 = X1(:, 1);
b = Phi\x1;
%% DMD reconstruction
mm1 = size(X1, 2)+1; % mm1 = m - 1
time_dynamics = zeros(r, mm1);
t = (0:mm1-1)*dt; % time vector
for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = Phi * time_dynamics;

XDMD3=zeros(r,size(X1,2));
XDMD3(:,1)=U_r'*X1(:,1);
Atilde=U_r'*Phi*D*INVPhi*U_r;
for i=2:mm1
    XDMD3(:,i)=Atilde*XDMD3(:,i-1);
end
ZZ2=((abs(Velocity1-real(Xdmd))));
ZZZ2=sum(sum(ZZ2))/(size(Velocity1,1)*size(Velocity1,2));
SI1=load('Co');
SI1=SI1.MEANCOV;
cov1=diag(SI1);

SI2=(1*5.5064)^2;
cov2=SI2*eye(size(Velocity2,1));
A=(X2* V_r / S_r)*(U_r');

[Velocity4,Velocity5]=KalMapRE(XDMD3,Velocity2,Atilde,cov1,cov2,H,U_r);
Velocity4=U_r*Velocity4;
Velocity5=U_r*Velocity5;


%KALMAN FILTER
ERK=(real(Velocity4)-VelocityP).^2;
ER2K=((sum(ERK))/(size(ERK,1))).^0.5;
ER3K=ER2K./max(abs(VelocityP));


%KALMAN SMOOTHER
ER=(real(Velocity5)-VelocityP).^2;
ER2=((sum(ER))/(size(ER,1))).^0.5;
ER3=ER2./max(abs(VelocityP));


% CFD
ERCFD=(real(Velocity1)-VelocityP).^2;
ER2CFD=((sum(ERCFD))/(size(ERCFD,1))).^0.5;
ER3CFD=ER2CFD./max(abs(VelocityP));


% MEASURMENT
ER4d=(real(Velocity2)-H*VelocityP).^2;
ER24d=((sum(ER4d))/(size(ER4d,1))).^0.5;
ER34d=ER24d./max(abs(VelocityP));

