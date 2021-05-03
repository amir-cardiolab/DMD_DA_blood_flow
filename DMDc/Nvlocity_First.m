%Data-Driven Pulsatile Blood Flow Physics with Dynamic Mode Decomposition
%by Milad Habibi, Scott Dawson, Amirhossein Arzani
%Paper: https://www.mdpi.com/2311-5521/5/3/111
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code used  coronary artery stenosis data of section 2.1.2: 
%Velocity data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the code generates the raw data for DMD modes and could be converted to VTK or other formats for visualization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code Courtesy of Mr. Milad Habibi (from Arzani lab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need the optimal_SVHT_coef.m code to run

clc
clear all;
close all;

data1=load ('V');
Velocity=data1.velocity;
format long
Ups=load('U');
Ups=Ups.Ups;
Velocity1=Velocity(:,1:end);
X   = (Velocity1(:,1:end-1));
Xp  = (Velocity1(:,2:end));
Ups1=Ups(1:(end-1));
Omega = [X;Ups1];
dt=0.95/1500;
[U,Sig,V] = svd(Omega,'econ');
SIG1=diag(Sig);
beta=(size(Sig,2)/size(Sig,1));
thresh1=optimal_SVHT_coef(beta, 0)*median(SIG1);
rtil =  length(find(diag(Sig)>thresh1));% p in the algorithm
Util    = U(:,1:rtil); 
Sigtil  = Sig(1:rtil,1:rtil);
Vtil    = V(:,1:rtil); 
[U,Sig,V] = svd(Xp,'econ');
SIG2=diag(Sig);
beta=(size(Sig,2)/size(Sig,1));
thresh2=optimal_SVHT_coef(beta, 0)*median(SIG2);
r=rtil;
Uhat    = U(:,1:r); 
Sighat  = Sig(1:r,1:r);
Vbar    = V(:,1:r); 
n = size(X,1); 
n2=size(X,2);
q = size(Ups1,1);
U_1 = Util(1:n,:);
U_2 = Util(n+1:end,:);
inv_Sigtil = inv(Sigtil);
approxA = Uhat'*(Xp)*Vtil*inv_Sigtil*U_1'*Uhat;
approxB =(Xp)*Vtil*inv_Sigtil*U_2';
[W,D] = eig(approxA);
landa=diag(D);
omega2=log(landa)/dt;
array_temp = U_1'*Uhat;
Phi = Xp * Vtil * inv_Sigtil * array_temp * W;%DMD modes
rank(Phi);
[nn,mm] = size(Phi);
mm1=size (X,2);
inv_Phi=pinv(Phi);

XDMD4=zeros(size(X));
XDMD5=zeros(size(X));
XDMD4(:,1)=X(:,1);
XDMD5(:,1)=X(:,1);
for iter=2:mm1
     XDMD4(:,iter)=((Phi*D)*((inv_Phi)*XDMD4(:,iter-1)))+approxB*Ups1(iter-1); 
     XDMD5(:,iter)=((Xp)*Vtil*inv_Sigtil)*(U_1'*XDMD5(:,iter-1))+approxB*Ups1(iter-1);
end
% ZZ=((abs(X-real(XDMD4))));
% ZZ2=((abs(X-real(XDMD5))));
% ZZZ=sum(sum(ZZ))/(n*n2);
% ZZZ2=sum(sum(ZZ2))/(n*n2);



