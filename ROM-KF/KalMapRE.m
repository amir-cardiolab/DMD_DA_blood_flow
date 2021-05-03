function [Xbar1,Xbar2]=KalMapRE(X,Xhat,A,Q,R,H,U)
% Q=White noise of process, R= White noise of measurment, Xhat=measurment ,
% X=system
tic
H=H*U;
Q=U'*Q*U;
%X=U'*X;
MPp=cell(size(X,2),1);
MPc=cell(size(X,2),1);
Qz=Q;
Rz=R;
Az=A;
Zp(:,1)=X(:,1);
Pp=Qz;
K=(Pp*H')/(H*Pp*H'+Rz);
Zc(:,1)=Zp(:,1)+K*(Xhat(:,1)-H*Zp(:,1));
Pc=(eye(size(K*H))-K*H)*Pp;
MPp{1}=Pp;
MPc{1}=Pc;
for i=2:size(X,2)
    Zp(:,i)=Az*Zc(:,i-1);
    Pp=Az*Pc*(Az')+Qz; 
    K=(Pp*H')/(H*Pp*H'+Rz);
    Zc(:,i)=Zp(:,i)+K*(Xhat(:,i)-H*Zp(:,i));
    Pc=(eye(size(K*H))-K*H)*Pp;
    MPp{i}=Pp;
    MPc{i}=Pc;
end
%Kalman Smoother
Zs=zeros(size(Zc));
Zs(:,size(X,2))=Zc(:,size(X,2));
%Zs(:,1)=Zp(:,1);
Ps=MPc{size(X,2)};
for i=size(X,2)-1:-1:1
    Ks=MPc{i}*(A')*((inv(MPp{i+1})));
    Ps=MPc{i}-Ks*(MPp{i+1}-Ps)*(Ks');
    Zs(:,i)=Zc(:,i)-Ks*(Zp(:,i+1)-Zs(:,i+1));
end
Xbar2=Zs;
Xbar1=Zc;
toc