%Integrating multi-fidelity blood flow data with reduced-order data assimilation
%by Milad Habibi, Roshan M D'Souza, Scott Dawson, Amirhossein Arzani
%Paper: https://arxiv.org/abs/2104.01971
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code used in Section 2.1 of the paper: 
%Test case 1: Womersley?s analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code Courtesy of Mr. Milad Habibi (from Arzani lab)


clc
clear all
close all
Cn=load('Cn');
ab=Cn.Cn;
fr=load('fr');
w=fr.fr;
a=ab(:,1);
b=ab(:,2);
a0=a(1);
b0=b(1);
a1=a(2);
b1=b(2);
a2=a(3);
b2=b(3);
a3=a(4);
b3=b(4);
a4=a(5);
b4=b(5);
a5=a(6);
b5=b(6);
a6=a(7);
b6=b(7);
a7=a(8);
b7=b(8);
a8=a(9);
b8=b(9);
mu=0.0035;
ru=1060;
nu=3.302e-6;
R=0.012;
r=linspace(0,R,101);
VV=zeros(101,101);
alpha=R*sqrt(w/nu);
t=linspace(0,1,101);
for k=1:101
    Vs=a0*(r(k)^2-(R)^2)/(4*mu);
    lamda1=sqrt(w*(1i^3)/nu);
    lamda2=sqrt(2*w*(1i^3)/nu);
    lamda3=sqrt(3*w*(1i^3)/nu);
    lamda4=sqrt(4*w*(1i^3)/nu);
    lamda5=sqrt(5*w*(1i^3)/nu);
    lamda6=sqrt(6*w*(1i^3)/nu);
    lamda7=sqrt(7*w*(1i^3)/nu);
    lamda8=sqrt(8*w*(1i^3)/nu);
    for kk=1:101
    Vu1=((a1-b1*1i)/(1i*w*ru))*((besselj(0,r(k)*lamda1)/(besselj(0,R*lamda1))-1))*exp(1i*w*t(kk));
    Vu2=((a2-b2*1i)/(1i*2*w*ru))*((besselj(0,r(k)*lamda2)/(besselj(0,R*lamda2))-1))*exp(1i*2*w*t(kk));
    Vu3=((a3-b3*1i)/(1i*3*w*ru))*((besselj(0,r(k)*lamda3)/(besselj(0,R*lamda3))-1))*exp(1i*3*w*t(kk));
    Vu4=((a4-b4*1i)/(1i*4*w*ru))*((besselj(0,r(k)*lamda4)/(besselj(0,R*lamda4))-1))*exp(1i*4*w*t(kk));
    Vu5=((a5-b5*1i)/(1i*5*w*ru))*((besselj(0,r(k)*lamda5)/(besselj(0,R*lamda5))-1))*exp(1i*5*w*t(kk));
    Vu6=((a6-b6*1i)/(1i*2*w*ru))*((besselj(0,r(k)*lamda6)/(besselj(0,R*lamda6))-1))*exp(1i*6*w*t(kk));
    Vu7=((a7-b7*1i)/(1i*7*w*ru))*((besselj(0,r(k)*lamda7)/(besselj(0,R*lamda7))-1))*exp(1i*7*w*t(kk));
    Vu8=((a8-b8*1i)/(1i*8*w*ru))*((besselj(0,r(k)*lamda8)/(besselj(0,R*lamda8))-1))*exp(1i*8*w*t(kk));
    Vus(k,kk)=real(Vu1)+real(Vu2)+real(Vu3)+real(Vu4)+real(Vu5)+real(Vu6)+real(Vu7)+real(Vu8);
    VV(k,kk)=(Vs+Vus(k,kk));
     end
end