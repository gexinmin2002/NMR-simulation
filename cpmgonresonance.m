%Simulate the CPMG pulse sequence at on resonance condition
clc;
clear;
tau = 0.05;		
T1=150;%ms
T2=100;%ms
[A,B] = freeprecess(tau,T1,T2);
M0=[0 0 1]';
alpha=pi/2;%90¡ãpulse
beta=2*alpha;%180¡ãpulse

n=1000;
M1=zeros(3,n);
T=[exp(-tau/T2) 0 0;0 exp(-tau/T2) 0;0 0 exp(-tau/T1)];

A1=[1 0 0;0 cos(alpha) sin(alpha);0 -sin(alpha) cos(alpha)];
A2=[cos(beta) 0 -sin(beta);0 1 0;sin(beta) 0 cos(beta)];%
B=[0;0;1-exp(-tau/T1)];
MP=T*A1*M0+B;
M1(:,1)=T*A2*MP+B;
    %%%start the iteration
 for i=2:n
     M1(:,i)=T*M1(:,i-1)+B;
     M1(:,i)=T*A2*M1(:,i)+B;
 end