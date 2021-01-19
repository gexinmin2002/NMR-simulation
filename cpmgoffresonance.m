%%simulate the infleunce of off resonance
% Give the intial normalized magnetization and the relaxation matrix T for tau=0.1ms, and the relaxation matrix T_2 for tau=0.2ms
%M_1 is the magnetization after the first 90-x RF pulse
%A is the rotation frame for the 180-y RF pulse

clc;
clear;
tau = 0.05;		% half echo spacing.
T1=150;%ms
T2=100;%ms
% Define the relaxation matrix
[A,B] = freeprecess(tau,T1,T2);
% Define the intial magnetization
M0=[0 0 1]';

w1=1;% frequency of the RF pulse 
dw=0.5*w1;% off resonance frequency（ contributed by the dB in z-axis，w=2*pi*f，dw and w1 are relative values)
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef
% introduce the 90-x RF pulse
%x的相位角为0，phi=0;
phi1=0;% phase angle of the excitation pulse
flip=pi/2;
t90=flip/w1;% time duration of the 90-x RF pulse

R=rotxn(flip,phi1,w1,dw);
M0=R*M0;
gama=dw*tau;
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];
% Relaxation after tau
M0=RTZ*A*M0+B;
%Give the 180-y RF pulse
flip2=pi;% refocusing angle
t180=flip/w1;%duration of the refocusing pulse
phi2=0;% phase angle for the refocusing angle
R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
%Measuring the first echo after the tau relaxation
M1=RTZ*A*M0+B;
%%Start the iteration
%define the number of echos
n=1000;%Number of echoes
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    M(:,i)=RTZ*A*M(:,i-1)+B;
    M(:,i)=R2*M(:,i);
        M(:,i)=RTZ*A*M(:,i)+B;
end
M=M';