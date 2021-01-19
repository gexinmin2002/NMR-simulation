function M=offcpmg(n,T2,T1,tau,dw,w1,flip1,flip2,phi1,phi2,ispulse)
if nargin<11
    ispulse=0;
end
if nargin<10
    phi2=0/180*pi;
end
if nargin<9
    phi1=0/180*pi;
end
if nargin<8
    flip2=pi;
end
if nargin<7
    flip1=pi/2;
end
if nargin<6
    w1=1;
end
if nargin<5
    dw=0;
end
if nargin<4
    tau=0.05;
end
if nargin<3
    T1=150;
end
if nargin<2
    T2=100;
end
if nargin<1
    n=1000;
end
%%simulate the infleunce of off resonances

M0=[0 0 1]';
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef 
tp1=flip1/w1;
tp2=flip2/w1;


[A1,B1] = freeprecess(tau+ispulse*tp1,T1,T2);

[A2,B2] = freeprecess(tau+ispulse*tp2,T1,T2);

R1=rotxn(flip1,phi1,w1,dw);
M0=R1*M0;
gama1=dw*(tau+ispulse*(tp1));
RTZ1=[cos(gama1) sin(gama1) 0; -sin(gama1) cos(gama1) 0; 0 0 1];

gama2=dw*(tau+ispulse*(tp2));
RTZ2=[cos(gama2) sin(gama2) 0; -sin(gama2) cos(gama2) 0; 0 0 1];
gama=dw*tau;
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];
M0=RTZ*A1*M0+B1;

R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
M1=RTZ2*A2*M0+B2;

[A,B] = freeprecess(tau,T1,T2);
%define the number of echos
n=1000;
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    
    M(:,i)=RTZ*A*M(:,i-1)+B;
    
    M(:,i)=R2*M(:,i);
    
    M(:,i)=RTZ2*A2*M(:,i)+B2;
end
