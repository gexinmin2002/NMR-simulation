clear;
clc;
N = 128;        %inversion data point
aa = 0.1;      % left relaxation time
bb = 10000;     % right relaxation time 
x = logspace( log10( aa ), log10( bb ), N);     %define the T2 time series
%%%% synthetize the echo train %%%%
NE = 4096;
TE = 0.1;
x1 = 3;         %mean value of the first peak
x2 = 60;        %mean value of the second peak
%%%% construct the gaussian function for better peak, s should be 4 times
%%%% lower than the mean value
%%%% creat peaks
xx  = log10(x);
xx1 = log10(x1);
xx2 = log10(x2);
s1  = xx1 / 3;
s2  = xx2 / 4;
p1  = normpdf( xx, xx1, s1);
p2  = normpdf( xx, xx2, s2);
%%%% cutoff value determined  %%%% 
for i = 1 : N
    if (xx(i) < xx1 - s1 * 3) || (xx(i) > xx1 + s1 * 3)
        p1(i) = 0;
    end
end
for i = 1 : N
    if (xx(i) < xx2 - s2 * 3) || (xx(i) > xx2 + s2 * 3)
        p2(i) = 0;
    end
end
%%%% normalization %%%%
pp1 = sum(p1);
p1  = p1 / pp1;
pp2 = sum(p2);
p2  = p2 / pp2;
%%%% define the ratio of each peak %%%%
a1 =1;
a2 =0;
%%%% synthetize the echo signal %%%%
ap = a1 * p1 + a2 * p2;
%%%% define the relaxation time %%%%
for i = 1 : NE
    t(i,1) = TE * i;
end
%%%% define the echo trains with two echo spacings %%%%
Echo = zeros(NE,1);
for i = 1 : NE
    for j = 1 : N
        Echo(i) = Echo(i) + ap(j) * exp(-t(i) / x(j));
    end
end
%%%% construct the coefficient matrix A  %%%%
A = zeros( NE, N);
for i = 1 : NE 
    for j = 1 : N
        A(i,j) = exp(-t(i) / x(j));
    end 
end 

for i=1: N
    if xx(i)<0
       xx(i)=0;
     end
end
ys1=[t Echo];
%%%% noise level %%%%    
SNR =100;         % change the signal to nosie ratio here
Necho = awgn( Echo, SNR, 'measured');
ys=[t Necho];
%%%% plot 1 %%%%
xx=x';
app=ap';
m=NE;
T2=x;
T1=3*T2;
dw=0;
flip1=pi/2;
flip2=pi;
phi1=0;
phi2=0;
ispulse=0;
tau=TE/2;

ratio=100;
tp=TE/2;%ms
w1=flip1/tp;
%NE number of echos
M1=zeros(3,NE,N);

for i=1:N
    M1(:,:,i)=offcpmg(NE,T2(i),T1(i),tau,dw,w1,flip1,flip2,phi1,phi2,ispulse);
end

[a,b,c]=size(M1);
M2=zeros(a,b);
for i=1:N
    M2(:,:)=M2(:,:)+M1(:,:,i)*app(i);
end
pq=[xx app];
M22=M2';
SNR=5;
Necho1 = awgn(M22(:,2), SNR, 'measured');
