%%simulate the infleunce of off resonance
%%%gexinmin2002 in usa 2019/06/11 程序正确无误！！
%{给出初始磁化矢量M_0并令其归一化为1；,时间间隔为tau=0.1的驰豫矩阵T，时间间隔为2*tau=0.2的驰豫矩阵T_2，沿X轴扳转90°后得到的磁化矢量M_1,A是绕Y轴旋转180°的旋转矩阵%}
clc;
clear;
tau = 0.05;		% 半回波间隔.
T1=150;%ms
T2=100;%ms
% 引用函数计算
[A,B] = freeprecess(tau,T1,T2);
%初始平衡状态
M0=[0 0 1]';
w1=1;%射频脉冲的频率！！
dw=0.5*w1;%脉冲偏移频率（z方向的dB导致，w=2*pi*f,可以算出来，dw和w1都是相对值)
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef
%x方向施加90°射频脉冲 theta=90
%x的相位角为0，phi=0;
phi1=0;
flip=pi/2;
t90=flip/w1;
%引用函数
R=rotxn(flip,phi1,w1,dw);
M0=R*M0;
%z轴旋转矩阵
gama=dw*tau;
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];
% 驰豫，tau时刻
M0=RTZ*A*M0+B;
%施加y方向的180°射频脉冲 相位为y
flip2=pi;
t180=flip/w1;
phi2=0;%y的相位角为0
R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
%%%经过tau开始驰豫，记录第一个回波幅度
M1=RTZ*A*M0+B;
%%开始循环
%define the number of echos
n=1000;%回波个数
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    %%驰豫tau时刻
    M(:,i)=RTZ*A*M(:,i-1)+B;
    %%加上y方向180°脉冲，绕有效轴旋转
    M(:,i)=R2*M(:,i);
    %%%驰豫tau时刻，记录下值
    M(:,i)=RTZ*A*M(:,i)+B;
end
M=M';