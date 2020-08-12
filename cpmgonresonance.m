%gexinnmin 2019/06/011 houston 程序正确无误！
%%%数值模拟CPMG脉冲序列
%M1是在不同扳倒角误差下的磁化矢量矩阵，是三维矩阵，size(3,m,n);
%T1 T2是核磁共振的纵横向驰豫时间；
%r1 r2是最小和最大扳倒角误差；
%n是r1至r2之间的线性布点数
%m是采样点数，也就是回波个数
%tau是半回波间隔
clc;
clear;
tau = 0.05;		% 半回波间隔
T1=150;%ms
T2=100;%ms
% 引用函数计算
[A,B] = freeprecess(tau,T1,T2);
%初始平衡状态
M0=[0 0 1]';
alpha=pi/2;%90°脉冲
beta=2*alpha;%180°脉冲
%%重新赋值为90度，只研究y轴扳倒的影响！
n=1000;%表示回波个数
M1=zeros(3,n);%矩阵，用于存储驰豫数组；3表示x y z;m表示m个回波；
T=[exp(-tau/T2) 0 0;0 exp(-tau/T2) 0;0 0 exp(-tau/T1)];%驰豫矩阵
%引用函数计算
A1=[1 0 0;0 cos(alpha) sin(alpha);0 -sin(alpha) cos(alpha)];%x旋转矩阵
A2=[cos(beta) 0 -sin(beta);0 1 0;sin(beta) 0 cos(beta)];%y旋转矩阵
B=[0;0;1-exp(-tau/T1)];%z轴驰豫矩阵，on resonance，不受影响，一直存在
%经过x轴alpha角扳倒并经过tau驰豫后的信号；
MP=T*A1*M0+B;
%再经过y轴beta旋转后，再驰豫tau时刻，记录到的第一回波;
M1(:,1)=T*A2*MP+B;
    %%%开始循环
 for i=2:n
     M1(:,i)=T*M1(:,i-1)+B;%驰豫了tau
     M1(:,i)=T*A2*M1(:,i)+B;%绕y旋转了beta角，再经过tau,采集得到回波；
 end