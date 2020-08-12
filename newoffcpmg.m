function M=newoffcpmg(n,T2,T1,tau,dw,w1,flip1,flip2,phi1,phi2,ispulse)
%%%gexinmin2002 in usa 2019/06/28 程序正确无误！！改进，考虑到了z轴方向旋转的影响
%%% off resonance的影响 2019/07/03新改进 function is ok
%%%注意这里面的是角频率
%%在有off resonance情况下CPMG脉冲序列数值模拟,考虑到不同的脉冲扳倒角，以及脉宽
%n number of echos 默认为4096
%T1 纵向驰豫时间
%T2 横向驰豫时间
%dw 表示z方向的 delta(B0)的off resonance 频率 detla(w0)，dw=2*pi*(df),单位可以转换成Hz
%w1 表示射频脉冲频率,w=2*pi*f，单位可以转换成Hz%  射频脉冲的频率！
%%脉冲偏移频率（z方向的dB导致，w=2*pi*f,可以算出来，dw和w1都是相对值)
%phi1 x轴的相位，表示xoy平面中从x轴算起的夹角，一般为0;
%phi2 y轴的相位，表示xoy平面中从y轴算起的夹角，一般为0;
%tau 半回波间隔 tau=TE/2 ，默认为0.1ms
%flip1 x轴的脉冲角度，如果是CPMG，一般为90°，flip1=pi;
%flip2 y轴的脉冲角度，如果是CPMG,一般为180°，flip=2*pi
%ispulse 是否考虑脉宽影响，ispulse=0 不考虑，=1考虑，默认不考虑
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
    n=4096;
end
%%simulate the infleunce of off resonances

%{给出初始磁化矢量M_0并令其归一化为1；,时间间隔为tau=0.1的驰豫矩阵T，时间间隔为2*tau=0.2的驰豫矩阵T_2，沿X轴扳转90°后得到的磁化矢量M_1,A是绕Y轴旋转180°的旋转矩阵%}
M0=[0 0 1]';
% 引用函数计算c驰豫矩阵，考虑到脉宽的影响
%初始平衡状态
%{研究DeltaB0=10%*B0时，即较弱非均匀场的影响%}
%w1=1;%射频脉冲的频率！！
%dw=0.1*w1;%脉冲偏移频率（z方向的dB导致，w=2*pi*f,可以算出来，dw和w1都是相对值)
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef 计算有效磁场的频率
%x方向施加90°射频脉冲 theta=90
%x的相位角为0，phi=0;
tp1=flip1/w1;%x轴旋转的脉宽
tp2=flip2/w1;%y轴旋转的脉宽，假定射频脉冲强度和频率不变B1 w1
%%%或者 tp1 tp2可以自己定义！
% 定义x-flip1脉冲过后经过tau+pl1的驰豫矩阵
[A1,B1] = freeprecess(tau+ispulse*tp1,T1,T2);
% 定义y-flip2脉冲过后经过tau+pl1的驰豫矩阵
[A2,B2] = freeprecess(tau+ispulse*tp2,T1,T2);
%引用函数
R1=rotxn(flip1,phi1,w1,dw);
M0=R1*M0;
%z轴旋转矩阵
%有效旋转角度
gama1=dw*(tau+ispulse*(tp1));
RTZ1=[cos(gama1) sin(gama1) 0; -sin(gama1) cos(gama1) 0; 0 0 1];
%有效旋转角度
gama2=dw*(tau+ispulse*(tp2));
RTZ2=[cos(gama2) sin(gama2) 0; -sin(gama2) cos(gama2) 0; 0 0 1];
gama=dw*(tau+ispulse*tp1);%90度脉冲+tau
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];%旋转矩阵
% 驰豫，tau时刻+pulse length
%off resonace的影响 用RTZ表示 作用时间为90°脉冲+半回波间隔
M0=RTZ*A1*M0+B1;
%施加y方向的flip2射频脉冲 相位为y
R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
%计算第一个点之前的off resoance的角度
tt=dw*(ispulse*tp2+tau);
%%计算沿着z轴的旋转
RTZZ=[cos(tt) sin(tt) 0; -sin(tt) cos(tt) 0; 0 0 1]; 
%%再次计算 off resonance角度
ttt=dw*tau;
RTZZZ=[cos(ttt) sin(ttt) 0; -sin(ttt) cos(ttt) 0; 0 0 1]; 
%%%经过tau开始驰豫，记录第一个回波幅度
M1=RTZZ*RTZ2*A2*M0+B2;
%%开始循环
[A,B] = freeprecess(tau,T1,T2);
%define the number of echos
%n=1000;%回波个数
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    %%驰豫tau时刻+pulse length
    %off resonance的时间是 tp2+tau（180°脉冲加上半回波间隔，其影响用RTZZ表示）
    M(:,i)=RTZZ*A*M(:,i-1)+B;
   %%加上y方向flip2脉冲，绕有效轴旋转
    M(:,i)=R2*M(:,i);
    %%%驰豫tau时刻，记录下值
    %off resonance的时间是 tau（半回波间隔，其影响用RTZZZ表示）
    M(:,i)=RTZZZ*A2*M(:,i)+B2;
end
