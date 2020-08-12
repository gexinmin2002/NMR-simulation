%%simulate the infleunce of off resonance
%%%gexinmin2002 in usa 2019/06/11 ������ȷ���󣡣�
%{������ʼ�Ż�ʸ��M_0�������һ��Ϊ1��,ʱ����Ϊtau=0.1�ĳ�ԥ����T��ʱ����Ϊ2*tau=0.2�ĳ�ԥ����T_2����X���ת90���õ��ĴŻ�ʸ��M_1,A����Y����ת180�����ת����%}
clc;
clear;
tau = 0.05;		% ��ز����.
T1=150;%ms
T2=100;%ms
% ���ú�������
[A,B] = freeprecess(tau,T1,T2);
%��ʼƽ��״̬
M0=[0 0 1]';
w1=1;%��Ƶ�����Ƶ�ʣ���
dw=0.5*w1;%����ƫ��Ƶ�ʣ�z�����dB���£�w=2*pi*f,�����������dw��w1�������ֵ)
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef
%x����ʩ��90����Ƶ���� theta=90
%x����λ��Ϊ0��phi=0;
phi1=0;
flip=pi/2;
t90=flip/w1;
%���ú���
R=rotxn(flip,phi1,w1,dw);
M0=R*M0;
%z����ת����
gama=dw*tau;
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];
% ��ԥ��tauʱ��
M0=RTZ*A*M0+B;
%ʩ��y�����180����Ƶ���� ��λΪy
flip2=pi;
t180=flip/w1;
phi2=0;%y����λ��Ϊ0
R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
%%%����tau��ʼ��ԥ����¼��һ���ز�����
M1=RTZ*A*M0+B;
%%��ʼѭ��
%define the number of echos
n=1000;%�ز�����
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    %%��ԥtauʱ��
    M(:,i)=RTZ*A*M(:,i-1)+B;
    %%����y����180�����壬����Ч����ת
    M(:,i)=R2*M(:,i);
    %%%��ԥtauʱ�̣���¼��ֵ
    M(:,i)=RTZ*A*M(:,i)+B;
end
M=M';