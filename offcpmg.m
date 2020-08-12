function M=offcpmg(n,T2,T1,tau,dw,w1,flip1,flip2,phi1,phi2,ispulse)
%%%gexinmin2002 in usa 2019/06/11 ������ȷ���󣡣�
%%����off resonance�����CPMG����������ֵģ��,���ǵ���ͬ������⵹�ǣ��Լ�����
%n number of echos Ĭ��Ϊ4096
%T1 �����ԥʱ��
%T2 �����ԥʱ��
%dw ��ʾz����� delta(B0)��off resonance Ƶ�� detla(w0)��dw=2*pi*(df),��λ����ת����Hz
%w1 ��ʾ��Ƶ����Ƶ��,w=2*pi*f����λ����ת����Hz%  ��Ƶ�����Ƶ�ʣ�
%%����ƫ��Ƶ�ʣ�z�����dB���£�w=2*pi*f,�����������dw��w1�������ֵ)
%phi1 x�����λ����ʾxoyƽ���д�x������ļнǣ�һ��Ϊ0;
%phi2 y�����λ����ʾxoyƽ���д�y������ļнǣ�һ��Ϊ0;
%tau ��ز���� tau=TE/2 ��Ĭ��Ϊ0.1ms
%flip1 x�������Ƕȣ������CPMG��һ��Ϊ90�㣬flip1=pi;
%flip2 y�������Ƕȣ������CPMG,һ��Ϊ180�㣬flip=2*pi
%ispulse �Ƿ�������Ӱ�죬ispulse=0 �����ǣ�=1���ǣ�Ĭ�ϲ�����
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

%{������ʼ�Ż�ʸ��M_0�������һ��Ϊ1��,ʱ����Ϊtau=0.1�ĳ�ԥ����T��ʱ����Ϊ2*tau=0.2�ĳ�ԥ����T_2����X���ת90���õ��ĴŻ�ʸ��M_1,A����Y����ת180�����ת����%}
M0=[0 0 1]';
% ���ú�������c��ԥ���󣬿��ǵ������Ӱ��
%��ʼƽ��״̬
%{�о�DeltaB0=10%*B0ʱ���������Ǿ��ȳ���Ӱ��%}
%w1=1;%��Ƶ�����Ƶ�ʣ���
%dw=0.1*w1;%����ƫ��Ƶ�ʣ�z�����dB���£�w=2*pi*f,�����������dw��w1�������ֵ)
w=sqrt(dw.^2+w1.^2);% effective frequence of Bef ������Ч�ų���Ƶ��
%x����ʩ��90����Ƶ���� theta=90
%x����λ��Ϊ0��phi=0;
tp1=flip1/w1;%x����ת������
tp2=flip2/w1;%y����ת�������ٶ���Ƶ����ǿ�Ⱥ�Ƶ�ʲ���B1 w1
%%%���� tp1 tp2�����Լ����壡
% ����x-flip1������󾭹�tau+pl1�ĳ�ԥ����
[A1,B1] = freeprecess(tau+ispulse*tp1,T1,T2);
% ����y-flip2������󾭹�tau+pl1�ĳ�ԥ����
[A2,B2] = freeprecess(tau+ispulse*tp2,T1,T2);
%���ú���
R1=rotxn(flip1,phi1,w1,dw);
M0=R1*M0;
%z����ת����
%��Ч��ת�Ƕ�
gama1=dw*(tau+ispulse*(tp1));
RTZ1=[cos(gama1) sin(gama1) 0; -sin(gama1) cos(gama1) 0; 0 0 1];
%��Ч��ת�Ƕ�
gama2=dw*(tau+ispulse*(tp2));
RTZ2=[cos(gama2) sin(gama2) 0; -sin(gama2) cos(gama2) 0; 0 0 1];
gama=dw*tau;
RTZ=[cos(gama) sin(gama) 0; -sin(gama) cos(gama) 0; 0 0 1];
% ��ԥ��tauʱ��+pulse length
M0=RTZ*A1*M0+B1;
%ʩ��y�����flip2��Ƶ���� ��λΪy
R2=rotyn(flip2,phi2,w1,dw);
M0=R2*M0;
%%%����tau��ʼ��ԥ����¼��һ���ز�����
M1=RTZ2*A2*M0+B2;
%%��ʼѭ��
[A,B] = freeprecess(tau,T1,T2);
%define the number of echos
%n=1000;%�ز�����
M=zeros(3,n);
M(:,1)=M1;
for i=2:n
    %%��ԥtauʱ��+pulse length
    M(:,i)=RTZ*A*M(:,i-1)+B;
    %%����y����flip2���壬����Ч����ת
    M(:,i)=R2*M(:,i);
    %%%��ԥtauʱ�̣���¼��ֵ
    M(:,i)=RTZ2*A2*M(:,i)+B2;
end
