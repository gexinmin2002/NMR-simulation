%function for the lognornal distribution fitting for bitumen by Elton Yang
clear;
clc;
clf;
load t_and_g dat
%第一个点存储Mo信号（也就是FID信号）

x=t_and_g(2:end,1);%采集时间，不加FID信号
y=t_and_g(2:end,2);%信号幅度，不加FID信号
y=abs(y);%取模
Mo=t_and_g(1,2);%取第一个点
a0=[2,0];%initial guess %a0=[S,0] M=ln(T2),即有T2=exp(M），M等于0,T2=1ms;
%S是标准差，初始为2
options=optimst('TolX',0.000001);%最优化函数
options.Display='off';%？
a=fminsearch('lognorm',a0,options,x,y);%调用最优化函数
%For plotting
load t2_and_f_water.txt%纯水的T2谱
t2_wat=t2_and_f_water(:,1);%纯水的横向驰豫时间分布
f_water=t2_and_f_water(:,2);%纯水的驰豫时间对应幅度
n=length(x) %number of CPMG signal points for fitting %CPMG信号采样点数
f0=Mo-sum(f_wat); %sum of f for bitumen part%沥青含量
S=a(1);%sigma in lognormal distribution
M=a(2);%Mu in lognornal distribution
Ln_t2_bit(1)=M-2.5*S;
%lognormal distribution fitting
for i=2:11
    Ln_t2_bit(i)=Ln_t2_bit(i-1)+S/2;
    g(i)=1/(S*sqrt(2*pi)*exp(-Ln_t2_bit(i)-M)^2/(2*S^2))*abs(S/2);
    t2_bit(i)=exp(Ln_t2_bit(i));
end
g=g';
t2_bit=t2_bit';
%generate new x for plotting
spacing=t_and_g(3,1);
add_x=x(1):0.01:x(round(3/spacing));
x_plot=[add_x',x((round(3/spacing)+1)end)];
n=length(x_plot);
for j=1:n
    %bitumen part
    for k=1:length(t2_bit)
        amp_bit(k)=f0*g(k)*exp(-x_plot(j)/t2_bit(k)+eps));
    end
    sig_bit=sum(amp_bit);
    %water part
    for kk=1:length(t2_wat)
        amp_wat(kk)=f_wat(kk)*exp(-x_plot(j)/t2_wat(kk));
    end
    sig_wat=sum(amp_wat);
    %CPMG signal=bitumen part+water part
    h(j,1)=sig_bit+sig_wat;
end
E=lognorm(a,x,y);
semilogx(x,y,'b+');
hold on
semilogx(x_plot,h,'r');
axis([0,20,0,20]);
xlabel('Time(msec)');
ylabel('Fitting curve');
fit_dat=[x_plot h];


    
    

