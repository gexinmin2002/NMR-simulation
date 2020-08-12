function E=lognorm(a,x,y);
%sum of the square residual during log normal nonlinear regression
%x是时间 0.1 0.2 0.3等
%y是对应的驰豫谱幅度
%a是初始值
load t2_and_f_water.txt %纯水的信号
load t_and_g.dat
t2_wat=t2_and_f_water(:,1);%水的T2布点数
f_wat=t2_and_f_water(:,2);%水的分布谱幅度
Mo=t2_and_g(1,2);
n=length(x);%number of CPMG signal points for fitting
f0=Mo-sum(f_wat); %sum of f for bitumen part 
S=a(1);%signal in lognormal distribution
M=a(2)%Mu in lognormal distribution
Ln_t2_bit(1)=M-2.5*S;
%Lognormal distribution fitting
for i=2:11
    Ln_t2_bit(i)=Ln_t2_bit(i-1)+S/2;
    g(i)=1/((S+eps)*sqrt(2*pi))*exp(-(Ln_t2_bit(i)-M)^2/(2*(S+eps)^2))*abs(S/2);
    t2_bit(i)=exp(Ln_t2_bit(i));
end
g=g';
t2_bit=t2_bit';
for j=1:n
    %bitumen part
    for k=1:length(t2_bit)
        amp_bit(k)=f0*g(k)*exp(-x(j)/(t2_bit(k)+eps));
    end
    sig_bit=sum(amp_bit);
    %water part
    for kk=1:length(t2_wat);
        amp_wat(kk)=f_wat(kk)*exp(-x(j)/t2_wat(kk));
    end
    sig_wat=sum(amp_wat);
    %CPMG signal=bitumen part+water part
    h(j,1)=sig_bit+sig_wat;
end
E=sum((y-h)^2);



