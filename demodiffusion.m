%%反演扩散系数 
%gexinmin 2017/07/16 Houston finished

%load into Matlab ascii file called data
%This containing one row, the measured intensities as a function of
%incrementing gradient strength
nbl=size(data);
I=data(1:nbl);
%I=[1.6044e+09 1.4715e+09 1.3535e+09 1.1901e+09 1.1901e+09 8.1347e+08
%6.3767e+08 4.8984e+08 3.6859e+08 2.7482e+08 2.0493e+08 1.7096e+08
%1.3426e+08 9.9259e+07 8.1729e+07 6.9994e+07]';
gam=7.16e+8; % the square of the gyro magnetic ratio for the proton 旋磁比的平方
 G=input('What is the gradient strength in G/cm A ? ');%G=25;
 Istep =input('What is the minimum current on your first gradient pulse ? ');%Istep = 2;
 Step =input('What is the step current in the experiment ? '); %step= 2;
 a1 =input('What is the Z-storage time in milliseconds ? '); %a1=250;
 a2 =input('What is the Tau value in milliseconds ? '); %a2 = 4;
 delta =input('What is the length of the gradient pulse in milliseconds ? '); %delta = 3;
 %计算梯度 predict the gradient 
 for I=1:nbl
     x1(i) = (i-1)*step+Istep; % calculating the applied gradient strengths
 end
 x=x1';
 %change the time of ms to s
 a1=a1*1e-3;%换成s
 a2=a2*1e-3;%换成s
 delta=delta*1e-3;%换成s
 B=-(a1+1.5*a2-delta/6)*(x.^2)*(G).^2*ga*4*delta.^2;
 figure(1)
 plot(x.^2,I,'r+');
 hold on
 Y=log(I); % Linearzation to the model Y=A*B
 A=[ones(nbl,1),B];
 %weighting of the data
 %W=ones(size(Y)); %equal weight
 W=I.^3; %weight propotional to the some power of the data-intensdities
 %call of the function mlinreg
 [b,E,Stdb,Stdfit]=mlinreg(Y,A,W);
 Y_modell=A*b;
 I_0=exp(b(1)); %Fitted inintial intensity
 Diffusion_coefficient=(b(2)); %fitted diffusion coefficient
 I_modell=exp(Y_modell); % plot of exponential attenuation
 plot(x.^2,I_modell,'b-');
 hold off;
 
 figure(2)
 plot(x.^2,log(I/I_null),'r+'); %plot of the linear attenuation
 hold on;
 plot(x.^2,log(I_modell/I_null),'b-');
 hold off;
 