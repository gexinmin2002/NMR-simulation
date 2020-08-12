function [Afp,Bfp]=freeprecess(T,T1,T2)
%%程序可用
%%%gexinmin usa 2019/06/11
%	 函数计算进动和时间间隔T的衰变
%    给出T1，T2和偏共振频率的值
%	 no时间单位是ms, 偏共振频率单位是 Hz.

E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1];
Bfp = [0 0 1-E1]';

