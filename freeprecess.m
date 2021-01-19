function [Afp,Bfp]=freeprecess(T,T1,T2)

E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1];
Bfp = [0 0 1-E1]';

