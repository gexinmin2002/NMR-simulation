%%%%%mlinreg.m
function [b,E,Stdb,Stdfit] = mlinreg(y,Vin,W)
%一种新的优化方法
% Weighted multilinear regression
% y = Vin * b + E
% (Jx1) (JxP) (Px1) (Jx1)
% where W (Jx1) is the weight of each data point,
% ordinary least squares: use W = ones(size(y)); 权重
%
% b (Px1) = fitting coefficients 反演参数
% E (Jx1) = absolute error in each fitting point 测量误差
% Stdb (Px1) = standard deviation in percent of fitting coefficients 标准差
% Stdfit (1x1) = standard deviation in percent of fit
% beware that this uses the relative
J = length(y);
P = length(Vin(1,:));
V = (W*ones(1,P) ).*Vin;
VTV = inv(V'*V);
b = VTV*V'*(W.* y );
E = y - Vin*b;
RelE = E./y;
Var = (E'*E)/(J-P);
Stdb = 100.*sqrt(Var.*diag(VTV))./abs(b);
Stdfit = 100.*sqrt((RelE'*RelE)./(J-P));
%in page 132 Dynamic Dynamic Pulsed-Field-Gradient NMR