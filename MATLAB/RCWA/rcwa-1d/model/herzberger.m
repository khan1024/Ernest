%c=[1,-0.00309946,-9.61396e-006,1.72005,0,0.00561194,0,0,0,0,0.028,0,0,0,0,-1.09862e-005];
%lambda=1;

function [n_model]=herzberger(lambda,c)
% kontrola velikosti
kontrola=size(c);
if kontrola(2)==16
else
    error('spatna velikost')
end

soucet=0;
for i=0:1:4
    soucet=soucet+(c(7+i-1))/(lambda^2-c(12+i-1))+(c(17-1))/(lambda^2-0.0028);
end

soucet=soucet+c(3-1)*lambda^3+c(5-1)+c(4-1)*lambda^4+c(6-1)/lambda^2;
n_model=soucet^c(2-1);