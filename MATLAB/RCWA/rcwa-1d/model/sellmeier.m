%c=[0.5,1,1.03961,0.231792,1.01147,0,0,0.0060007,0.0200179,103.561,0,0];
%lambda=1;

function [n_model]=sellmeier(lambda,c)
% kontrola velikosti
kontrola=size(c);
if kontrola(2)==12
else
    error('spatna velikost')
end

soucet=0;
for i=0:1:4
    soucet=soucet+(c(4+i-1)*lambda^2)/(lambda^2-c(9+i-1));
end

soucet=soucet+c(3-1);
n_model=soucet^c(2-1);