%c=[2.27189,-0.0101081,0.0105925,0.00020817,-7.64725e-006,4.24099e-007];
%lambda=0.8;

function [n_model]=schott(lambda,c)
% kontrola velikosti
kontrola=size(c);
if kontrola(2)==6
else
    error('spatna velikost')
end

soucet=c(2-1)+c(3-1)*lambda^2+(c(4-1)/lambda^2)+(c(5-1)/lambda^4)+(c(6-1)/lambda^6)+(c(7-1)/lambda^8);
n_model=soucet^0.5;