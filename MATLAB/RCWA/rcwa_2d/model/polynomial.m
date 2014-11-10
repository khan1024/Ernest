%c=[1.5,0.2,0,0,0,0,0,0,0,0,-0.1,0.095,-0.1,0,0,0,0,0,0,0;];
%lambda=0.8;
%j=sqrt(-1);

function [n_model]=polynomial(lambda,c)
% kontrola velikosti

kontrola=size(c);
if kontrola(2)==20
else
    error('spatna velikost')
end

soucet_n=0;
for i=0:1:9
    soucet_n=soucet_n+c(i+1)*lambda^(i);
end
soucet_n=abs(soucet_n);

soucet_k=0;
for i=0:1:9
    soucet_k=soucet_k+c(10+i+1)*lambda^(i);
end
soucet_k=abs(soucet_k);


n_model=soucet_n-soucet_k*j;