%{
clear all
pocet_vrstev=1;
sirka_celek=1;
sirka=sirka_celek*1E-6;
l=1;
pocet_radu=201;
M=pocet_radu-1;
data=[1,-5,3,-1,1,1,3,5];
d=10;
%}

function [toeplitz_matrix]=function_toeplitz_matrix(data,pocet_radu,d)

l=1;
M=pocet_radu-1;

velikost_dat=size(data);
velikost_dat_x=velikost_dat(2);


index=1:1:velikost_dat_x; % pomocné pole

% kde se nachází NaN
% !!!!! v aktuální vrstvě se NaN nemusí nacházet

polohaNaN=0;
for i=2:2:velikost_dat_x
    if isnan(data(l,index(i)))==1
        polohaNaN=i;
        posledni_eps=i-2;
        break % najde první NaN
    end
end

if polohaNaN==0
    posledni_eps=velikost_dat_x-1;
end

fourier_coefficients=0*[-M:M];

for i=2:2:posledni_eps
    x_1_fourier=data(l,index(i-1));
    x_2_fourier=data(l,index(i+1));
    eps_fourier=data(l,index(i));
    a_fourier=x_1_fourier+x_2_fourier;
    
    fourier_coefficients=fourier_coefficients+2.*eps_fourier.*(x_2_fourier-1/2.*a_fourier).*sinc(2.*[-M:M].*(x_2_fourier-...
        1/2.*a_fourier)./d)./d./exp(1i.*pi.*[-M:M].*a_fourier./d);
end
toeplitz_matrix=toeplitz(fliplr(fourier_coefficients(1:pocet_radu)),fourier_coefficients(pocet_radu:2*pocet_radu-1));

