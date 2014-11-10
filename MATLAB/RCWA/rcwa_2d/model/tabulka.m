%clear all
%format long
%lambda=1.0;
%nahraj_material='gold.txt';
%j=sqrt(-1);

function [n_model]=tabulka(lambda,nahraj_material)

filename=[nahraj_material];
data=dlmread(filename);

velikost=size(data);

% kontrola zadání
if lambda>=data(1,1) & lambda<=data(velikost(1),1)
else
    error('mimo rozsah')
end

% musím nalézt nejbližšší větší prvek, poté cyklus končí
% pokud jsem se hodnotou trefil do hodnoty v tabulce, přečtou se další data a výpočet končí
for i=1:1:velikost(1)
    if data(i,1)>lambda
        index=i;
        zasah=1;
        break
    elseif data(i,1)==lambda
        n_model=data(i,2)-data(i,3)*j;
        zasah=2;
        break
    end
end

% mezi dvěma známými body se hodnoty linearizují -> používám rovnici přímky y=kx+b
% -> musím zjistit ze dvou známých bodů konstanty k, b
if zasah==1
    k_r=(data(index,2)-data(index-1,2))/(data(index,1)-data(index-1,1));
    k_i=(data(index,3)-data(index-1,3))/(data(index,1)-data(index-1,1));
    b_r=(data(index-1,1)*data(index,2)-data(index,1)*data(index-1,2))/(data(index-1,1)-data(index,1));
    b_i=(data(index-1,1)*data(index,3)-data(index,1)*data(index-1,3))/(data(index-1,1)-data(index,1));
    n_r=k_r*lambda+b_r;
    n_i=k_i*lambda+b_i;
    n_model=n_r-n_i*j;
end