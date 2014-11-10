% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function [alpha,beta]=alpha_beta_mex(n_x,alpha_m,beta_n)

n_sq=n_x*n_x;

alpha=zeros(n_sq,n_sq);
beta=zeros(n_sq,n_sq);

for m=1:1:n_x
    for n=1:1:n_x
        for j=1:1:n_x
            for l=1:1:n_x
                u=n_x*(m-1)+n;
                v=n_x*(j-1)+l;                
                    
                if m==j && n==l
                    alpha(u,v)=alpha_m(m);
                    beta(u,v)=beta_n(n);
                end
            end
        end
    end
end