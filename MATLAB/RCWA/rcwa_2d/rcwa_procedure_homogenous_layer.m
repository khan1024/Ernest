% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

current_gamma=zeros(n_sq,n_sq);
current_gamma=sqrt(I_mat*(2*pi/lambda*current_hom_n(layer))^2-sec(zeta)^2*(alpha^2+beta^2-2*alpha*beta*sin(zeta)));
for i=1:1:n_sq
    %current_gamma(i,i)=sqrt((2*pi/lambda*current_hom_n(layer))^2-sec(zeta)^2*(alpha(i,i)^2+beta(i,i)^2-2*alpha(i,i)*beta(i,i)*sin(zeta)));
    
    % change of signs, equation 16 [1]
    
    % if imag(gamma)<0 (? nektere clanky maji obracene)
    if (imag(current_gamma(i,i))>0)
        current_gamma(i,i)=-1*current_gamma(i,i);
    end
    
    % vlastni uprava --- equation 16 [1]
    if abs(imag(current_gamma(i,i)))<abs(real(current_gamma(i,i)))
        if (real(current_gamma(i,i))<0)
            current_gamma(i,i)=-1*current_gamma(i,i);
        end
    end
end 

A_mn_h=sparse((I_mat*(k_0*current_hom_n(layer))^2-alpha*alpha)/(mu*k_0*cos(zeta)*current_gamma));
B_mn_h=sparse((I_mat*(k_0*current_hom_n(layer))^2-beta*beta)/(mu*k_0*cos(zeta)*current_gamma));
C_mn_h=sparse((alpha*beta-I_mat*(k_0*current_hom_n(layer))^2*sin(zeta))/(mu*k_0*cos(zeta)*current_gamma));

eigvec_E(:,:,layer)=sparse([I_mat zero ; zero I_mat]);
eigvec_H(:,:,layer)=[-C_mn_h, -A_mn_h; B_mn_h, C_mn_h];

% ?
current_gamma=diag([diag(current_gamma);diag(current_gamma)]);

eval(['gamma_2_' num2str(layer), '=sparse(current_gamma);']);