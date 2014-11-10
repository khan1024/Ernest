% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% pruchod celou strukturou

p_eigvec_E=eigvec_E;
p_eigvec_H=eigvec_H;

for i=1:1:number_of_layers
    current_gamma=eval(['gamma_2_' num2str(i)]);
    eval(['p_gamma_' num2str(i), '=current_gamma;']);
end

number_of_layers=original_number_of_layers;

%gamma=zeros(2*n_sq,2*n_sq,number_of_layers);
eigvec_E=zeros(2*n_sq,2*n_sq,number_of_layers);
eigvec_H=zeros(2*n_sq,2*n_sq,number_of_layers);

% zpětné přiřazení vypočtených matic
for i=1:1:number_of_layers
    eigvec_E(:,:,i)=p_eigvec_E(:,:,matice_stejnosti_final(i,3));
    eigvec_H(:,:,i)=p_eigvec_H(:,:,matice_stejnosti_final(i,3));
    
    current_gamma=eval(['p_gamma_' num2str(matice_stejnosti_final(i,3))]);
    eval(['gamma_2_' num2str(i), '=current_gamma;']);
end

clear p_eigvec_E p_eigvec_H p_gamma
