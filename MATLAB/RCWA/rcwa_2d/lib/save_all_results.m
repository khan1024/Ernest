% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

record_p_0=[cycle; D_R_0; D_T_0];
record_0=[record_0,record_p_0];

record_p_90=[cycle; D_R_90; D_T_90];
record_90=[record_90,record_p_90];

record_p_psi=[cycle; D_R; D_T];
record_psi=[record_psi,record_p_psi];

savefile=[save_m_as '.mat'];
save(savefile, 'record_0','record_90','record_psi','N_X','N_Y')

if grating==0
    %for i=1:1:number_of_layers
    %    save(savefile,'-append', 'eval(['coordinate_x_' num2str(i)])','eval(['coordinate_x_' num2str(i)])')
    %end
else
    save(savefile,'-append', 'a_ellipse','b_ellipse','n_layer','n_cylinder')
end

save(savefile,'-append', 'Lambda_x','Lambda_y','thickness');