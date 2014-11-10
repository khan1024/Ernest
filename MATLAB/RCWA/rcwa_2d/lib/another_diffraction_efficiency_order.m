% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function output=another_diffraction_efficiency_order(input)

load temp.mat;

diffraction_efficiency_order(1)=input(1);
diffraction_efficiency_order(2)=input(2);

poloha_dif_ucinnost=1+(diffraction_efficiency_order(1)+N_Y)+(2*N_Y+1)*(diffraction_efficiency_order(2)+N_Y);

output=[diffraction_efficiency_order(1) diffraction_efficiency_order(2) D_R_0(poloha_dif_ucinnost) D_T_0(poloha_dif_ucinnost) sum(D_R_0) sum(D_T_0) D_R_90(poloha_dif_ucinnost) D_T_90(poloha_dif_ucinnost) sum(D_R_90) sum(D_T_90) D_R(poloha_dif_ucinnost) D_T(poloha_dif_ucinnost) sum(D_R) sum(D_T)];
fprintf(1,'order  [%5.2f,%6.2f] R_0=%12.10f T_0=%12.10f sumR=%12.10f sumT=%12.10f R_90=%12.10f T_90=%12.10f sumR=%12.10f sumT=%12.10f R_psi=%12.10f T_psi=%12.10f sumR=%12.10f sumT=%12.10f\n', output);