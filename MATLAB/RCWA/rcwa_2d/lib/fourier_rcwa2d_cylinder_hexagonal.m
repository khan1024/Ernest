% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

coefficient_E=0;
coefficient_A=0;
subtotal_E=0;
subtotal_A=0;
subsubtotal_E=0;
subsubtotal_A=0;

rho=sqrt((coeff_i/(b_ellipse*Lambda_x))^2+(coeff_j/(a_ellipse*Lambda_y))^2);
if coeff_i==0 && coeff_j==0         
    subsubtotal_E=epsilon_1+(epsilon_2-epsilon_1)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y);
    subsubtotal_A=1/epsilon_1+(1/epsilon_2-1/epsilon_1)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y);
    subtotal_E=subsubtotal_E+(epsilon_1d+(epsilon_2d-epsilon_1d)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y))*exp(-1i*2*pi*0.5*(coeff_i+coeff_j));
    subtotal_A=subsubtotal_A+(epsilon_1i+(epsilon_2i-epsilon_1i)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y))*exp(-1i*2*pi*0.5*(coeff_i+coeff_j));
else
    besselka=besselj(1,2*pi*rho*a_ellipse*b_ellipse);
    subsubtotal_E=(epsilon_2-epsilon_1)*besselka/(rho*Lambda_x*Lambda_y);
    subsubtotal_A=(1/epsilon_2-1/epsilon_1)*besselka/(rho*Lambda_x*Lambda_y);
    subtotal_E=subsubtotal_E+(epsilon_2d-epsilon_1d)*besselka/(rho*Lambda_x*Lambda_y)*exp(-1i*2*pi*0.5*(coeff_i+coeff_j));
    subtotal_A=subsubtotal_A+(epsilon_2i-epsilon_1i)*besselka/(rho*Lambda_x*Lambda_y)*exp(-1i*2*pi*0.5*(coeff_i+coeff_j));
end

coefficient_E=coefficient_E+subtotal_E;    
coefficient_A=coefficient_A+subtotal_A;
