% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic


coefficient_E=0;
coefficient_A=0;

subtotal_E=0;
subtotal_A=0;

if measurement==2
    
end

rho=sqrt((coeff_i/(b_ellipse*Lambda_x))^2+(coeff_j/(a_ellipse*Lambda_y))^2);
if coeff_i==0 && coeff_j==0         
    subtotal_E=epsilon_1+(epsilon_2-epsilon_1)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y);
    subtotal_A=1/epsilon_1+(1/epsilon_2-1/epsilon_1)*pi*a_ellipse*b_ellipse/(Lambda_x*Lambda_y);
else
    besselka=besselj(1,2*pi*rho*a_ellipse*b_ellipse)/(rho*Lambda_x*Lambda_y);
    subtotal_E=(epsilon_2-epsilon_1)*besselka;
    subtotal_A=(1/epsilon_2-1/epsilon_1)*besselka;
end

coefficient_E=coefficient_E+subtotal_E;    
coefficient_A=coefficient_A+subtotal_A;