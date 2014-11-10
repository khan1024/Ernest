% Function L_mat
% This function calculates the propagation matrix, L, through a material of
% complex dielectric constant n and thickness d for the wavelength lambda.
function L = L_mat_aoi(n,d,lambda,theta)
xi=2*pi*n/lambda*d;
costh=cos(asin(1*sin(theta)/n));
L=[exp(-1i*xi*costh) 0; 0 exp(1i*xi*costh)];