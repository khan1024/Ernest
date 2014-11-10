% Function I_mat
% This function calculates the transfer matrix, I, for reflection and
% transmission at an interface between materials with complex dielectric 
% constant n1 and n2.

function I = I_mat_s_polarized(n_1,n_2,theta)
% r=(n_1-n_2)/(n_1+n_2);
% t=2*n_1/(n_1+n_2);
costh1 = cos(asin((sin(theta))/n_1));
costh2 = cos(asin((sin(theta))/n_2));

r=(n_2*costh2-n_1*costh1)/(n_2*costh2+n_1*costh1);
t=1+r;
I=[1 r; r 1]/t;