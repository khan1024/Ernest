% To calculate the angle from normal of the light ray travel in a waveguide where the light is coupled in from a prism evanescently
% 
% n0 is the index of refraction of the surrounding medium, which is normally 1
% nP is the index of the prism
% nF is the index of the thin film waveguide (or whichever kind of waveguide)
% theta0 is the incoming angle of light, normal to the prism
% theta is the bottom angle of the prism, with respect to parallel to the waveguide
% theta1 is the refracted angle in the prism, relative to the normal of the slanted face of the prism
% thetaP0 is the angle that makes a triangle with theta and (90-theta1)
% thetaP is the angle of light "exiting" the prism, with resped to normal to the bottom of the prism and the waveguide
% thetaF is the angle of light entering the film
% m is the set of allowed modes in the waveguide
% d is the thickness of the waveguide
% lambda0 is the incoming wavelength of light
% N is the maximum number of modes allowed in the waveguide

% Created by Dr. Katherine Han 7/11/14
clear all

lambda0 = 1000; % incoming wavelength in nm
d = 1500; % thickness of waveguide
nF = 2; % index of the waveguide
nP = 2.5; %index of the prism
n0 = 1; %index of air
n3 = 1.5; %highest index of either cladding material, such as a substrate for the waveguide
N = 2 * d / lambda0 * (nF^2-n3^2)^0.5; % maximum number of allowed modes in the waveguide
Nfloor = floor(N);
m = 0:1:Nfloor-1;
thetaF = asin(m*lambda0/2/d/nF);

for z = 2:Nfloor
    thetaP (z-1) = asin(nF/nP*sin(thetaF(z))); %calculates outgoing prism angles to normal of the waveguide for allowed modes
end

thetaP0 = pi / 2 - thetaP;

syms theta
for n = 1:length(thetaP)

solve((pi/2 - asin(n0/nP*sin(pi/2-theta)))+theta+thetaP0(n)-pi==0)
theta_options(n) = ans;
end

theta_options_deg = theta_options .* (180 / pi)
