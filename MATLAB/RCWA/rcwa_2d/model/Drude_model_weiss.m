function epsM = Drude_model_weiss(lambda)

c = 299792458; % light velocity in vacuum
omega = 2*pi*c./lambda; % light frequency
omegaP = 1.37e16; % plasma frequency
gamma = 0.85e14; % colission frequency
eps_inf = 1; % high-frequency permitivity

epsM = eps_inf - omegaP^2./(omega.^2+1i*omega*gamma);
end
