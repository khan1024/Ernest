function epsM = Drude_model_fun(lambda)

c = 299792458; % light velocity in vacuum
omega = 2*pi*c./lambda; % light frequency
omegaP = 1.3388e16; % plasma frequency
gamma = 7.07592e13; % colission frequency
eps_inf = 3.36174; % high-frequency permitivity

epsM = eps_inf - omegaP^2./(omega.^2+1i*omega*gamma);
epsM = conj(sqrt(epsM)); % return refractive index
end
