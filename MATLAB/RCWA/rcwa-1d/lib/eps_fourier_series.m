% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% see Moharam et al., Formulation for stable and efficient implementation
% of the rigorous coupled-wave analysis of binary gratings, JOSA A, 12, p.
% 1068, 1995

% equation (2)
% zero fourier coefficient
% eps_0=n^2_{gr}(1-h)+n^2_{rd}h
% eps_n=\left(n^2_{rd}-n^2_{gr} \right)\frac{\sin(\pi mh)}{\pi m}

% zero coefficient
epsG=(1-duty_cycle(l))*epsg+duty_cycle(l)*epsr;
i_epsG=(1-duty_cycle(l))/epsg + duty_cycle(l)/epsr;

% non-zero coefficients
Sinc=sin(pi*duty_cycle(l)*(1:M))./(pi*(1:M));   
v_m=(epsr-epsg)*fliplr(Sinc);    
v_0=epsG;
v_p=(epsr-epsg)*Sinc;

% eqution (2) holds for -Lambda/2 to Lambda/2 -> move to 0 to Lambda    
v=[v_m v_0 v_p].*exp(-1i*2*pi*shift(l)*(-M:M));

    
% same procedure for inverse coefficients
i_vm=(1/epsr-1/epsg)*fliplr(Sinc);    
i_v0=i_epsG;
i_vp=(1/epsr-1/epsg)*Sinc;
i_v=[i_vm i_v0 i_vp].*exp(-1i*2*pi*shift(l)*(-M:M));