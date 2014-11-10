clear all

%       material ==>    'Ag'  = silver
%                       'Al'  = aluminum
%                       'Au'  = gold
%                       'Cu'  = copper
%                       'Cr'  = chromium
%                       'Ni'  = nickel
%                       'W'   = tungsten
%                       'Ti'  = titanium
%                       'Be'  = beryllium
%                       'Pd'  = palladium
%                       'Pt'  = platinum
%                       'H2O' = pure water (triply distilled)

lambda=(0.2:0.01:2.8);

epsilon=LD(lambda.*1e-6,'Au','LD');

number_of_poles=6; % 1-Drude, 2-5(6)-Lorentz-Drude
damping=1; % 1-ano, 2-ne
imag_sign=1; % 1-imag=+, 2-imag=-
[real_epsilon,imag_epsilon,n]=new_LD(lambda,'Au',number_of_poles,damping,imag_sign);


plot(lambda,real(epsilon));
figure
plot(lambda,real_epsilon);