function Res = MaxCosSin(Phi)
% Determining maximum between abs(cos(Phi)) and abs(sin(Phi))
% Phi is matrix in general

Res = abs(cos(Phi)).*(abs(cos(Phi))>=cos(pi/4)) + abs(sin(Phi)).*(abs(cos(Phi))<cos(pi/4));
