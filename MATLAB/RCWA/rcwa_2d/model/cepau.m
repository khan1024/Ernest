function [n_model]=cepau(lambda)

ni=1/(1e-4*lambda);
eps0=0.31;
OmegP0=71955.95;
OmegTau0=340.11;
OmegOO1=39164.40;
OmegP1=81797.21;
OmegTau1=0.;
OmegOO2=21861.01;
OmegP2=15793.08;
OmegTau2=3851.46;
OmegOO3=25104.10;
OmegP3=27053.08;
OmegTau3=5481.97;
% evaluation
epsr = eps0-OmegP0^2/(ni^2+OmegTau0)+(OmegP1^2*(OmegOO1^2-ni^2))/...
  ((OmegOO1^2-ni^2)^2+(OmegTau1*ni)^2)+(OmegP2^2*(OmegOO2^2-ni^2))/...
  ((OmegOO2^2-ni^2)^2+(OmegTau2*ni)^2)+(OmegP3^2*(OmegOO3^2-ni^2))/...
  ((OmegOO3^2-ni^2)^2+(OmegTau3*ni)^2); 
epsi = OmegTau0*OmegP0^2*(ni^-1)/(ni^2+OmegTau0^2)+OmegP1^2*OmegTau1*ni/...
  ((OmegOO1^2-ni^2)^2+(OmegTau1*ni)^2)+OmegP2^2*OmegTau2*ni/...
  ((OmegOO2^2-ni^2)^2+(OmegTau2*ni)^2)+OmegP3^2*OmegTau3*ni/...
  ((OmegOO3^2-ni^2)^2+(OmegTau3*ni)^2);
qn = epsr + i*epsi;
n_model=conj(sqrt(qn));