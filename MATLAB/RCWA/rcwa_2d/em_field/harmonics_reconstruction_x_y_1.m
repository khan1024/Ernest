function matrix=harmonics_reconstruction_x_y_1(N_X,comp,eigvec)

% 1 prvni složka = 1 členu E_0_0
% 2 druhá složka = N_X členů E_0_j
% 3 třetí složka = N_X členů E_i_0
% 4 čtvrtá složka = N_X*N_X Ex_i_j
% 5 pátá složeka = N_X*N_X Ey_i_j

if comp==1
    E_harmonics=eigvec(1);
    matrix=E_harmonics;
elseif comp==2
    E_harmonics=eigvec(2:N_X+1);
    matrix=conj(E_harmonics)';
elseif comp==3
    E_harmonics=eigvec(N_X+2:2*N_X+1);
    matrix=E_harmonics;
elseif comp==4
    E_harmonics=eigvec(2*N_X+2:end-N_X*N_X);
    matrix=reshape(E_harmonics',N_X,N_X)';
elseif comp==5
    E_harmonics=eigvec(end-N_X*N_X+1:end);
    matrix=reshape(E_harmonics',N_X,N_X)';
end