function matrix=harmonics_reconstruction_x_y_2(N_X,comp,eigvec)

% 1 druhá složka = N_X členů E_0_j
% 2 třetí složka = N_X členů E_i_0
% 3 čtvrtá složka = N_X*N_X Ex_i_j
% 4 pátá složeka = N_X*N_X Ey_i_j

if comp==1
    E_harmonics=eigvec(1:N_X);
    matrix=conj(E_harmonics)';
elseif comp==2
    E_harmonics=eigvec(N_X+1:2*N_X);
    matrix=E_harmonics;
elseif comp==3
    E_harmonics=eigvec(2*N_X+1:end-N_X*N_X);
    matrix=reshape(E_harmonics',N_X,N_X)';
elseif comp==4
    E_harmonics=eigvec(end-N_X*N_X+1:end);
    matrix=reshape(E_harmonics',N_X,N_X)';
end