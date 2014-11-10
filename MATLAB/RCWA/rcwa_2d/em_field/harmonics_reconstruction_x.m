function matrix=harmonics_reconstruction_x(N_X,comp,eigvec)

n_ords=2*N_X+1;

% 1 prvni složka = (2*N_X+1) členu
% 2 druhá složka = (2*N_X+1)*N_X členů
% 3 třetí složka = (2*N_X+1)*N_X členů

if comp==1
    E_harmonics=eigvec(1:(2*N_X+1));
    matrix=conj(E_harmonics'); % ze sloupce udělán řádek
elseif comp==2
    E_harmonics=eigvec((2*N_X+1)+1:(end-(2*N_X+1)*(N_X)));
    matrix=reshape(E_harmonics',n_ords,(n_ords+1)/2-1)';
elseif comp==3
    E_harmonics=eigvec((end-(2*N_X+1)*(N_X)+1):end);
    matrix=reshape(E_harmonics',n_ords,(n_ords+1)/2-1)';
end