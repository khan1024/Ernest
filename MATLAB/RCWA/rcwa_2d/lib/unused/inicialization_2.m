% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%if index_zrychleni==1
E=zeros(n_sq,n_sq);
A=zeros(n_sq,n_sq);

inv_E=zeros(n_sq,n_sq);
inv_A=zeros(n_sq,n_sq);
%end

A_d=zeros(n_sq,n_sq);
B=zeros(n_sq,n_sq);
C=zeros(n_sq,n_sq);
D=zeros(n_sq,n_sq);

inv_B=zeros(n_sq,n_sq);
inv_D=zeros(n_sq,n_sq);

G_1=zeros(n_sq,n_sq);
G_2=zeros(n_sq,n_sq);

G=zeros(n_sq,n_sq);
F=zeros(n_sq,n_sq);

inv_G=zeros(n_sq,n_sq);
inv_F=zeros(n_sq,n_sq);

if maticovy_algoritmus==1
    Q=zeros(2*n_sq, 2*n_sq, number_of_layers);
    W_1=zeros(n_sq, 2*n_sq, number_of_layers);
    W_2=zeros(n_sq, 2*n_sq, number_of_layers);
    V_1=zeros(n_sq, 2*n_sq, number_of_layers);
    V_2=zeros(n_sq, 2*n_sq, number_of_layers);
else
    Q=zeros(2*n_sq,2*n_sq);
    W_1=zeros(n_sq,2*n_sq);
    W_2=zeros(n_sq,2*n_sq);
    V_1=zeros(n_sq,2*n_sq);
    V_2=zeros(n_sq,2*n_sq);
end