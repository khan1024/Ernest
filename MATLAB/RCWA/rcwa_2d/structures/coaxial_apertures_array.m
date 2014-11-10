number_of_layers=1;

mu=1;
lambda=0.54;
n_1=1;
n_3=1.45;
theta0=0;

l=1; % layer 1
thickness(l)=0.15;

Lambda_x=0.3;
Lambda_y=Lambda_x;

w_1=0.155;
w_2=0.105;

coordinate_x_1=[0,Lambda_x/2-w_1/2,Lambda_x/2-w_2/2,Lambda_x/2+w_2/2,Lambda_x/2+w_1/2,Lambda_x];
coordinate_y_1=coordinate_x_1;

n_ag=0.0187-3.8069*1i;

r_index_1=[n_ag, n_ag, n_ag, n_ag, n_ag;
    n_ag, 1, 1, 1, n_ag;
    n_ag, 1, n_ag, 1, n_ag;
    n_ag, 1, 1, 1, n_ag;
    n_ag, n_ag, n_ag, n_ag, n_ag];