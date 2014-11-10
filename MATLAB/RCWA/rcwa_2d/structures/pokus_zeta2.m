number_of_layers=1;

mu=1;
lambda=1;
n_1=1.0;
n_3=1.5;
theta0=1e-5;
phi0=0;

l=1; % layer 1
thickness(l)=lambda;

Lambda_x=2*lambda;
Lambda_y=Lambda_x;


coordinate_x_1=[0,lambda,Lambda_x];
coordinate_y_1=[0,Lambda_y];

r_index_1=[1.5, 1.0];