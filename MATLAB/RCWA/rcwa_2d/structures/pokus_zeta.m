number_of_layers=1;

mu=1;
lambda=1;
n_1=1;
n_3=1.5;
theta0=30;
phi0=30;

l=1; % layer 1
thickness(l)=lambda;

Lambda_x=2*lambda;
Lambda_y=Lambda_x;


coordinate_x_1=[0,lambda/2,Lambda_x-lambda/2, Lambda_x];
coordinate_y_1=[0,lambda/2,Lambda_y-lambda/2, Lambda_y];

n_ag=conj(sqrt(5-1*1i)*1i);

r_index_1=[n_ag, n_ag, n_ag;
    n_ag, 1, n_ag;
    n_ag, n_ag, n_ag];