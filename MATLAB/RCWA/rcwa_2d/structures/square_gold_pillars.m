number_of_layers=1;

mu=1;
%lambda=0.54;
n_l=1.33;
n_1=n_l;
n_3=sqrt(1+1.03961212*lambda^2/(lambda^2-0.00600069867) + 0.231792344*lambda^2/(lambda^2-0.0200179144) + 1.01046945*lambda^2/(lambda^2-103.560653));

theta0=20;

l=1; % layer 1
thickness(l)=0.063;

Lambda_x=0.507;
Lambda_y=Lambda_x;

w_x=0.146;
w_y=0.129;

coordinate_x_1=[0,Lambda_x/2-w_x/2,Lambda_x/2+w_x/2,Lambda_x];
coordinate_y_1=[0,Lambda_y/2-w_y/2,Lambda_y/2+w_y/2,Lambda_y];

n_g=gold;

r_index_1=[n_l, n_l, n_l; n_l, n_g, n_l; n_l, n_l, n_l];