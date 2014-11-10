%clear all
%gold=1.25;
%lambda=0.8;

number_of_layers=2;

mu=1;
%lambda=0.54;
n_l=1.33;
n_1=n_l;
n_3=sqrt(1+1.03961212*lambda^2/(lambda^2-0.00600069867) + 0.231792344*lambda^2/(lambda^2-0.0200179144) + 1.01046945*lambda^2/(lambda^2-103.560653));

theta0=0;

d_ch=0.01;
n_dch=1.5;

Lambda_x=0.507;
Lambda_y=Lambda_x;

w_x=0.146;
w_y=0.129;
n_g=gold;

l=1; % layer 1
thickness(l)=d_ch;

coordinate_x_1=[0,Lambda_x/2-w_x/2-d_ch,Lambda_x/2-w_x/2,Lambda_x/2+w_x/2,Lambda_x/2+w_x/2+d_ch,Lambda_x];
coordinate_y_1=[0,Lambda_y/2-w_y/2-d_ch,Lambda_y/2-w_y/2,Lambda_y/2+w_y/2,Lambda_y/2+w_y/2+d_ch,Lambda_y];

r_index_1=[n_l, n_l, n_l, n_l, n_l; n_l, n_dch, n_dch, n_dch, n_l; n_l, n_dch, n_dch, n_dch, n_l; n_l, n_dch, n_dch, n_dch, n_l; n_l, n_l, n_l, n_l, n_l];

l=2; % layer 1
thickness(l)=0.063;

coordinate_x_2=[0,Lambda_x/2-w_x/2-d_ch,Lambda_x/2-w_x/2,Lambda_x/2+w_x/2,Lambda_x/2+w_x/2+d_ch,Lambda_x];
coordinate_y_2=[0,Lambda_y/2-w_y/2-d_ch,Lambda_y/2-w_y/2,Lambda_y/2+w_y/2,Lambda_y/2+w_y/2+d_ch,Lambda_y];

r_index_2=[n_l, n_l, n_l, n_l, n_l; n_l, n_dch, n_dch, n_dch, n_l;n_l, n_dch, n_g, n_dch, n_l; n_l, n_dch, n_dch, n_dch, n_l; n_l, n_l, n_l, n_l, n_l];