n_1=1;
n_3=0.22-6.71*1i;
lambda=1;
theta0=30.;
Lambda=1;

%data=[0.5,0, n_1, Lambda/4, n_1, 3*Lambda/4, n_1, Lambda;
%    1,0, n_1, Lambda/4, n_3, 3*Lambda/4, n_1, Lambda;
%    0.5,0, n_3, Lambda/4, n_3, 3*Lambda/4, n_3, Lambda;];

data=[0.2,-Lambda/2, n_1, -Lambda/4, n_1, Lambda/4, n_1, Lambda/2;
    1,-Lambda/2, n_1, -Lambda/4, n_3, Lambda/4, n_1, Lambda/2;
    0.7,-Lambda/2, n_3, -Lambda/4, n_3, Lambda/4, n_3, Lambda/2];
