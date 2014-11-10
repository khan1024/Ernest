function [n_model]=drude(lambda,p)

lambda_=10000/lambda;

e_1=p(2)^2/(lambda_^2+p(1)^2);
e_2=(p(2)^2*p(1))/(lambda_*(lambda_^2+p(1)^2));

n_model=-j*sqrt(e_1+j*e_2);