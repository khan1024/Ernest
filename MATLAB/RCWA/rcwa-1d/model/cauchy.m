
function [n_model]=cauchy(lambda,A,B,C)
% kontrola velikosti

n_model=A+B/lambda^2+C/lambda^4;