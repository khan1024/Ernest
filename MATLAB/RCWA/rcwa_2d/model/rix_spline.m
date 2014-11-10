function [n_model]=rix_spline(lambda,nahraj_material)

filename=[nahraj_material];
data_s=dlmread(filename);

velikost=size(data_s);

% kontrola zadání
if lambda>=data_s(1,1) && lambda<=data_s(velikost(1),1)
else
    error('lambda is out of range')
end

data_lambda=data_s(:,1);
real_rix=data_s(:,2);
imag_rix=data_s(:,3);

n_r=spline(data_lambda,real_rix,lambda);
n_i=spline(data_lambda,imag_rix,lambda);

n_model=n_r-n_i*1i;