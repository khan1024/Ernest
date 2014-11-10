% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%clear all
%number_of_layers=1;
%thickness_total=1;
%thickness=thickness_total*1E-6;
%l=1;
%M=10;
%addpath('gratings');
%input_grating_file='data_m4.txt';

%filename=[input_grating_file];
%data=dlmread(filename); 

size_data=size(data);
size_data_x=size_data(2);
size_data_y=size_data(1);

%layer_thickness_pom=ones(1,number_of_layers);

layer_thickness=1e-6*data(:,1);

%for i=1:1:number_of_layers
%    layer_thickness(i)=layer_thickness_pom(i)*1E-6*data(i,1);
%end

index=1:1:size_data_x;

%------------------------------------------------------------------------
% zero fourier coefficient
epsG=0;
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        epsG=epsG+data(l,index(i))^2*(data(l,index(i+1))-data(l,index(i-1)))/Lambda_um;
    end
end

% inverse zero fourier coefficient
i_epsG=0;
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        i_epsG=i_epsG+(1/data(l,index(i))^2)*(data(l,index(i+1))-data(l,index(i-1)))/Lambda_um;
    end
end
%------------------------------------------------------------------------
% plus non-zero fourier coefficients
v_p=0*(1:M);
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        v_p=v_p+data(l,index(i))^2*(-exp(-1i*2*pi*(1:M)*data(l,index(i+1))/Lambda_um)./(1i*2*pi*(1:M))+ ...
            exp(-1i*2*pi*(1:M)*data(l,index(i-1))/Lambda_um)./(1i*2*pi*(1:M)));
    end
end

v_0=epsG;

% minus non-zero fourier coefficients
v_m=0*(-M:-1);
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        v_m=v_m+data(l,index(i))^2*(-exp(-1i*2*pi*(-M:-1)*data(l,index(i+1))/Lambda_um)./(1i*2*pi*(-M:-1))+ ...
        exp(-1i*2*pi*(-M:-1)*data(l,index(i-1))/Lambda_um)./(1i*2*pi*(-M:-1)));
    end
end

v=[v_m v_0 v_p]; % fourier coefficients
%------------------------------------------------------------------------
% inverse plus non-zero fourier coefficients
i_vp=0*(1:M);
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        i_vp=i_vp+(1/data(l,index(i))^2)*(-exp(-1i*2*pi*(1:M)*data(l,index(i+1))/Lambda_um)./(1i*2*pi*(1:M))+ ...
        exp(-1i*2*pi*(1:M)*data(l,index(i-1))/Lambda_um)./(1i*2*pi*(1:M)));
    end
end

i_v0=i_epsG;

% inverse minus non-zero fourier coefficients
i_vm=0*(-M:-1);
for i=3:2:size_data_x
    if isnan(data(l,index(i)))==1
    else
        i_vm=i_vm+(1/data(l,index(i))^2)*(-exp(-1i*2*pi*(-M:-1)*data(l,index(i+1))/Lambda_um)./(1i*2*pi*(-M:-1))+ ...
        exp(-1i*2*pi*(-M:-1)*data(l,index(i-1))/Lambda_um)./(1i*2*pi*(-M:-1)));
    end
end

i_v=[i_vm i_v0 i_vp]; % inverse fourier coefficients
