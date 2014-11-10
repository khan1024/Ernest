% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%function [data,number_of_layers,Lambda_um,thickness_total_a]=kolik_ma_vrstvicek(input_grating_file)

% zatím nevím, jak jinak nahrát *.m soubor, dočasně vytvořím kopii zadaného souboru,
% ten nahraju a následně ho vymažu
cd('gratings');
copyfile(input_grating_file,'temporary.m')
temporary
delete('temporary.m')
cd('..');

size_data=size(data);
number_of_layers=size_data(1);

%-------------------------------------------------------------------------
% if grating is described by absolute size -> we need to know 1) Lambda, 2) total thickness

% 1) Lambda
size_data_x=size_data(2);
index=1:1:size_data_x; % temporary field

for i=3:2:size_data_x
    if isnan(data(1,index(i)))==1
    else
        if data(1,2)<0 || data(1,2)>0
            Lambda_um=data(1,index(i+1))-data(1,index(2));
        else
            Lambda_um=data(1,index(i+1));
        end
    end
end

% 2) total thickness
thickness_total_a=sum(data(:,1));

clear size_data size_data_x index
