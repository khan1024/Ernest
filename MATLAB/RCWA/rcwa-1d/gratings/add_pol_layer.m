function data_layer=add_pol_layer(Lambda,x,y,ref_index_in,ref_index_out,n_layers)

%{
clear all

x=[0.1 0.4 0.6 0.9 0.6 0.6 0.4 0.1];
y=[0 0 0 0 1 0 0 1];

ref_index_out=1;
ref_index_in=3.5;

n_layers=10;

fill(x,y,'r')

Lambda=1;

%}

x=[x x(1)];
y=[y y(1)];

Lambda_object=max(x)-min(x);
thickness_object=max(y)-min(y);

d=thickness_object;
sirka=thickness_object*1E-6;

rozdeleni=zeros(1,n_layers+1);
new_z=zeros(1,n_layers);
pomocna=zeros(1,n_layers);
duty_cycle=zeros(1,n_layers);
shift=zeros(1,n_layers);

layer_thickness=sirka/n_layers*ones(1,n_layers);

rozdeleni=linspace(0,d,n_layers+1);%rozdělím zadaný interval
for i=1:1:n_layers %udělám průměr dvou sousedních y
    new_z(i)=(rozdeleni(i)+rozdeleni(i+1))/2;
end

% test

boundary_test=polyxpoly([0,Lambda],[thickness_object/2,thickness_object/2],x,y);
length_boundary_test=length(boundary_test);

boundaries=zeros(n_layers,length_boundary_test);

for i=1:1:n_layers
    % new_z(i)~y
    if isempty(find(new_z(i)==y, 1))==0
        new_z(i)=new_z(i)+1e-10;
    end
    [boundaries(i,:),b]=polyxpoly([0,Lambda],[new_z(i),new_z(i)],x,y);
end

boundaries=rot90(boundaries,2);

data_layer=zeros(n_layers,4+2*length_boundary_test);

% add thickness
data_layer(:,1)=layer_thickness*1e6;

% add origin of coordinate
data_layer(:,2)=0;

index_data_layer=3;
for i=1:1:length_boundary_test/2
% add outer refractive index
data_layer(:,4*i-1)=ref_index_out;

% add first boundary
data_layer(:,4*i)=boundaries(:,2*i-1);

% add inner refractive index
data_layer(:,4*i+1)=ref_index_in;

% add first boundary
data_layer(:,4*i+2)=boundaries(:,2*i);
end

% add outer refractive index
data_layer(:,end-1)=ref_index_out;

data_layer(:,end)=Lambda;


