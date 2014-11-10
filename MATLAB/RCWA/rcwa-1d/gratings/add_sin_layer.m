% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function data_layer=add_sin_layer(Lambda,Lambda_object,layer_thickness_,ref_index_in,ref_index_out,n_layers)

%{
clear all


Lambda=1.5;
Lambda_object=1.5/3;
ref_index_in=3;
ref_index_out=1;
n_layers=21;
layer_thickness_=4;
%}


plot_grating=2;

number_of_repetition=Lambda/Lambda_object;

if rem(number_of_repetition,1)==0
else
    error('Lambda/Lambda_object is not integer')
end

d=layer_thickness_;
sirka=layer_thickness_*1E-6;

rozdeleni=zeros(1,n_layers+1);
new_z=zeros(1,n_layers);
pomocna=zeros(1,n_layers);
duty_cycle=zeros(1,n_layers);
shift=zeros(1,n_layers);

rozdeleni=linspace(-d/2,d/2,n_layers+1);%rozdělím zadaný interval
for i=1:1:n_layers %udělám průměr dvou sousedních y
    new_z(i)=(rozdeleni(i)+rozdeleni(i+1))/2;
end
pomocna=Lambda_object/(2*pi)*acos(2*new_z/d);%inverzní fuknkce

duty_cycle=fliplr(2*pomocna/Lambda_object);
layer_thickness=sirka/n_layers*ones(1,n_layers);%případně zakomentovat
shift=.5*ones(1,n_layers);

switch plot_grating
    case 1
        figure;
        x=-Lambda_object/2:Lambda_object/100:Lambda_object/2;
        y1=(d/2)*cos(2*pi/Lambda_object*x);
        hold on
        plot(x,y1,'r','LineWidth',2);
        line([-Lambda_object/2,Lambda_object/2],[-d/2,-d/2],'Color','r','LineWidth',2);
        xlabel('\Lambda_object [\mu{}m]','FontSize',18);
        ylabel('d [\mu{}m]','FontSize',18);
        for k=1:1:n_layers
            z=-pomocna(k):Lambda_object/100:pomocna(k);
            plot(z,0*z+rozdeleni(k))
            plot(z,0*z+rozdeleni(k)+d/n_layers)
            line([pomocna(k);pomocna(k)],[new_z(k)-d/(2*n_layers);new_z(k)+d/(2*n_layers)])
            line([-pomocna(k);-pomocna(k)],[new_z(k)-d/(2*n_layers);new_z(k)+d/(2*n_layers)])
        end
    case 2
end

for i=1:1:n_layers
    boundaries(i,1:2)=[Lambda_object/2-Lambda_object*duty_cycle(i)/2,Lambda_object/2+Lambda_object*duty_cycle(i)/2];
end


data_layer=zeros(n_layers,8);

% add thickness
data_layer(:,1)=layer_thickness*1e6;

% add origin of coordinate
data_layer(:,2)=0;

% add outer refractive index
data_layer(:,3)=ref_index_out;

% add first boundary
data_layer(:,4)=boundaries(:,1);

% add inner refractive index
data_layer(:,5)=ref_index_in;

% add first boundary
data_layer(:,6)=boundaries(:,2);

% add outer refractive index
data_layer(:,7)=ref_index_out;

% add Lambda_ coordinate
data_layer(:,8)=Lambda_object;



if number_of_repetition>1
    add_data_layer=data_layer(:,3:8);
    for i=1:1:number_of_repetition
        add_data_layer_=add_data_layer;
        add_data_layer(:,2)=add_data_layer_(:,2)+Lambda_object;
        add_data_layer(:,4)=add_data_layer_(:,4)+Lambda_object;
        add_data_layer(:,6)=add_data_layer_(:,6)+Lambda_object;
        data_layer=[data_layer,add_data_layer_];
    end
end