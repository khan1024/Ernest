% coordinate_z is easy

temp_coordinate_z=layer_thickness*1e6;

coordinate_z=[0];

for i=1:1:length(temp_coordinate_z)
    coordinate_z=[coordinate_z temp_coordinate_z(i)+coordinate_z(end)];
end
%

% pridat region 1 a region 3
coordinate_z=[-size_region_1, coordinate_z, coordinate_z(end)+size_region_3];

temp_coordinate_x=[0, Lambda];

data=zeros(length(duty_cycle),8);
data(:,1)=temp_coordinate_z;
grating_0_inversion=2;

for i=1:1:length(duty_cycle);
    coord_1=shift(i)*Lambda-duty_cycle(i)*Lambda/2;
    coord_2=shift(i)*Lambda+duty_cycle(i)*Lambda/2;
    
    if coord_1==0
        coord_1=coord_1+1e-10;
    end
    
    if coord_2==Lambda
        coord_2=coord_2-1e-10;
    end
    
    if grating==0 % lze muze shift+duty_cycle přesahovat periodu -> zaporné hodnoty
        
        if coord_1<0
            coord_1=Lambda+coord_1;
            grating_0_inversion=1;
        end
        
        if coord_2>Lambda
            coord_2=coord_2-Lambda;
            grating_0_inversion=1;
        end
        
        if grating_0_inversion==1
            data(i,3)=nr;        
            data(i,4)=coord_2;        
            data(i,5)=ng;        
            data(i,6)=coord_1;        
            data(i,7)=nr;
        else
            data(i,3)=ng;        
            data(i,4)=coord_1;        
            data(i,5)=nr;        
            data(i,6)=coord_2;        
            data(i,7)=ng;
        end
        data(i,8)=Lambda;
    else
        data(i,3)=ng;
        data(i,4)=coord_1;
        data(i,5)=nr;
        data(i,6)=coord_2;
        data(i,7)=ng;
        data(i,8)=Lambda;
    end
    
    temp_coordinate_x=[temp_coordinate_x,coord_1];
    temp_coordinate_x=[temp_coordinate_x,coord_2];
end


sort_temp_coordinate_x=unique(sort(temp_coordinate_x));

coordinate_x=sort_temp_coordinate_x;

ref_index_matrix=zeros(length(coordinate_z)-1,length(coordinate_x)-1);

% prvni radek = region_1
ref_index_matrix(1,:)=n1;

% prvni radek = region_1
ref_index_matrix(end,:)=n3;

% prvni ref index je jasný
for i=1:1:length(duty_cycle)
    current_layer_coordinate_x=data(i,2:2:end);
    
    if grating_0_inversion==1
        ref_index_matrix(i+1,1)=nr;
    else
        ref_index_matrix(i+1,1)=ng;
    end
    
    for j=2:1:length(coordinate_x)-1
        ind=find(current_layer_coordinate_x>coordinate_x(j));
        if isempty(ind)==1
        else
            %[ind(1) data(i,2*ind(1)-1)]
            ref_index_matrix(i+1,j)=data(i,2*ind(1)-1);
        end 
    end
end