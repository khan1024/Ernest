% je potřeba vytvořit zcela novou matici, která bude obsahovat všechny
% souřadnice

% coordinate_z is easy
temp_coordinate_z=data(:,1);
coordinate_z=[0];

for i=1:1:length(temp_coordinate_z)
    coordinate_z=[coordinate_z temp_coordinate_z(i)+coordinate_z(end)];
end
%

% pridat region 1 a region 3
coordinate_z=[-size_region_1, coordinate_z, coordinate_z(end)+size_region_3];

size_data=size(data);

% nashromáždit všechny souřadnice x, seřadit je a opakující se souřadnice
% vyřadit

temp_coordinate_x=[];

for i=1:1:size_data(1);
    for j=2:2:size_data(2)
        temp_coordinate_x=[temp_coordinate_x,data(i,j)];
    end
end

sort_temp_coordinate_x=unique(sort(temp_coordinate_x));

% ještě odstranit NaN
coordinate_x=[];

for i=1:1:length(sort_temp_coordinate_x)
    if isnan(sort_temp_coordinate_x(i))==1
    else
        coordinate_x=[coordinate_x sort_temp_coordinate_x(i)];
    end
end

% vytvorit vyslednou matici

ref_index_matrix=zeros(length(coordinate_z)-1,length(coordinate_x)-1);

% prvni radek = region_1
ref_index_matrix(1,:)=n1;

% prvni radek = region_1
ref_index_matrix(end,:)=n3;


% prvni ref index je jasný
for i=1:1:size_data(1)
    current_layer_coordinate_x=data(i,2:2:end);
    ref_index_matrix(i+1,1)=data(i,3);
    
    for j=2:1:length(coordinate_x)-1
        ind=find(current_layer_coordinate_x>coordinate_x(j));
        if isempty(ind)==1
        else
            %[ind(1) data(i,2*ind(1)-1)]
            ref_index_matrix(i+1,j)=data(i,2*ind(1)-1);
        end 
    end
end