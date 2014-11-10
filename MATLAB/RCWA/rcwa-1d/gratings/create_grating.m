
% ? number of layers

number_of_layers=0;
next_layer=1;

while next_layer==1
    number_of_layers=number_of_layers+1;
    if exist(['layer_' num2str(number_of_layers)],'var')==1
        next_layer=1;
    else
        next_layer=0;
        number_of_layers=number_of_layers-1;
    end
end

% joint layers

data=[];
for i=1:1:number_of_layers
    
    add_layer_data=eval(['layer_' num2str(i)]);
    
    %Dimensions of matrices must be same -> add NaN
    % 1) size_data(2)<size_add_layer_data(2) -> add NaN to data
    % 2) size_data(2)>size_add_layer_data(2) -> add NaN to add_layer_data  
    
    data_size=size(data);
    add_layer_data_size=size(add_layer_data);
    
    if data_size(2)<add_layer_data_size(2)
        new_data=fill_NaN(data_size(1),add_layer_data_size(2));
        new_data(:,1:data_size(2))=data;
        data=new_data;
    elseif data_size(2)>add_layer_data_size(2)
        new_add_layer_data=fill_NaN(add_layer_data_size(1),data_size(2));
        new_add_layer_data(:,1:add_layer_data_size(2))=add_layer_data;
        add_layer_data=new_add_layer_data; 
    else
        %dimensions are same -> OK
    end
    
    data=[data;add_layer_data];
end

if plot_grating==1
    figure;
    plot_grating_file;
end
