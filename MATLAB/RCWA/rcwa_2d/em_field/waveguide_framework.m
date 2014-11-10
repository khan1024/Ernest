
if calculate_field==1
    r_index_matrix=eval(['r_index_' num2str(xy_slice_section_to_view)]); % nahraje aktualni matici eps
    coord_x=eval(['coordinate_x_' num2str(xy_slice_section_to_view)]);
    coord_y=eval(['coordinate_y_' num2str(xy_slice_section_to_view)]);
end


r_index_matrix=rot90(r_index_matrix,3);
coordinate_y_l=fliplr(coord_x);
coordinate_x_l=fliplr(coord_y);

size_eps=size(r_index_matrix);

n_of_rows=size_eps(1);
n_of_columns=size_eps(2);

% yes boundary =1
% no boundary =0

matrix_geometry=zeros(n_of_rows,n_of_columns,4);

% left boundary
matrix_geometry(:,1,1)=2;

% top boundary
matrix_geometry(1,:,2)=2;

% right boundary
matrix_geometry(:,end,3)=2;

% bottom edge
matrix_geometry(end,:,4)=2;

Z=[100,100];

for i=1:1:n_of_rows
    for j=1:1:n_of_columns
        
        % left
        if j==1 % již vyřešeno
        else
            if r_index_matrix(i,j)==r_index_matrix(i,j-1)
                matrix_geometry(i,j,1)=2;
            else
                matrix_geometry(i,j,1)=1;
            end
            
        end
        
        % top
        if i==1 % již vyřešeno
        else
            if r_index_matrix(i,j)==r_index_matrix(i-1,j)
                matrix_geometry(i,j,2)=2;
            else
                matrix_geometry(i,j,2)=1;
            end
            
        end
        
        % right
        if j==n_of_columns % již vyřešeno
        else
            if r_index_matrix(i,j)==r_index_matrix(i,j+1)
                matrix_geometry(i,j,3)=2;
            else
                matrix_geometry(i,j,3)=1;
            end
            
        end
        
        % bottom
        if i==n_of_rows
        else
            if r_index_matrix(i,j)==r_index_matrix(i+1,j)
                matrix_geometry(i,j,4)=2;
            else
                matrix_geometry(i,j,4)=1;
            end
            
        end
        
        
    end
end

%figure;

for i=1:1:n_of_rows
    for j=1:1:n_of_columns
        % left
        
        if matrix_geometry(i,j,1)==2  %nic
        else
            x_1=coordinate_x_l(j);
            y_1=coordinate_y_l(i);
            
            x_2=coordinate_x_l(j);
            y_2=coordinate_y_l(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[0 0 0]);
            
        end
        
        % top
        
        if matrix_geometry(i,j,2)==2 %nic
        else
            
            x_1=coordinate_x_l(j);
            y_1=coordinate_y_l(i);
            
            x_2=coordinate_x_l(j+1);
            y_2=coordinate_y_l(i);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[0 0 0]);
            
        end
        
        
        % right
        
        if matrix_geometry(i,j,3)==2 %nic
        else
            x_1=coordinate_x_l(j+1);
            y_1=coordinate_y_l(i);
            
            x_2=coordinate_x_l(j+1);
            y_2=coordinate_y_l(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[0 0 0]);
            
        end
        
        
        % bottom
        
        if matrix_geometry(i,j,4)==2 %nic
        else
            
            x_1=coordinate_x_l(j+1);
            y_1=coordinate_y_l(i+1);
            
            x_2=coordinate_x_l(j);
            y_2=coordinate_y_l(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[0 0 0]);
            
            
        end
        
    end
end