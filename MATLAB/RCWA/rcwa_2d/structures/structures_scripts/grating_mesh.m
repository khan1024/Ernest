% yes boundary =1
% no boundary =0

matrix_geometry=zeros(n_of_rows,n_of_columns,4);

% left boundary
matrix_geometry(:,1,1)=1;

% top boundary
matrix_geometry(1,:,2)=1;

% right boundary
matrix_geometry(:,end,3)=1;

% bottom edge
matrix_geometry(end,:,4)=1;

for i=1:1:n_of_rows
    for j=1:1:n_of_columns
        
        % left
        if j==1 % již vyřešeno
        else
            if mesh_n_mat(i,j)==mesh_n_mat(i,j-1)
                matrix_geometry(i,j,1)=2;
            else
                matrix_geometry(i,j,1)=1;
            end
            
        end
        
        % top
        if i==1 % již vyřešeno
        else
            if mesh_n_mat(i,j)==mesh_n_mat(i-1,j)
                matrix_geometry(i,j,2)=2;
            else
                matrix_geometry(i,j,2)=1;
            end
            
        end
        
        % right
        if j==n_of_columns % již vyřešeno
        else
            if mesh_n_mat(i,j)==mesh_n_mat(i,j+1)
                matrix_geometry(i,j,3)=2;
            else
                matrix_geometry(i,j,3)=1;
            end
            
        end
        
        % bottom
        if i==n_of_rows
        else
            if mesh_n_mat(i,j)==mesh_n_mat(i+1,j)
                matrix_geometry(i,j,4)=2;
            else
                matrix_geometry(i,j,4)=1;
            end
            
        end
        
        
    end
end

figure;

for i=1:1:n_of_rows
    for j=1:1:n_of_columns
        % left
        
        if matrix_geometry(i,j,1)==2 %nic
        else
            x_1=p_x(j);
            y_1=p_y(i);
            
            x_2=p_x(j);
            y_2=p_y(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y);
            
        end
        
        % top
        
        if matrix_geometry(i,j,2)==2 %nic
        else
            
            x_1=p_x(j);
            y_1=p_y(i);
            
            x_2=p_x(j+1);
            y_2=p_y(i);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y);
            
        end
        
        
        % right
        
        if matrix_geometry(i,j,3)==2 %nic
        else
            x_1=p_x(j+1);
            y_1=p_y(i);
            
            x_2=p_x(j+1);
            y_2=p_y(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y);
            
        end
        
        
        % bottom
        
        if matrix_geometry(i,j,4)==2 %nic
        else
            
            x_1=p_x(j+1);
            y_1=p_y(i+1);
            
            x_2=p_x(j);
            y_2=p_y(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y);
            
            
        end
        
    end
end

view([0 -90])