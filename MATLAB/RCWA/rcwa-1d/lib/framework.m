
if grating==10
    framework_grating_10;
else
    framework_grating_0_9;
end

if grating==10
    epsilon_1=flipud(rot90(ref_index_matrix,3));
else
    epsilon_1=rot90(ref_index_matrix,3);
end

coordinate_y_1=fliplr(coordinate_x);
coordinate_x_1=fliplr(coordinate_z);

%coordinate_y_1=[2,1.75,1.25,1,0.75,0.25,0];
%coordinate_y_1=[(coordinate_y_1+Lambda) coordinate_y_1];

epsilon_1_0=epsilon_1;
coordinate_y_1_0=coordinate_y_1;
if number_of_plotted_period>=2
    for i=2:1:number_of_plotted_period
        epsilon_1=[epsilon_1_0;epsilon_1];
        coordinate_y_1=[(coordinate_y_1_0(1:end-1)+(i-1)*Lambda) coordinate_y_1];
    end
end

size_eps=size(epsilon_1);

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

Z=[2000,2000];

for i=1:1:n_of_rows
    for j=1:1:n_of_columns
        
        % left
        if j==1 % již vyřešeno
        else
            if epsilon_1(i,j)==epsilon_1(i,j-1)
                matrix_geometry(i,j,1)=2;
            else
                matrix_geometry(i,j,1)=1;
            end
            
        end
        
        % top
        if i==1 % již vyřešeno
        else
            if epsilon_1(i,j)==epsilon_1(i-1,j)
                matrix_geometry(i,j,2)=2;
            else
                matrix_geometry(i,j,2)=1;
            end
            
        end
        
        % right
        if j==n_of_columns % již vyřešeno
        else
            if epsilon_1(i,j)==epsilon_1(i,j+1)
                matrix_geometry(i,j,3)=2;
            else
                matrix_geometry(i,j,3)=1;
            end
            
        end
        
        % bottom
        if i==n_of_rows
        else
            if epsilon_1(i,j)==epsilon_1(i+1,j)
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
            x_1=coordinate_x_1(j);
            y_1=coordinate_y_1(i);
            
            x_2=coordinate_x_1(j);
            y_2=coordinate_y_1(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[1 1 1]);
            
        end
        
        % top
        
        if matrix_geometry(i,j,2)==2 %nic
        else
            
            x_1=coordinate_x_1(j);
            y_1=coordinate_y_1(i);
            
            x_2=coordinate_x_1(j+1);
            y_2=coordinate_y_1(i);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[1 1 1]);
            
        end
        
        
        % right
        
        if matrix_geometry(i,j,3)==2 %nic
        else
            x_1=coordinate_x_1(j+1);
            y_1=coordinate_y_1(i);
            
            x_2=coordinate_x_1(j+1);
            y_2=coordinate_y_1(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[1 1 1]);
            
        end
        
        
        % bottom
        
        if matrix_geometry(i,j,4)==2 %nic
        else
            
            x_1=coordinate_x_1(j+1);
            y_1=coordinate_y_1(i+1);
            
            x_2=coordinate_x_1(j);
            y_2=coordinate_y_1(i+1);
            
            X=[x_1,x_2];
            Y=[y_1,y_2];
            
            line(X,Y,Z,'Color',[1 1 1]);
            
            
        end
        
    end
end