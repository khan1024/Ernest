clear all

phi=13; %90 --- K vector is parallel to the axis x, -> 0 --- K vector is parallel to the axis z
Lambda=0.35228;  % (um)

% check or change a profile => n(x,z) of the refractive index/epsilon distribution
%epsilon_0=2.25;
%d_epsilon=0.64;
n_0=1.5;
dn_1=0.07;
dn_2=0.000;
thickness_total=10;  % (um)

save_grating=1; % 1-yes, 2-no
save_grating_as='modulated2';
plot_grating=2; % 1-yes, 2-no

discretization_x=100; % number of discretization points along x-axis 
discretization_z=300; % number of discretization points along z-axis


grating_relative_dimensions=2; % 1- yes relative_dimensions -> Lambda=1, thickness_total=1; -> in the file main.m -> grating_relative_dimensions=1 + set new Lambda and thickness_total 
                               % 2- no -> absolute dimensions -> in the file main.m -> grating_relative_dimensions=2

%--------------------------------------------------------------
if phi==0
    error('ilegal value')
end
phi_rad=-(phi)*pi/180;

if grating_relative_dimensions==1
    thickness_total=1;
    Lambda=1;
end

K=2*pi/Lambda;
KX=abs(K*sin(phi_rad));
KZ=K*cos(phi_rad);

Lambda_x=2*pi/KX;
axis_x=0:Lambda_x/discretization_x:Lambda_x;
axis_z=0:thickness_total/discretization_z:thickness_total;

size_axis_x=size(axis_x);
size_axis_z=size(axis_z);


for i=1:1:size_axis_x(2)-1
    new_axis_x(i)=(axis_x(i)+axis_x(i+1))/2;
end

for i=1:1:size_axis_z(2)-1
    new_axis_z(i)=(axis_z(i)+axis_z(i+1))/2;
end

for x=1:1:size_axis_x(2)-1
    for z=1:1:size_axis_z(2)-1;
        %epsilon=2.25*(1+0.33*cos(K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad))));
        %n(x,z)=sqrt(epsilon);
        n(x,z)=n_0+dn_1*cos(K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad)))+dn_2*cos(2*K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad)));
    end
end

min_n=min(n);
max_n=max(n);

figure
surf(axis_z(1:size_axis_z(2)-1),axis_x(1:size_axis_x(2)-1),n);
shading flat;
axis([0 axis_z(size_axis_z(2)-1) 0 axis_x(size_axis_x(2)-1)]);
view([-90 -90]);


refractive_index_matrix=rot90(n,-1);

result_matrix_1=zeros(size_axis_z(2)-1,2*size_axis_x(2)-1);

for j=1:1:size_axis_x(2)
    for i=1:1:size_axis_z(2)-1
        result_matrix_1(i,2*j-1)=axis_x(j);
    end
end

for j=1:1:size_axis_x(2)-1
    for i=1:1:size_axis_z(2)-1
        result_matrix_1(i,2*j)=refractive_index_matrix(i,j);
    end
end

result_matrix=[axis_z(2:size_axis_z(2))' result_matrix_1];

if save_grating==1
    fid = fopen([save_grating_as '.m'], 'wt');
    if grating_relative_dimensions==2
        fprintf(fid, '\n%%absolute size%0f ');
    else
        fprintf(fid, '\n%%relative size%0f ');
    end
    fprintf(fid, '\n%%pocet_vrstev=%0f ',size_axis_z-1);
    fprintf(fid, '\n%%Lambda_x=%0f ',Lambda_x);
    fprintf(fid, '\n%%Lambda=%0f ',Lambda);
    fprintf(fid, '\n%%thickness_total=%0f ',thickness_total);
    fprintf(fid, '\ndata=[ ');
    if (size_axis_z(2)-1)==1
        fprintf(fid, '\n');
        for j=1:1:2*size_axis_x(2)
            if j==1
                fprintf(fid, '%0f', thickness_total/discretization_z);
            else
                if imag(result_matrix(1,j))>0 || imag(result_matrix(1,j))<0
                    fprintf(fid, ', %0f', real(result_matrix(1,j)));                        
                    if imag(result_matrix(1,j))>0                    
                        fprintf(fid, '+%0f', imag(result_matrix(1,j)));                        
                        fprintf(fid, '*1i');                        
                    else
                        fprintf(fid, '%0f', imag(result_matrix(1,j)));                        
                        fprintf(fid, '*1i');                        
                    end
                else
                    fprintf(fid, ', %0f', real(result_matrix(1,j)));                    
                end
            end
        end
    else
        for i=1:1:(size_axis_z(2)-1)
            fprintf(fid, '\n');
            for j=1:1:2*size_axis_x(2)
                if j==1
                    fprintf(fid, '%0f', thickness_total/discretization_z);
                else
                    if imag(result_matrix(i,j))>0 || imag(result_matrix(i,j))<0                        
                        fprintf(fid, ', %0f', real(result_matrix(i,j)));
                        if imag(result_matrix(i,j))>0
                            fprintf(fid, '+%0f', imag(result_matrix(i,j)));
                            fprintf(fid, '*1i');
                        else
                            fprintf(fid, '%0f', imag(result_matrix(i,j)));
                            fprintf(fid, '*1i');
                        end
                    else
                        fprintf(fid, ', %0f', real(result_matrix(i,j)));
                    end
                end
            end
            if i<(size_axis_z(2)-1)
                fprintf(fid, ';');
            end
        end
    end
    fprintf(fid, '\n ];');
    fclose(fid);
end
break
if plot_grating==1;
    figure
    for i=1:1:(size_axis_z(2)-1)
        result_matrix(i,1)=thickness_total/discretization_z;
    end
    
    data=result_matrix;
    size_result_matrix_m=size(data);
    result_matrix_m=data;

    y_max=0;
    for i=1:1:size_result_matrix_m(1)    
        y_max=y_max+result_matrix_m(i,1);
    end

    y_predesla=0;
    for i=1:1:size_result_matrix_m(1)    
        for j=2:2:size_result_matrix_m(2)-1        
            if isnan(result_matrix_m(i,j+2))==1            
        
            else
                x_pocatek=result_matrix_m(i,j);            
                x_konec=result_matrix_m(i,j+2);            
                x_width=x_konec-x_pocatek;            
                y_pocatek=y_max-result_matrix_m(i,1)-y_predesla;            
                y_height=result_matrix_m(i,1);            
                barva=[real(result_matrix_m(i,j+1)/2),0,0];            
                if x_width>0                
                    rectangle('Position',[x_pocatek,y_pocatek,x_width,y_height],'FaceColor',barva)            
                end
            end
        end
        y_predesla=y_predesla+result_matrix_m(i,1);
    end

    xlabel('\Lambda [\mu{}m]','FontSize',17);
    ylabel('d [\mu{}m]','FontSize',17);
end