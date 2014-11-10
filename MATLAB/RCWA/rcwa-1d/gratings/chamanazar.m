phi=-30; %90 --- K vector is parallel to the axis x, -> 0 --- K vector is parallel to the axis z
Lambda=1;  % (um)

% check or change a profile => n(x,z) of the refractive index/epsilon distribution
%epsilon_0=2.25;
%d_epsilon=0.64;
n_0=1.5;
dn_1=0.021;
dn_2=0.000;
thickness_total=loop;  % (um)

save_grating=1; % 1-yes, 2-no
save_grating_as='modulated';
plot_grating=2; % 1-yes, 2-no

discretization_x=100; % number of discretization points along x-axis 
discretization_z=100; % number of discretization points along z-axis


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
        epsilon=2.25*(1+0.33*cos(K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad))));
        n(x,z)=sqrt(epsilon);
        %n(x,z)=n_0+dn_1*cos(K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad)))+dn_2*cos(2*K*(new_axis_x(x)*sin(phi_rad)+new_axis_z(z)*cos(phi_rad)));
    end
end

min_n=min(n);
max_n=max(n);

%surf(axis_z(1:size_axis_z(2)-1),axis_x(1:size_axis_x(2)-1),n);
%shading flat;
%axis([0 axis_z(size_axis_z(2)-1) 0 axis_x(size_axis_z(2)-1)]);
%view([-90 -90]);


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
    fidd = fopen([save_grating_as '.m'], 'wt');
    if grating_relative_dimensions==2
        fprintf(fidd, '\n%%absolute size%0f ');
    else
        fprintf(fidd, '\n%%relative size%0f ');
    end
    fprintf(fidd, '\n%%pocet_vrstev=%0f ',size_axis_z-1);
    fprintf(fidd, '\n%%Lambda_x=%0f ',Lambda_x);
    fprintf(fidd, '\n%%Lambda=%0f ',Lambda);
    fprintf(fidd, '\n%%thickness_total=%0f ',thickness_total);
    fprintf(fidd, '\ndata=[ ');
    if (size_axis_z(2)-1)==1
        fprintf(fidd, '\n');
        for j=1:1:2*size_axis_x(2)
            if j==1
                fprintf(fidd, '%0f', thickness_total/discretization_z);
            else
                if imag(result_matrix(1,j))>0 || imag(result_matrix(1,j))<0
                    fprintf(fidd, ', %0f', real(result_matrix(1,j)));                        
                    if imag(result_matrix(1,j))>0                    
                        fprintf(fidd, '+%0f', imag(result_matrix(1,j)));                        
                        fprintf(fidd, '*1i');                        
                    else
                        fprintf(fidd, '%0f', imag(result_matrix(1,j)));                        
                        fprintf(fidd, '*1i');                        
                    end
                else
                    fprintf(fidd, ', %0f', real(result_matrix(1,j)));                    
                end
            end
        end
    else
        for i=1:1:(size_axis_z(2)-1)
            fprintf(fidd, '\n');
            for j=1:1:2*size_axis_x(2)
                if j==1
                    fprintf(fidd, '%0f', thickness_total/discretization_z);
                else
                    if imag(result_matrix(i,j))>0 || imag(result_matrix(i,j))<0                        
                        fprintf(fidd, ', %0f', real(result_matrix(i,j)));
                        if imag(result_matrix(i,j))>0
                            fprintf(fidd, '+%0f', imag(result_matrix(i,j)));
                            fprintf(fidd, '*1i');
                        else
                            fprintf(fidd, '%0f', imag(result_matrix(i,j)));
                            fprintf(fidd, '*1i');
                        end
                    else
                        fprintf(fidd, ', %0f', real(result_matrix(i,j)));
                    end
                end
            end
            if i<(size_axis_z(2)-1)
                fprintf(fidd, ';');
            end
        end
    end
    fprintf(fidd, '\n ];');
    fclose(fidd);
end
movefile('modulated.m','gratings/modulated.m')