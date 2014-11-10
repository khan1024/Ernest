clear all

% this scipt generates refractive index distribution of the volume grating as a input file for main.m 
% check or change a profile => profile_n of the refractive index/epsilon distribution

n_average=1.5;
n_amplitude_1=0.01;
%n_amplitude_1=1E-3;
Lambda=1.0; % (um)
thickness_total=50; % (um)
discretization_x=400; % number of discretization points along x-axis

save_as='volume_grating_example';

grating_relative_dimensions=2; % 1- yes relative_dimensions -> Lambda=1, thickness_total=1; -> in the file main.m -> grating_relative_dimensions=1 + set new Lambda and thickness_total 
                               % 2- no -> absolute dimensions -> in the file main.m -> grating_relative_dimensions=2


x=0:Lambda/discretization_x:Lambda;
size_x=size(x);
number_of_elements=size_x(2);

for i=1:1:number_of_elements-1
    new_x(i)=(x(i)+x(i+1))/2;
end

% !!!!!
profile_n=n_average+n_amplitude_1*sin(2*pi/Lambda*new_x);%+n_amplitude_2*sin(2*2*pi/Lambda*new_x);

% sometimes you need somethink like:
% profile_epsilon=eps_average+eps_amplitude_1*sin(2*pi/Lambda*new_x);
% profile_n=sqrt(profile_epsilon);
%-----------------------------------------------------------------------


if grating_relative_dimensions==1
    for i=1:1:number_of_elements-1    
        data(2*i)=x(i)/Lambda;    
        data(2*i+1)=profile_n(i);
    end

    size_data=size(data);
    data(1)=1;
    data(size_data(2)+1)=1;
else
    for i=1:1:number_of_elements-1    
        data(2*i)=x(i);    
        data(2*i+1)=profile_n(i);
    end

    size_data=size(data);
    data(1)=thickness_total;
    data(size_data(2)+1)=Lambda;
end



size_data=size(data);

fid = fopen([save_as '.m'], 'wt');
fprintf(fid, '\n%%Lambda=%0f ',Lambda);
fprintf(fid, '\ndata=[ ');
for i=1:1:size_data(2)
    if i==1
        fprintf(fid, '\n%0f', data(i));
    else
        fprintf(fid, ', %0f', data(i));
    end

end

fprintf(fid, '\n ];');
fclose(fid);



