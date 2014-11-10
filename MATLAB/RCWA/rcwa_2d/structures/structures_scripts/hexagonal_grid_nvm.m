%{
clear all

Lambda_x=0.3;
Lambda_y=Lambda_x*sqrt(3);

dx=0.01;
r_cylinder=0.1;
eps_cylinder=1;
eps_layer=3;

show_grating_mesh=1; % 1-yes, 2-no
%}

% nutnost přenastavit Lambda_x a Lambda_y podle aktualni rozliseni
function [N_x_input,N_y_input]=hexagonal_grid_nvm(Lambda_x,Lambda_y,resolution_dx,factorization_method)

dx=Lambda_x/resolution_dx;

Lambda_x_old=Lambda_x;
Lambda_y_old=Lambda_y;

% 1) Lambda_x
fix_Lambda_x=fix(Lambda_x/dx);
Lambda_x=fix_Lambda_x*dx;

% 1) Lambda_y
fix_Lambda_y=fix(Lambda_y/dx);
Lambda_y=fix_Lambda_y*dx;


p_x=0:dx:Lambda_x;
p_y=0:dx:Lambda_y;



p_x_mat=p_x;
p_y_mat=p_y;

switch factorization_method
    case 0 % Li factorization
        N_x_input=[];
        N_y_input=[];
    case 1 % moje normal vector method
        period_x=Lambda_x;
        period_y=Lambda_y;
        step=dx;
        x=-period_x/2+step/2:step:period_x/2-step/2;
        y=-period_y/2+step/2:step:period_y/2-step/2;
        [x,y]=meshgrid(x,y);

        z_mm=-sqrt(x.^2+y.^2);
        z_dl=-sqrt((x+period_x*0.5).^2+(y+period_x*sqrt(3)/2).^2);
        z_dr=-sqrt((x-period_x*0.5).^2+(y+period_x*sqrt(3)/2).^2);

        z_ur=-sqrt((x-period_x*0.5).^2+(y-period_x*sqrt(3)/2).^2);
        z_ul=-sqrt((x+period_x*0.5).^2+(y-period_x*sqrt(3)/2).^2);


        % 1) vytvoreni ctverce o strane 2*a/sqrt(3)

        % 2) oriznutí levych a pravych částí



        size_x_y=size(x);
        for j=1:1:size_x_y(2)
    
            %for i=1:1:size_x_y(1)        
                % z_mm        
                %%{        
                delete=equation_line2(x(:,j),y(:,j),0,period_x*1/sqrt(3),period_x*0.5,period_x*1/sqrt(3)/2,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_mm(i,j)=0;            
                    end                    
                end
                
                delete=equation_line2(x(:,j),y(:,j),0,period_x*1/sqrt(3),-0.5*period_x,1/sqrt(3)/2*period_x,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_mm(i,j)=0;            
                    end                    
                end
                               
                delete=equation_line2(x(:,j),y(:,j),0,-1/sqrt(3)*period_x,-0.5*period_x,-1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_mm(i,j)=0;            
                    end                    
                end
                               
                delete=equation_line2(x(:,j),y(:,j),0,-1/sqrt(3)*period_x,0.5*period_x,-1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_mm(i,j)=0;            
                    end                    
                end
                               
                %%}
                
                % z_dl        
                %%{        
                delete=equation_line2(x(:,j),y(:,j),0,-1/sqrt(3)*period_x,-0.5*period_x,-1/sqrt(3)/2*period_x,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_dl(i,j)=0;            
                    end                    
                end
                               
                
                delete=equation_line2(x(:,j),y(:,j),0,0,0,-1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_dl(i,j)=0;            
                    end                    
                end
                               
                %%}        
        
                % z_dr        
                %%{        
                delete=equation_line2(x(:,j),y(:,j),0,-1/sqrt(3)*period_x,0.5*period_x,-1/sqrt(3)/2*period_x,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_dr(i,j)=0;            
                    end                    
                end               
                
                delete=equation_line2(x(:,j),y(:,j),0,0,0,-1/sqrt(3)/2*period_x,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_dr(i,j)=0;            
                    end                    
                end               
                %%}        
        
                % z_ur        
                %%{        
                delete=equation_line2(x(:,j),y(:,j),0,1/sqrt(3)*period_x,0.5*period_x,1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_ur(i,j)=0;            
                    end                    
                end               
                
                delete=equation_line2(x(:,j),y(:,j),0,0,0,-1/sqrt(3)/2*period_x,1);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_ur(i,j)=0;            
                    end                    
                end                
                %%}        
        
                % z_ul        
                %%{        
                delete=equation_line2(x(:,j),y(:,j),0,1/sqrt(3)*period_x,-0.5*period_x,1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_ul(i,j)=0;            
                    end                    
                end               
                
                delete=equation_line2(x(:,j),y(:,j),0,0,0,-1/sqrt(3)/2*period_x,2);        
                
                for i=1:1:length(delete)            
                    if delete(i)==1                
                        z_ul(i,j)=0;            
                    end                    
                end                
                %%}    
            %end            
        end        
        
        z=z_mm+z_dl+z_dr+z_ul+z_ur;

        [px,py] = gradient(z,step,step);
        
        
        N_x_input=px;%fliplr(N_x_input);
        N_y_input=py;
end

%savefile = [save_as '.mat'];
%save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat');