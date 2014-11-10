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

%{
Lambda_x_old=Lambda_x;
Lambda_y_old=Lambda_y;

% 1) Lambda_x
fix_Lambda_x=fix(Lambda_x/dx);
Lambda_x=fix_Lambda_x*dx;

% 1) Lambda_y
fix_Lambda_y=fix(Lambda_y/dx);
Lambda_y=fix_Lambda_y*dx;
%}

object(:,:,1)=[Lambda_x/2,Lambda_y/2,r_cylinder,n_cylinder];

save_as='rectangular_grating_mat';

p_x=linspace(0,Lambda_x,resolution_dx+1);
p_y=linspace(0,Lambda_y,resolution_dx+1);

%p_x=0:dx:Lambda_x;
%p_y=0:dx:Lambda_y;

mesh_n_mat=zeros(length(p_y)-1,length(p_x)-1)+n_layer;
pixel_center=zeros(length(p_y)-1,length(p_x)-1,2); % 1-x coordinatex, 2-y coordinate

% nalézt středy jednotlivých pixelů
for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        pixel_center(i,j,1)=(p_x(j)+p_x(j+1))/2;
        pixel_center(i,j,2)=(p_y(i)+p_y(i+1))/2;
    end
end

% pokud střed pixelu bude uvnitř nejakého prvku, permitivita pixelu se
% označí podle permitivity prvku

[m,n,number_of_objects] = size(object);

for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        
        x_coordinate=pixel_center(i,j,1);
        y_coordinate=pixel_center(i,j,2);
        
        for k=1:1:number_of_objects
            
            distance_origin_pixel=sqrt((x_coordinate-object(1,1,k))^2+(y_coordinate-object(1,2,k))^2);
            
            if distance_origin_pixel<=object(1,3,k)
                mesh_n_mat(i,j)=object(1,4,k);
            end
        
        end
        
    end
end
%

size_eps=size(mesh_n_mat);

n_of_rows=size_eps(1);
n_of_columns=size_eps(2);

if show_grating_mesh==1
    grating_mesh;
end

p_x_mat=p_x;
p_y_mat=p_y;

switch factorization_method
    case 0 % Li factorization
    case 1 % moje normal vector method
        [cone_x,cone_y]=meshgrid(pixel_center(1,:,1),pixel_center(:,1,2));
        cone_z=-sqrt((cone_x-Lambda_x/2).^2+(cone_y-Lambda_y/2).^2);
        
        dx=pixel_center(1,2,1)-pixel_center(1,1,1);
        dy=pixel_center(2,1,2)-pixel_center(1,1,2);
        [N_x_input,N_y_input]=gradient(cone_z,dx,dy);
        N_x_input=fliplr(N_x_input);
        N_y_input=flipud(N_y_input);
        %figure;contour(cone_x,cone_y,cone_z), hold on, quiver(cone_x,cone_y,N_x_input,N_y_input)
    case 2
        dx=pixel_center(1,2,1)-pixel_center(1,1,1);
        dy=pixel_center(2,1,2)-pixel_center(1,1,2);
        
        Nf=length(dx/2:dx:Lambda_x-dx/2); 
        
        Nx=(1:Nf)'*ones(1,Nf); 
        Ny=ones(Nf,1)*(1:Nf);
        
        phi_n=angle(Nx-Nf/2+1i*(Ny-Nf/2));
        
        th=phi_n;
        a=cos(th);
        b=sin(th);
        
        % obraceně zadané 
        N_x_input=b; %a
        N_y_input=a; %b

    case 3
        dx=pixel_center(1,2,1)-pixel_center(1,1,1);
        
        wx=r_cylinder*2/Lambda_x;
        wy=wx;
        
        Nf=length(dx/2:dx:Lambda_x-dx/2); Nx=(1:Nf)'*ones(1,Nf); Ny=ones(Nf,1)*(1:Nf);
        Nfx=floor(wx*Nf/2+0.5); Nfy=floor(wy*Nf/2+0.5);
        
        phi_n=angle(Nx-Nf/2+1i*(Ny-Nf/2));
        rr=sqrt((Nx-Nf/2).^2+(Ny-Nf/2).^2);
        th=phi_n;
        mcs = MaxCosSin(phi_n);
        elliOut = (pi/8)*(1 + cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        elliIn = (pi/8)*(1+cos(pi*rr/Nfx));
        elli = elliOut + (elliIn - elliOut).*(((Nx-Nf/2).^2/Nfx^2 + (Ny-Nf/2).^2/Nfy^2)<=1);
        
        a = exp(1i*th).*(cos(th).*cos(elli) - 1i*sin(th).*sin(elli));
        b = exp(1i*th).*(sin(th).*cos(elli) + 1i*cos(th).*sin(elli));
        
        % obraceně zadané 
        N_x_input=b; %a
        N_y_input=a; %b
    case 4
        wx=r_cylinder*2/Lambda_x;
        wy=wx;
        
        dx=pixel_center(1,2,1)-pixel_center(1,1,1);
        
        Nf=length(dx/2:dx:Lambda_x-dx/2); Nx=(1:Nf)'*ones(1,Nf); Ny=ones(Nf,1)*(1:Nf);
        Nfx=floor(wx*Nf/2+0.5); Nfy=floor(wy*Nf/2+0.5);
        
        phi_n=angle(Nx-Nf/2+i*(Ny-Nf/2));
        rr=sqrt((Nx-Nf/2).^2+(Ny-Nf/2).^2);
        
        mcs = MaxCosSin(phi_n);
        thB = pi/2*floor(phi_n/(pi/2)+0.5);
        elliB = pi/8*(1 - cos(4*phi_n));
        
        thIn = phi_n;
        thOut = 0.5*(thB+phi_n + (thB-phi_n).*cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        
        th = thOut + (thIn - thOut).*(((Nx-Nf/2).^2/Nfx^2 + (Ny-Nf/2).^2/Nfy^2)<=1);
        
        elliIn = (pi/8)*(1+cos(pi*rr/Nfx));
        elliOut = (elliB/2).*(1 + cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        
        elli = elliOut + (elliIn - elliOut).*(((Nx-Nf/2).^2/Nfx^2 + (Ny-Nf/2).^2/Nfy^2)<=1);
        
        a = exp(1i*th).*(cos(th).*cos(elli) - 1i*sin(th).*sin(elli));
        b = exp(1i*th).*(sin(th).*cos(elli) + 1i*cos(th).*sin(elli));
        
        % obraceně zadané 
        N_x_input=b; %a
        N_y_input=a; %b
end


savefile = [save_as '.mat'];
if factorization_method==0
    save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat');
else
    save(savefile, 'p_x_mat', 'p_y_mat', 'mesh_n_mat', 'N_x_input', 'N_y_input');
end