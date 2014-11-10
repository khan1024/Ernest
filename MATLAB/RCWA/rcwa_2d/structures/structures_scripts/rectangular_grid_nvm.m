function [N_x_input,N_y_input]=rectangular_grid_nvm(Lambda_x,Lambda_y,resolution_dx,factorization_method,r_cylinder)

p_x=linspace(0,Lambda_x,resolution_dx+1);
p_y=linspace(0,Lambda_y,resolution_dx+1);

%p_x=0:dx:Lambda_x;
%p_y=0:dx:Lambda_y;

pixel_center=zeros(length(p_y)-1,length(p_x)-1,2); % 1-x coordinatex, 2-y coordinate

center_px=(p_x(1:end-1)+p_x(2:end))/2;
center_py=(p_y(1:end-1)+p_y(2:end))/2;

for i=1:1:length(p_y)-1
    pixel_center(i,:,1)=center_px;
end

for j=1:1:length(p_x)-1
    pixel_center(:,j,2)=center_py;
end

%{
% nalézt středy jednotlivých pixelů
for i=1:1:length(p_y)-1
    for j=1:1:length(p_x)-1
        pixel_center(i,j,1)=(p_x(j)+p_x(j+1))/2;
        pixel_center(i,j,2)=(p_y(i)+p_y(i+1))/2;
    end
end
%}

if mod(resolution_dx,2)==0
    % even _number
    nvm_shift=0.5;
else
    % odd number
    nvm_shift=0;
end

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
        
        phi_n=angle(Nx-Nf/2-nvm_shift+1i*(Ny-Nf/2-nvm_shift));
        
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
        
        phi_n=angle(Nx-Nf/2-nvm_shift+1i*(Ny-Nf/2-nvm_shift));
        rr=sqrt((Nx-Nf/2-nvm_shift).^2+(Ny-Nf/2-nvm_shift).^2);
        
        th=phi_n;
        mcs = MaxCosSin(phi_n);
        elliOut = (pi/8)*(1 + cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        elliIn = (pi/8)*(1+cos(pi*rr/Nfx));
        elli = elliOut + (elliIn - elliOut).*(((Nx-Nf/2-nvm_shift).^2/Nfx^2 + (Ny-Nf/2-nvm_shift).^2/Nfy^2)<=1);
        
        
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
        
        phi_n=angle(Nx-Nf/2-nvm_shift+1i*(Ny-Nf/2-nvm_shift));
        rr=sqrt((Nx-Nf/2-nvm_shift).^2+(Ny-Nf/2-nvm_shift).^2);
        
        mcs = MaxCosSin(phi_n);
        thB = pi/2*floor(phi_n/(pi/2)+0.5);
        elliB = pi/8*(1 - cos(4*phi_n));
        
        thIn = phi_n;
        thOut = 0.5*(thB+phi_n + (thB-phi_n).*cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        
        th = thOut + (thIn - thOut).*(((Nx-Nf/2-nvm_shift).^2/Nfx^2 + (Ny-Nf/2-nvm_shift).^2/Nfy^2)<=1);
        
        elliIn = (pi/8)*(1+cos(pi*rr/Nfx));
        elliOut = (elliB/2).*(1 + cos(2*pi*(rr+Nf./(2*mcs)-2*Nfx)./(Nf./mcs-2*Nfx)));
        
        elli = elliOut + (elliIn - elliOut).*(((Nx-Nf/2-nvm_shift).^2/Nfx^2 + (Ny-Nf/2-nvm_shift).^2/Nfy^2)<=1);
        
        a = exp(1i*th).*(cos(th).*cos(elli) - 1i*sin(th).*sin(elli));
        b = exp(1i*th).*(sin(th).*cos(elli) + 1i*cos(th).*sin(elli));
        
        % obraceně zadané 
        N_x_input=b; %a
        N_y_input=a; %b
end