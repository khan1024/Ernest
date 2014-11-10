% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%D \installlanguage[cz][date={day,{.},\ ,month,\ ,year},]
%D \def\mat#1{\text{\bf #1}}

%D Nejprve si musíme nadefinovat permitivitu $\varepsilon$ a permeabilitu
%D $\mu$.
mu=1.2566370614E-6;
epsilon=8.854187818E-12;

thickness=thickness_total*1e-6;
thickness_1=-size_region_1*thickness;

%D Pokud máme navoleno \type{vypocet_field=1}, tak se provede následující 
%D výpočet.
if plot_field_region_1==1
    if grating==10
        x_min=data(1,2)*1e-6;
        x_max=data(1,2)*1e-6+number_of_plotted_period*perioda;
    else
        x_min=0;
        x_max=number_of_plotted_period*perioda;% kolik period mřížky si necháme vykreslit
    end
    z_min=0;
    x_res_=perioda/resolution_x;% rozlišeni x
    z_res_=abs(thickness_1)/(size_region_1*resolution_z);% rozlišení z
    z=thickness_1:z_res_:0; %pohybuju se uprostřed l-té vrstvičky
    z_11=thickness_1:z_res_:0;
    x=x_min:x_res_:x_max;
    
    if polarization==0 %conical            
        if conical_field_component==1                
            H_x=zeros(length(x),length(z));            
        elseif conical_field_component==2                
            H_y=zeros(length(x),length(z));            
        elseif conical_field_component==3                
            H_z=zeros(length(x),length(z));            
        elseif conical_field_component==4                
            E_x=zeros(length(x),length(z));            
        elseif conical_field_component==5                
            E_y=zeros(length(x),length(z));            
        elseif conical_field_component==6                
            E_z=zeros(length(x),length(z));            
        end        
    elseif polarization==1 %TM            
        if planar_field_component==1                
            E_x=zeros(length(x),length(z));            
        elseif planar_field_component==2                
            H_y=zeros(length(x),length(z));            
        elseif planar_field_component==3                
            E_z=zeros(length(x),length(z));            
        end        
    elseif polarization==2 %TE            
        if planar_field_component==1                
            H_x=zeros(length(x),length(z));            
        elseif planar_field_component==2                
            E_y=zeros(length(x),length(z));            
        elseif planar_field_component==3                
            H_z=zeros(length(x),length(z));            
        end
    end
        

    if polarization==0 % conical diffraction

        u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
        u_y=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);
        u_z=-cos(psi)*sin(theta);
        
        R_z_n=(R_x.*k_x + R_y.*k_y)./k_1_z_n';
        
        if sum(isnan(R_z_n))>=1
            sprintf('this is bug, change theta slightly')
        end
        
        y=0;

        for z_1=1:1:length(z)
            exp_field_R=exp(-1i*(bsxfun(@minus,k_x*x,k_1_z*z(z_1))));
            exp_plane_wave=exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x+sin(theta)*sin(phi)*y+cos(theta)*z(z_1)));
            if conical_field_component==1   
                H_x(:,z_1)=-1/k_0*((-k_0*n1*sin(theta)*sin(phi))*u_z.*exp_plane_wave+ ...
                        conj(k_y.*R_z_n)'*exp_field_R+(k_0*n1*cos(theta))*u_y.*exp_plane_wave- ...
                        conj(k_1_z.*R_y)'*exp_field_R);
            elseif conical_field_component==2
                H_y(:,z_1)=-1/k_0*((-k_0*n1*cos(theta))*u_x.*exp_plane_wave+ ...
                        conj(k_1_z.*R_x)'*exp_field_R+(k_0*n1*sin(theta)*cos(phi))*u_z.*exp_plane_wave+ ...
                        conj(-k_x.*R_z_n)'*exp_field_R);
            elseif conical_field_component==3   
                H_z(:,z_1)=1/k_0*((k_0*n1*sin(theta)*cos(phi))*u_y.*exp_plane_wave+ ...
                        conj(k_x.*R_y)'*exp_field_R-(k_0*n1*sin(theta)*sin(phi))*u_x.*exp_plane_wave- ...
                        conj(k_y.*R_x)'*exp_field_R);
            elseif conical_field_component==4
                E_x(:,z_1)=u_x.*exp_plane_wave+conj(R_x)'*exp_field_R;
            elseif conical_field_component==5
                E_y(:,z_1)=u_y.*exp_plane_wave+conj(R_y)'*exp_field_R;
            elseif conical_field_component==6
                E_z(:,z_1)=u_z.*exp_plane_wave-conj(R_z_n)'*exp_field_R;
            end
        end
        
        if conical_field_component==1
            field_H_1=H_x;
        elseif conical_field_component==2
            field_H_1=H_y;
        elseif conical_field_component==3
            field_H_1=H_z;
        elseif conical_field_component==4
            field_E_1=E_x;
        elseif conical_field_component==5
            field_E_1=E_y;
        elseif conical_field_component==6
            field_E_1=E_z;
        end
    elseif polarization==1 % TM polarization

        for z_1=1:1:length(z)
            exp_field_R=exp(-1i*(bsxfun(@minus,k_x*x,k_1_z*z(z_1))));
            if planar_field_component==1
                E_x(:,z_1)=sqrt(mu/epsilon)*cos(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))-sqrt(mu/epsilon)*(1/(n1^2*k_0))*conj(k_1_z.*R)'*exp_field_R;  
            elseif planar_field_component==2     
                H_y(:,z_1)=exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))+conj(R')*exp_field_R;
            elseif planar_field_component==3
                E_z(:,z_1)=(sqrt(mu/epsilon)*sin(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))+sqrt(mu/epsilon)*(1/(k_0*n1^2))*conj(k_x.*R)'*exp_field_R);
            end
        end
        
        if planar_field_component==1
            field_E_1=E_x;
        elseif planar_field_component==2
            field_H_1=H_y;
        elseif planar_field_component==3
            field_E_1=E_z;
        end
    elseif polarization==2 % TE polarization
        
        for z_1=1:1:length(z)
            exp_field_R=exp(-1i*(bsxfun(@minus,k_x*x,k_1_z*z(z_1))));
            if planar_field_component==1
                H_x(:,z_1)=sqrt(epsilon/mu)*(-n1)*cos(theta)*exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))+ ...
                    sqrt(epsilon/mu)/k_0*conj(k_1_z.*R)'*exp_field_R;
            elseif planar_field_component==2
                E_y(:,z_1)=exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))+conj(R)'*exp_field_R;   
            elseif planar_field_component==3
                H_z(:,z_1)=sqrt(epsilon/mu)*n1*sin(theta)*exp(-1i*k_0*n1*(sin(theta)*x+cos(theta)*z(z_1)))+ ...
                    sqrt(epsilon/mu)*(1/k_0)*conj(k_x.*R)'*exp_field_R;
            end
        end
        
        if planar_field_component==1
            field_H_1=H_x;
        elseif planar_field_component==2
            field_E_1=E_y;
        elseif planar_field_component==3
            field_H_1=H_z;
        end
    end
end