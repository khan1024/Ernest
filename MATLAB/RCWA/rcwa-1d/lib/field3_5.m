% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%D Nejprve si musíme nadefinovat permitivitu $\varepsilon$ a permeabilitu
%D $\mu$.
mu=1.2566370614E-6;
epsilon=8.854187818E-12;
thickness=thickness_total*1e-6; %patch
thickness_3=size_region_3*thickness;
D_l=thickness;

%D Pokud máme navoleno \type{vypocet_field=1}, tak se provede následující 
%D výpočet.
if plot_field_region_3==1
    if grating==10
        x_min=data(1,2)*1e-6;
        x_max=data(1,2)*1e-6+number_of_plotted_period*perioda;
    else
        x_min=0;
        x_max=number_of_plotted_period*perioda;% kolik period mřížky si necháme vykreslit
    end
    x_rozliseni=perioda/resolution_x;% rozlišeni x
    z_rozliseni=abs(thickness_3)/(size_region_3*resolution_z);% rozlišení z
    z=(D_l+z_rozliseni):z_rozliseni:(D_l+thickness_3); %pohybuju se uprostřed l-té vrstvičky
    z_33=z;
    x=x_min:x_rozliseni:x_max;
    
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
    
    if polarization==0 % kónická difrakce
        
        T_z_n=-(T_x.*k_x + T_y.*k_y)./k_3_z_n';
        
        y=0;
        for z_3=1:1:length(z)
            for x_3=1:1:length(x)
                if conical_field_component==1
                    H_x(x_3,z_3)=1/k_0*sqrt(epsilon/mu)*(-k_y*T_z'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)))+-k_3_z_n.*T_y'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l))));
                elseif conical_field_component==2
                    H_y(x_3,z_3)=-1/k_0*sqrt(epsilon/mu)*(-k_3_z_n.*T_x'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)))+k_x'.*T_z_n'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l))));
                elseif conical_field_component==3
                    H_z(x_3,z_3)=-1/k_0*sqrt(epsilon/mu)*(k_y.*T_x'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)))-k_x'.*T_y'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l))));
                elseif conical_field_component==4
                    E_x(x_3,z_3)=T_x'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)));
                elseif conical_field_component==5
                    E_y(x_3,z_3)=T_y'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)));
                elseif conical_field_component==6
                    E_z(x_3,z_3)=-T_z_n'*exp(-1i*(k_x*x(x_3)+k_y*y-k_3_z'*(z(z_3)-D_l)));
                end
            end
        end
        if conical_field_component==1
            field_H_3=H_x;
        elseif conical_field_component==2
            field_H_3=H_y;
        elseif conical_field_component==3
            field_H_3=H_z;
        elseif conical_field_component==4
            field_E_3=E_x;
        elseif conical_field_component==5
            field_E_3=E_y;
        elseif conical_field_component==6
            field_E_3=E_z;
        end
    elseif polarization==1 %TM
        for z_3=1:1:length(z)
            for x_3=1:1:length(x)
                if vykresleni_urciteho_radu==1
                    if planar_field_component==1
                        E_x(x_3,z_3)=-sqrt(mu/epsilon)*(1/(k_0*n3^2))*k_3_z_n(rad_index)*T(rad_index)*exp(-1i*(k_x(rad_index)*x(x_3)-k_3_z(rad_index)*(z(z_3)-D_l)));
                    elseif planar_field_component==2
                        H_y(x_3,z_3)=T(rad_index)*exp(-1i*(k_x(rad_index)*x(x_3)-k_3_z(rad_index)*(z(z_3)-D_l)));
                    elseif planar_field_component==3
                        E_z(x_3,z_3)=-sqrt(mu/epsilon)*1/(k_0*n3^2)*(k_x(rad_index)*T(rad_index)*exp(-1i*(k_x(rad_index)*x(x_3)-k_3_z(rad_index)*(z(z_3)-D_l))));
                    end
                elseif vykresleni_urciteho_radu==2
                    if planar_field_component==1
                        E_x(x_3,z_3)=-sqrt(mu/epsilon)*(1/(k_0*n3^2))*k_3_z_n.*T'*exp(-1i*(k_x*x(x_3)-k_3_z'*(z(z_3)-D_l)));
                    elseif planar_field_component==2
                        H_y(x_3,z_3)=T'*exp(-1i*(k_x*x(x_3)-(k_3_z'*(z(z_3)-D_l))));
                    elseif planar_field_component==3
                        E_z(x_3,z_3)=-sqrt(mu/epsilon)*1/(k_0*n3^2)*(k_x'.*T'*exp(-1i*(k_x*x(x_3)-k_3_z'*(z(z_3)-D_l))));
                    end
                end
            end
        end
        if planar_field_component==1
            field_E_3=E_x;
        elseif planar_field_component==2
            field_H_3=H_y;
        elseif planar_field_component==3
            field_E_3=E_z;
        end
    elseif polarization==2 %TE polarization

        for z_3=1:1:length(z)
            for x_3=1:1:length(x)
                if vykresleni_urciteho_radu==1
                    if planar_field_component==1
                        H_x(x_3,z_3)=sqrt(epsilon/mu)*(1/k_0)*k_3_z_n(rad_index)*T(rad_index)*exp(-1i*(k_x(rad_index)*x(x_3)+k_3_z(rad_index)*(z(z_3)-D_l)));
                    elseif planar_field_component==2
                        E_y(x_3,z_3)=T(rad_index)*exp(-1i*(k_x(rad_index)*x(x_3)+k_3_z(rad_index)*(z(z_3)-D_l)));
                    elseif planar_field_component==3
                        H_z(x_3,z_3)=sqrt(epsilon/mu)*(1/k_0)*k_x'.*T'*exp(-1i*(k_x*x(x_3)-(k_3_z'*(z(z_3)-D_l))));
                    end
                elseif vykresleni_urciteho_radu==2
                    if planar_field_component==1
                        H_x(x_3,z_3)=sqrt(epsilon/mu)/k_0*k_3_z_n.*T'*exp(-1i*(k_x*x(x_3)-(k_3_z'*(z(z_3)-D_l))));
                    elseif planar_field_component==2
                        E_y(x_3,z_3)=T'*exp(-1i*(k_x*x(x_3)-k_3_z'*(z(z_3)-D_l)));
                    elseif planar_field_component==3
                        H_z(x_3,z_3)=-sqrt(epsilon/mu)*(1/k_0)*k_x'.*T'*exp(-1i*(k_x*x(x_3)-(k_3_z'*(z(z_3)-D_l))));
                    end
                end
            end
        end
        

        if planar_field_component==1
            field_H_3=H_x;
        elseif planar_field_component==2
            field_E_3=E_y;
        elseif planar_field_component==3
            field_H_3=H_z;
        end
    end
end