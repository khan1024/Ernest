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
    x_rozliseni=perioda/resolution_x;% rozlišeni x
    z_rozliseni=abs(thickness_1)/(size_region_1*resolution_z);% rozlišení z
    z=thickness_1:z_rozliseni:0; %pohybuju se uprostřed l-té vrstvičky
    z_11=thickness_1:z_rozliseni:0;
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

        u_x=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
        u_y=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);
        u_z=-cos(psi)*sin(theta);
        
        R_z_n=(R_x.*k_x + R_y.*k_y)./k_1_z_n';
        
        %
        y=0;
        for z_1=1:1:length(z)
            for x_1=1:1:length(x)
                if conical_field_component==1
                    H_x(x_1,z_1)=1/k_0*sqrt(epsilon/mu)*((-k_0*n1*sin(theta)*sin(phi))*u_z.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))- ...
                        (-k_y).*R_z_n'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)))+(k_0*n1*cos(theta))*u_y.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))+ ...
                        k_1_z_n.*R_y'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1))));
                elseif conical_field_component==2
                    H_y(x_1,z_1)=1/k_0*sqrt(epsilon/mu)*((-k_0*n1*cos(theta))*u_x.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))- ...
                        k_1_z_n.*R_x'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)))+(k_0*n1*sin(theta)*cos(phi))*u_z.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))- ...
                        k_x'.*R_z_n'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1))));
                elseif conical_field_component==3
                    H_z(x_1,z_1)=1/k_0*sqrt(epsilon/mu)*((k_0*n1*sin(theta)*cos(phi))*u_y.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))+ ...
                        k_x'.*R_y'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)))-(k_0*n1*sin(theta)*sin(phi))*u_x.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))- ...
                        k_y.*R_x'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1))));
                elseif conical_field_component==4
                    E_x(x_1,z_1)=u_x.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))+R_x'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)));
                elseif conical_field_component==5
                    E_y(x_1,z_1)=u_y.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))+R_y'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)));
                elseif conical_field_component==6
                    E_z(x_1,z_1)=u_z.*exp(-1i*k_0*n1*(sin(theta)*cos(phi)*x(x_1)+sin(theta)*sin(phi)*y-cos(theta)*z(z_1)))-R_z_n'*exp(-1i*(k_x*x(x_1)+k_y*y+k_1_z'*z(z_1)));
                end
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
            for x_1=1:1:length(x)
                if vykresleni_urciteho_radu==1
                    if planar_field_component==1
                        E_x(x_1,z_1)=sqrt(mu/epsilon)*cos(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+sqrt(mu/epsilon)*(1/(k_0*n1^2))*k_1_z_n(rad_index)*R(rad_index)*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)*z(z_1)));
                    elseif planar_field_component==2
                        H_y(x_1,z_1)=exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+R(rad_index)*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)*z(z_1)));
                    elseif planar_field_component==3
                        E_z(x_1,z_1)=1i*sqrt(mu/epsilon)*sin(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))-sqrt(mu/epsilon)*(1/(k_0*n1^2))*k_x(rad_index)*R(rad_index)*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)*z(z_1)));
                    end
                elseif vykresleni_urciteho_radu==2
                    if planar_field_component==1
                        E_x(x_1,z_1)=sqrt(mu/epsilon)*cos(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+sqrt(mu/epsilon)*(1/(n1^2*k_0))*k_1_z_n.*R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1)));
                    elseif planar_field_component==2
                        H_y(x_1,z_1)=exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1)));
                    elseif planar_field_component==3
                        E_z(x_1,z_1)=-(sqrt(mu/epsilon)*sin(theta)/(n1)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+sqrt(mu/epsilon)*(1/(k_0*n1^2))*k_x'.*R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1))));
                    end
                end
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
            for x_1=1:1:length(x)
                if vykresleni_urciteho_radu==1
                    if planar_field_component==1
                        H_x(x_1,z_1)=sqrt(epsilon/mu)*(-n1)*cos(theta)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))- ...
                            sqrt(epsilon/mu)/k_0*k_1_z_n(rad_index).*R(rad_index)'*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)'*z(z_1)));
                    elseif planar_field_component==2
                        E_y(x_1,z_1)=exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+R(rad_index)*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)*z(z_1)));
                    elseif planar_field_component==3
                        H_z(x_1,z_1)=-sqrt(epsilon/mu)*n1*sin(theta)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))- ...
                                     sqrt(epsilon/mu)*(1/k_0)*k_x(rad_index)*R(rad_index)*exp(-1i*(k_x(rad_index)*x(x_1)+k_1_z(rad_index)*z(z_1)));
                    end
                %---------------------------------------------------------
                elseif vykresleni_urciteho_radu==2
                    if planar_field_component==1
                        H_x(x_1,z_1)=sqrt(epsilon/mu)*(-n1)*cos(theta)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))- ...
                            sqrt(epsilon/mu)/k_0*k_1_z_n.*R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1)));
                    elseif planar_field_component==2
                        E_y(x_1,z_1)=exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))+R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1)));
                    elseif planar_field_component==3
                        H_z(x_1,z_1)=-sqrt(epsilon/mu)*n1*sin(theta)*exp(-1i*k_0*n1*(sin(theta)*x(x_1)-cos(theta)*z(z_1)))- ...
                                     sqrt(epsilon/mu)*(1/k_0)*k_x'.*R'*exp(-1i*(k_x*x(x_1)+k_1_z'*z(z_1)));
                    end
                end
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