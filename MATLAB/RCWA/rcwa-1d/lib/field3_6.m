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
    x_res_=perioda/resolution_x;% rozlišeni x
    z_res_=abs(thickness_3)/(size_region_3*resolution_z);% rozlišení z
    z=(D_l+z_res_):z_res_:(D_l+thickness_3); %pohybuju se uprostřed l-té vrstvičky
    z_33=z;
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
    
    if polarization==0 % kónická difrakce
        
        T_z_n=-(T_x.*k_x + T_y.*k_y)./k_3_z_n';
        
        if sum(isnan(T_z_n))>=1
            sprintf('this is bug, change theta slightly')
        end
        
        y=0;
        
        for z_3=1:1:length(z)
            exp_field=exp(-1i*(bsxfun(@minus,k_x*x,k_3_z*(D_l-z(z_3)))));
        
            if conical_field_component==1
                H_x(:,z_3)=-1/k_0*(conj(k_y.*T_z_n)'+conj(k_3_z.*T_y)')*exp_field;
            elseif conical_field_component==2
                H_y(:,z_3)=1/k_0/n1*(conj(k_3_z.*T_x)'+conj(k_x.*T_z_n)')*exp_field;
            elseif conical_field_component==3 
                H_z(:,z_3)=-1/k_0*(conj(k_y.*T_x)'-conj(k_x.*T_y)')*exp_field;
            elseif conical_field_component==4 
                E_x(:,z_3)=conj(T_x)'*exp_field; % 1/n1
            elseif conical_field_component==5  
                E_y(:,z_3)=conj(T_y)'*exp_field;
            elseif conical_field_component==6
                E_z(:,z_3)=-conj(T_z_n)'*exp_field;
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
            exp_field=exp(-1i*(bsxfun(@minus,k_x*x,k_3_z*(D_l-z(z_3)))));
            if planar_field_component==1
                E_x(:,z_3)=sqrt(mu/epsilon)*(1/(k_0*n3^2))*conj(k_3_z.*T)'*exp_field;
            elseif planar_field_component==2
                H_y(:,z_3)=conj(T')*exp_field;
            elseif planar_field_component==3
                E_z(:,z_3)=sqrt(mu/epsilon)*1/(k_0*n3^2)*(conj(k_x.*T)'*exp_field);
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
            exp_field=exp(-1i*(bsxfun(@minus,k_x*x,k_3_z*(D_l-z(z_3)))));
            if planar_field_component==1                        
                H_x(:,z_3)=-sqrt(epsilon/mu)/k_0*conj(k_3_z.*T)'*exp_field;                    
            elseif planar_field_component==2                        
                E_y(:,z_3)=conj(T)'*exp_field;                   
            elseif planar_field_component==3                        
                H_z(:,z_3)=sqrt(epsilon/mu)*(1/k_0)*conj(k_x.*T)'*exp_field;                    
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