% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

% P. Lalanne and M.P. Jurek, J. Mod. Opt. 45, 1357-1374 (1998)

%D Pro normální mřížku
if switch_permittivity==0 % rozvoj permitivity pomocí Sincu
    E_x(x_l,z_l)=D_x(x_l,z_l)./(epsilon*epsg); % vyplním celé field touto hodnotou
    for i=0:1:number_of_plotted_period-1 % pro vykreslení více period
        if x_l>=(i+shift(l)-duty_cycle(l)/2)*resolution_x && x_l<=(i+shift(l)+duty_cycle(l)/2)*resolution_x
            E_x(x_l,z_l)=D_x(x_l,z_l)./(epsilon*epsr); % pouze vybrané intervaly přeměním na tuto hodnotu
        end
    end
end

%D Pro mřížku ze souboru
if grating==10
    if  grating_relative_dimensions==1
        if switch_permittivity==1 % přímo z vypočítaného integrálu 
            for p=0:1:number_of_plotted_period-1 % pro vykreslení více period 
                for i=3:2:size_data_x % projíždím zadané hodnoty 
                    if isnan(data(l,index(i)))==1 % hlídím si konec dat, ty s NaN vynechám 
                    elseif perioda*(p+data(l,index(i-1)))<=x(x_l) && perioda*(p+data(l,index(i+1)))>=x(x_l) 
                        E_x(x_l,z_l)=D_x(x_l,z_l)./(epsilon*(data(l,index(i))^2)); % podle intervalů dělím příslušnou permitivitou
                    end
                end
            end
        end
    elseif  grating_relative_dimensions==2
        if switch_permittivity==1 % přímo z vypočítaného integrálu 
            for p=0:1:number_of_plotted_period-1 % pro vykreslení více period   
                for i=3:2:size_data_x % projíždím zadané hodnoty 
                    if isnan(data(l,index(i)))==1 % hlídím si konec dat, ty s NaN vynechám  
                    elseif perioda*(p+data(l,index(i-1))/Lambda_um)<=x(x_l) && perioda*(p+data(l,index(i+1))/Lambda_um)>=x(x_l)
                        E_x(x_l,z_l)=D_x(x_l,z_l)./(epsilon*(data(l,index(i))^2)); % podle intervalů dělím příslušnou permitivitou
                    end
                end
            end
        end
    end
end