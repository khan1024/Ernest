% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

%{
old_transformation=2;
 
 if old_transformation==1 
     for r=0:1:n-1
         for s=0:1:n-1;
             koef_i=(fix(r/(2*N_Y+1)))-(fix(s/(2*N_Y+1))); % koeficient i
             koef_j=(fix(rem(r,(2*N_Y+1))))-(fix(rem(s,(2*N_Y+1)))); % koeficent j         
             switch grating
                 case 1
                     fourier_rcwa2d_cylinder;
                 case 2
                     fourier_rcwa2d_cylinder_hexagonal;
             end
             E(r+1,s+1)=koeficient_E;
             A(r+1,s+1)=koeficient_A;
         end
     end
 elseif old_transformation==2;
%}
% eq[3] [[epsilon]]_(mn,jl)=
%
     
for m=1:1:n_x
    for n=1:1:n_x
        for j=1:1:n_x
            for l=1:1:n_x
                u=n_x*(m-1)+n;
                v=n_x*(j-1)+l;
                coeff_i=m-j;
                coeff_j=n-l;
                     
                switch grating
                    case 1   
                        fourier_rcwa2d_cylinder;
                    case 2   
                        fourier_rcwa2d_cylinder_hexagonal;
                end
                
                E(u,v)=coefficient_E;
                A(u,v)=coefficient_A;
            end
            
        end
        
    end
    
end
     
 %end