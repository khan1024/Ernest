
switch_permittivity=0; % do not change, gratings 0-9
switch grating
    case 0 % binary grating
        number_of_layers=1;
        thickness=thickness_total*1E-6;
        layer_thickness=thickness;
    case 1 % blazed grating, type 1 <
        blazed_grating_1;
    case 2 % blazed grating, typ 2 >
        blazed_grating_2;
    case 3 % sinusoidal grating
        sinusoidal_grating_y;
    case 4 % triangular grating
        triangular_grating;
    case 5 % parabolic grating
        parabolic_grating;
    case 6 % general triangular grating
        general_triangular_grating;
    case 7 % trapezoidal grating
        trapezoidal_grating;
    case 8 % general tetragonal grating
        tetragonal_grating;
    case 9 % gaussian grating
        gaussian_grating;
  
    case 10 % grating from file
            
        if use_dispersion==1                
            setup_dispersion;            
        end        
        
        how_many_layers;            
        if grating_relative_dimensions==2                
            thickness_total=thickness_total_a;                
            Lambda=Lambda_um;            
        end        
        
        thickness=thickness_total*1E-6;            
        switch_permittivity=1;
     case 11
        Kat_grating_1;

end

if see_diffraction_layout==1
    diffraction_angles=diffraction_layout(lambda,n1,n3,Lambda,theta0);
end