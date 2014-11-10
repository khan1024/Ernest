if symmetry_switch==5 % even-even mode

    if plot_field_component==1 % E_1 
        % (E^x)_(i,j)= (E^x)_(-i,j)= (E^x)_(i,-j)=(E^x)_(-i,-j),
        % (E^x)_(i,j)->4
        % (E^x)_(i,0) -> 3,
        % (E^x)_(0,j) -> 2
        
        E_1_22_h=harmonics_reconstruction_x_y_1(N_X,1,E_1_harmonics);
        E_1_23_h=harmonics_reconstruction_x_y_1(N_X,2,E_1_harmonics);
        E_1_32_h=harmonics_reconstruction_x_y_1(N_X,3,E_1_harmonics);
        E_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,E_1_harmonics);
        
        E_1_11_h=rot90(E_1_33_h,2);
        E_1_12_h=flipud(E_1_32_h);
        E_1_13_h=flipud(E_1_33_h);
        E_1_21_h=fliplr(E_1_23_h);
        E_1_31_h=fliplr(E_1_33_h);
        
        E_1_harmonics2=[E_1_11_h, E_1_12_h, E_1_13_h; E_1_21_h, E_1_22_h, E_1_23_h; E_1_31_h, E_1_32_h, E_1_33_h]; 
    elseif plot_field_component==2 % E_2
        % (E^y)_(i,j)= -(E^y)_(-i,j)= -(E^y)_(i,-j)=(E^y)_(-i,-j)
        % (E^y)_(i,j)->5
        % (E^y)_(i,0)=(E^y)_(0,j)=0
        
        E_2_22_h=0;
        E_2_23_h=zeros(1,N_X);
        E_2_32_h=zeros(N_X,1);
        E_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,E_2_harmonics);
        
        E_2_11_h=rot90(E_2_33_h,2);
        E_2_12_h=E_2_32_h;
        E_2_13_h=-flipud(E_2_33_h);
        E_2_21_h=E_2_23_h;
        E_2_31_h=-fliplr(E_2_33_h);
        
        E_2_harmonics2=[E_2_11_h, E_2_12_h, E_2_13_h; E_2_21_h, E_2_22_h, E_2_23_h; E_2_31_h, E_2_32_h, E_2_33_h];
        

    elseif plot_field_component==3 % E_3, S_3=1i*inv(eps_33)*(K_x*U_2-K_y*U_1)
        
        H_1_22_h=0;
        H_1_23_h=zeros(1,N_X);
        H_1_32_h=zeros(N_X,1);
        H_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,H_1_harmonics);
        
        H_1_11_h=rot90(H_1_33_h,2);
        H_1_12_h=H_1_32_h;
        H_1_13_h=-flipud(H_1_33_h);
        H_1_21_h=H_1_23_h;
        H_1_31_h=-fliplr(H_1_33_h);
        
        H_1_harmonics2=[H_1_11_h, H_1_12_h, H_1_13_h; H_1_21_h, H_1_22_h, H_1_23_h; H_1_31_h, H_1_32_h, H_1_33_h];
        
        H_2_22_h=harmonics_reconstruction_x_y_1(N_X,1,H_2_harmonics);
        H_2_23_h=harmonics_reconstruction_x_y_1(N_X,2,H_2_harmonics);
        H_2_32_h=harmonics_reconstruction_x_y_1(N_X,3,H_2_harmonics);
        H_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,H_2_harmonics);
        
        H_2_11_h=rot90(H_2_33_h,2);
        H_2_12_h=flipud(H_2_32_h);
        H_2_13_h=flipud(H_2_33_h);
        H_2_21_h=fliplr(H_2_23_h);
        H_2_31_h=fliplr(H_2_33_h);
        
        H_2_harmonics2=[H_2_11_h, H_2_12_h, H_2_13_h; H_2_21_h, H_2_22_h, H_2_23_h; H_2_31_h, H_2_32_h, H_2_33_h];
        
        H_1_harmonics=conj(reshape(H_1_harmonics2',n_x*n_x,1));
        H_2_harmonics=conj(reshape(H_2_harmonics2',n_x*n_x,1));
        
        %F_1_h=inv_eps_33*alpha;
        %F_2_h=inv_eps_33*beta;
        
        H_1_harmonics=F_2_h*H_1_harmonics;    
        H_1_harmonics2=reshape(H_1_harmonics',n_x,n_x)';       
       
        H_2_harmonics=F_1_h*H_2_harmonics;    
        H_2_harmonics2=reshape(H_2_harmonics',n_x,n_x)';
    
        E_3_harmonics2=1i*(H_2_harmonics2-H_1_harmonics2);
    elseif plot_field_component==4 % H_1
        % (H^x)_(i,j)= -(H^x)_(-i,j)= -(H^x)_(i,-j)=(H^x)_(-i,-j)
        % (H^x)_(i,j)->5
        % (H^x)_(i,0)=(H^x)_(0,j)=0
        
        H_1_22_h=0;
        H_1_23_h=zeros(1,N_X);
        H_1_32_h=zeros(N_X,1);
        H_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,H_1_harmonics);
        
        H_1_11_h=rot90(H_1_33_h,2);
        H_1_12_h=H_1_32_h;
        H_1_13_h=-flipud(H_1_33_h);
        H_1_21_h=H_1_23_h;
        H_1_31_h=-fliplr(H_1_33_h);
        
        H_1_harmonics2=[H_1_11_h, H_1_12_h, H_1_13_h; H_1_21_h, H_1_22_h, H_1_23_h; H_1_31_h, H_1_32_h, H_1_33_h];
    elseif plot_field_component==5 % H_2    
        % (H^y)_(i,j)= (H^y)_(-i,j)= (H^y)_(i,-j)=(H^y)_(-i,-j),
        % (H^y)_(i,j)->4
        % (H^y)_(i,0) -> 3,
        % (H^y)_(0,j) -> 2
        
        H_2_22_h=harmonics_reconstruction_x_y_1(N_X,1,H_2_harmonics);
        H_2_23_h=harmonics_reconstruction_x_y_1(N_X,2,H_2_harmonics);
        H_2_32_h=harmonics_reconstruction_x_y_1(N_X,3,H_2_harmonics);
        H_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,H_2_harmonics);
        
        H_2_11_h=rot90(H_2_33_h,2);
        H_2_12_h=flipud(H_2_32_h);
        H_2_13_h=flipud(H_2_33_h);
        H_2_21_h=fliplr(H_2_23_h);
        H_2_31_h=fliplr(H_2_33_h);
        
        H_2_harmonics2=[H_2_11_h, H_2_12_h, H_2_13_h; H_2_21_h, H_2_22_h, H_2_23_h; H_2_31_h, H_2_32_h, H_2_33_h];
        
    elseif plot_field_component==6 % H_3, U_3=1i*inv(mu_33)*(K_x*S_2-K_y*S_1) 
        
        E_1_22_h=harmonics_reconstruction_x_y_1(N_X,1,E_1_harmonics);
        E_1_23_h=harmonics_reconstruction_x_y_1(N_X,2,E_1_harmonics);
        E_1_32_h=harmonics_reconstruction_x_y_1(N_X,3,E_1_harmonics);
        E_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,E_1_harmonics);
        
        E_1_11_h=rot90(E_1_33_h,2);
        E_1_12_h=flipud(E_1_32_h);
        E_1_13_h=flipud(E_1_33_h);
        E_1_21_h=fliplr(E_1_23_h);
        E_1_31_h=fliplr(E_1_33_h);
        
        E_1_harmonics2=[E_1_11_h, E_1_12_h, E_1_13_h; E_1_21_h, E_1_22_h, E_1_23_h; E_1_31_h, E_1_32_h, E_1_33_h];
        
        E_2_22_h=0;
        E_2_23_h=zeros(1,N_X);
        E_2_32_h=zeros(N_X,1);
        E_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,E_2_harmonics);
        
        E_2_11_h=rot90(E_2_33_h,2);
        E_2_12_h=E_2_32_h;
        E_2_13_h=-flipud(E_2_33_h);
        E_2_21_h=E_2_23_h;
        E_2_31_h=-fliplr(E_2_33_h);
        
        E_2_harmonics2=[E_2_11_h, E_2_12_h, E_2_13_h; E_2_21_h, E_2_22_h, E_2_23_h; E_2_31_h, E_2_32_h, E_2_33_h];
        
        E_1_harmonics=conj(reshape(E_1_harmonics2',n_x*n_x,1));
        E_2_harmonics=conj(reshape(E_2_harmonics2',n_x*n_x,1));
        
        %G_1_h=inv_mu_33*beta;        
        %G_2_h=inv_mu_33*alpha;
        
        E_1_harmonics=G_1_h*E_1_harmonics;    
        E_1_harmonics2=reshape(E_1_harmonics',n_x,n_x)';
        
        E_2_harmonics=G_2_h*E_2_harmonics;    
        E_2_harmonics2=reshape(E_2_harmonics',n_x,n_x)';    
    
        H_3_harmonics2=1i*(E_2_harmonics2-E_1_harmonics2);
    end

elseif symmetry_switch==6 % odd-odd mode
    
    if plot_field_component==1 % E_1 
        % (E^x)_(i,j)= -(E^x)_(-i,j)= -(E^x)_(i,-j)=(E^x)_(-i,-j),
        % (E^x)_(i,j) -> 4
        % (E^x)_(i,0)=(E^x)_(0,j)=0
        
        E_1_22_h=0;
        E_1_23_h=zeros(1,N_X);
        E_1_32_h=zeros(N_X,1);
        E_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,E_1_harmonics);
        
        E_1_11_h=rot90(E_1_33_h,2);
        E_1_12_h=E_1_32_h;
        E_1_13_h=-flipud(E_1_33_h);
        E_1_21_h=E_1_23_h;
        E_1_31_h=-fliplr(E_1_33_h);
        
        E_1_harmonics2=[E_1_11_h, E_1_12_h, E_1_13_h; E_1_21_h, E_1_22_h, E_1_23_h; E_1_31_h, E_1_32_h, E_1_33_h];
    elseif plot_field_component==2 % E_2
        % (E^y)_(i,j)= (E^y)_(-i,j)= (E^y)_(i,-j)=(E^y)_(-i,-j)
        % (E^y)_(i,j)->5
        % (E^y)_(0,0) ->1, (E^y)_(i,0) -> 3,
        % (E^y)_(0,j) -> 2
        
        E_2_22_h=harmonics_reconstruction_x_y_1(N_X,1,E_2_harmonics);
        E_2_23_h=harmonics_reconstruction_x_y_1(N_X,2,E_2_harmonics);
        E_2_32_h=harmonics_reconstruction_x_y_1(N_X,3,E_2_harmonics);
        E_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,E_2_harmonics);
        
        E_2_11_h=rot90(E_2_33_h,2);
        E_2_12_h=flipud(E_2_32_h);
        E_2_13_h=flipud(E_2_33_h);
        E_2_21_h=fliplr(E_2_23_h);
        E_2_31_h=fliplr(E_2_33_h);
        
        E_2_harmonics2=[E_2_11_h, E_2_12_h, E_2_13_h; E_2_21_h, E_2_22_h, E_2_23_h; E_2_31_h, E_2_32_h, E_2_33_h]; 

    elseif plot_field_component==3 % E_3, S_3=1i*inv(eps_33)*(K_x*U_2-K_y*U_1) 
        
        H_1_22_h=harmonics_reconstruction_x_y_1(N_X,1,H_1_harmonics);
        H_1_23_h=harmonics_reconstruction_x_y_1(N_X,2,H_1_harmonics);
        H_1_32_h=harmonics_reconstruction_x_y_1(N_X,3,H_1_harmonics);
        H_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,H_1_harmonics);
        
        H_1_11_h=rot90(H_1_33_h,2);
        H_1_12_h=flipud(H_1_32_h);
        H_1_13_h=flipud(H_1_33_h);
        H_1_21_h=fliplr(H_1_23_h);
        H_1_31_h=fliplr(H_1_33_h);
        
        H_1_harmonics2=[H_1_11_h, H_1_12_h, H_1_13_h; H_1_21_h, H_1_22_h, H_1_23_h; H_1_31_h, H_1_32_h, H_1_33_h];
        
        H_2_22_h=0;
        H_2_23_h=zeros(1,N_X);
        H_2_32_h=zeros(N_X,1);
        H_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,H_2_harmonics);
        
        H_2_11_h=rot90(H_2_33_h,2);
        H_2_12_h=H_2_32_h;
        H_2_13_h=-flipud(H_2_33_h);
        H_2_21_h=H_2_23_h;
        H_2_31_h=-fliplr(H_2_33_h);
        
        H_2_harmonics2=[H_2_11_h, H_2_12_h, H_2_13_h; H_2_21_h, H_2_22_h, H_2_23_h; H_2_31_h, H_2_32_h, H_2_33_h];
        
        H_1_harmonics=conj(reshape(H_1_harmonics2',n_x*n_x,1));
        H_2_harmonics=conj(reshape(H_2_harmonics2',n_x*n_x,1));
        
        %F_1_h=inv_eps_33*alpha;
        %F_2_h=inv_eps_33*beta;
        
        H_1_harmonics=F_2_h*H_1_harmonics;    
        H_1_harmonics2=reshape(H_1_harmonics',n_x,n_x)';       
       
        H_2_harmonics=F_1_h*H_2_harmonics;    
        H_2_harmonics2=reshape(H_2_harmonics',n_x,n_x)';
    
        E_3_harmonics2=1i*(H_2_harmonics2-H_1_harmonics2);
    elseif plot_field_component==4 % H_1
       % (H^x)_(i,j)= (H^x)_(-i,j)= (H^x)_(i,-j)=(H^x)_(-i,-j)
        % (H^x)_(i,j)->5
        % (H^x)_(0,0) ->1, (H^x)_(i,0) -> 3,
        % (H^x)_(0,j) -> 2
        
        H_1_22_h=harmonics_reconstruction_x_y_1(N_X,1,H_1_harmonics);
        H_1_23_h=harmonics_reconstruction_x_y_1(N_X,2,H_1_harmonics);
        H_1_32_h=harmonics_reconstruction_x_y_1(N_X,3,H_1_harmonics);
        H_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,H_1_harmonics);
        
        H_1_11_h=rot90(H_1_33_h,2);
        H_1_12_h=flipud(H_1_32_h);
        H_1_13_h=flipud(H_1_33_h);
        H_1_21_h=fliplr(H_1_23_h);
        H_1_31_h=fliplr(H_1_33_h);
        
        H_1_harmonics2=[H_1_11_h, H_1_12_h, H_1_13_h; H_1_21_h, H_1_22_h, H_1_23_h; H_1_31_h, H_1_32_h, H_1_33_h]; 
    elseif plot_field_component==5 % H_2    
        % (H^y)_(i,j)= -(H^y)_(-i,j)= -(H^y)_(i,-j)=(H^y)_(-i,-j),
        % (H^y)_(i,j) -> 4
        % (H^y)_(i,0)=(H^y)_(0,j)=0
        
        H_2_22_h=0;
        H_2_23_h=zeros(1,N_X);
        H_2_32_h=zeros(N_X,1);
        H_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,H_2_harmonics);
        
        H_2_11_h=rot90(H_2_33_h,2);
        H_2_12_h=H_2_32_h;
        H_2_13_h=-flipud(H_2_33_h);
        H_2_21_h=H_2_23_h;
        H_2_31_h=-fliplr(H_2_33_h);
        
        H_2_harmonics2=[H_2_11_h, H_2_12_h, H_2_13_h; H_2_21_h, H_2_22_h, H_2_23_h; H_2_31_h, H_2_32_h, H_2_33_h];
        
    elseif plot_field_component==6 % H_3, U_3=1i*inv(mu_33)*(K_x*S_2-K_y*S_1) 
        
        E_1_22_h=0;
        E_1_23_h=zeros(1,N_X);
        E_1_32_h=zeros(N_X,1);
        E_1_33_h=harmonics_reconstruction_x_y_1(N_X,4,E_1_harmonics);
        
        E_1_11_h=rot90(E_1_33_h,2);
        E_1_12_h=E_1_32_h;
        E_1_13_h=-flipud(E_1_33_h);
        E_1_21_h=E_1_23_h;
        E_1_31_h=-fliplr(E_1_33_h);
        
        E_1_harmonics2=[E_1_11_h, E_1_12_h, E_1_13_h; E_1_21_h, E_1_22_h, E_1_23_h; E_1_31_h, E_1_32_h, E_1_33_h];
        
        E_2_22_h=harmonics_reconstruction_x_y_1(N_X,1,E_2_harmonics);
        E_2_23_h=harmonics_reconstruction_x_y_1(N_X,2,E_2_harmonics);
        E_2_32_h=harmonics_reconstruction_x_y_1(N_X,3,E_2_harmonics);
        E_2_33_h=harmonics_reconstruction_x_y_1(N_X,5,E_2_harmonics);
        
        E_2_11_h=rot90(E_2_33_h,2);
        E_2_12_h=flipud(E_2_32_h);
        E_2_13_h=flipud(E_2_33_h);
        E_2_21_h=fliplr(E_2_23_h);
        E_2_31_h=fliplr(E_2_33_h);
        
        E_2_harmonics2=[E_2_11_h, E_2_12_h, E_2_13_h; E_2_21_h, E_2_22_h, E_2_23_h; E_2_31_h, E_2_32_h, E_2_33_h];
        
        E_1_harmonics=conj(reshape(E_1_harmonics2',n_x*n_x,1));
        E_2_harmonics=conj(reshape(E_2_harmonics2',n_x*n_x,1));
        
        %G_1_h=inv_mu_33*beta;        
        %G_2_h=inv_mu_33*alpha;
        
        E_1_harmonics=G_1_h*E_1_harmonics;    
        E_1_harmonics2=reshape(E_1_harmonics',n_x,n_x)';
        
        E_2_harmonics=G_2_h*E_2_harmonics;    
        E_2_harmonics2=reshape(E_2_harmonics',n_x,n_x)';    
    
        H_3_harmonics2=1i*(E_2_harmonics2-E_1_harmonics2);
    end
end