if symmetry_switch==3 % even mode

    % (E^y)_(i,0)=0, (H^x)_(i,0)=0
    %( (E^x)_(i,0); (E^x)_(i,j); (E^y)_(i,j) )
    %( (H^y)_(i,0); (H^x)_(i,j); (H^y)_(i,j) )

    if plot_field_component==1 % E_1 
        % (E^x)_(i,j)=(E^x)_(i,-j), (E^x)_(i,j) -> 2, (E^x)_(i,0) -> 1
        
        E_1_harmonics_right=harmonics_reconstruction_y(N_X,2,E_1_harmonics);    
        E_1_harmonics_left=fliplr(E_1_harmonics_right);
        E_1_harmonics_zero=harmonics_reconstruction_y(N_X,1,E_1_harmonics);
        E_1_harmonics2=[E_1_harmonics_left,E_1_harmonics_zero,E_1_harmonics_right];
    elseif plot_field_component==2 % E_2
        % (E^y)_(i,0)=0, (E^y)_(i,j)=-(E^y)_(i,-j), (E^y)_(i,j) -> 3
        
        E_2_harmonics_right=harmonics_reconstruction_y(N_X,3,E_2_harmonics);    
        E_2_harmonics_left=-fliplr(E_2_harmonics_right);    
        E_2_harmonics_zero=zeros(2*N_X+1,1);    
        E_2_harmonics2=[E_2_harmonics_left,E_2_harmonics_zero,E_2_harmonics_right];   

    elseif plot_field_component==3 % E_3, S_3=1i*inv(eps_33)*(K_x*U_2-K_y*U_1)
        
        H_1_harmonics_right=harmonics_reconstruction_y(N_X,2,H_1_harmonics);    
        H_1_harmonics_left=-fliplr(H_1_harmonics_right);    
        H_1_harmonics_zero=zeros(2*N_X+1,1);    
        H_1_harmonics2=[H_1_harmonics_left,H_1_harmonics_zero,H_1_harmonics_right];
        
        H_2_harmonics_right=harmonics_reconstruction_y(N_X,3,H_2_harmonics);    
        H_2_harmonics_left=fliplr(H_2_harmonics_right);    
        H_2_harmonics_zero=harmonics_reconstruction_y(N_X,1,H_2_harmonics);
        H_2_harmonics2=[H_2_harmonics_left,H_2_harmonics_zero,H_2_harmonics_right];
        
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
        % (H^x)_(i,0)=0, (H^x)_(i,j)=-(H^x)_(i,-j), (H^x)_(i,j) -> 2
    
        H_1_harmonics_right=harmonics_reconstruction_y(N_X,2,H_1_harmonics);    
        H_1_harmonics_left=-fliplr(H_1_harmonics_right);    
        H_1_harmonics_zero=zeros(2*N_X+1,1);    
        H_1_harmonics2=[H_1_harmonics_left,H_1_harmonics_zero,H_1_harmonics_right];   
        
    elseif plot_field_component==5 % H_2    
        % (H^y)_(i,j)=(H^y)_(i,-j), (H^y)_(i,j) -> 3, (H^y)_(i,0) -> 1
        
        H_2_harmonics_right=harmonics_reconstruction_y(N_X,3,H_2_harmonics);    
        H_2_harmonics_left=fliplr(H_2_harmonics_right);    
        H_2_harmonics_zero=harmonics_reconstruction_y(N_X,1,H_2_harmonics);
        H_2_harmonics2=[H_2_harmonics_left,H_2_harmonics_zero,H_2_harmonics_right];
        
    elseif plot_field_component==6 % H_3, U_3=1i*inv(mu_33)*(K_x*S_2-K_y*S_1)
        
        E_1_harmonics_right=harmonics_reconstruction_y(N_X,2,E_1_harmonics);    
        E_1_harmonics_left=fliplr(E_1_harmonics_right);
        E_1_harmonics_zero=harmonics_reconstruction_y(N_X,1,E_1_harmonics);
        E_1_harmonics2=[E_1_harmonics_left,E_1_harmonics_zero,E_1_harmonics_right];
        
        E_2_harmonics_right=harmonics_reconstruction_y(N_X,3,E_2_harmonics);    
        E_2_harmonics_left=-fliplr(E_2_harmonics_right);    
        E_2_harmonics_zero=zeros(2*N_X+1,1);    
        E_2_harmonics2=[E_2_harmonics_left,E_2_harmonics_zero,E_2_harmonics_right];
        
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

elseif symmetry_switch==4 % odd mode
    
    % (E^x)_(i,0)=0, (H^y)_(i,0)=0
    %( (E^y)_(i,0); (E^x)_(i,j); (E^y)_(i,j) )
    %( (H^x)_(i,0); (H^x)_(i,j); (H^y)_(i,j) )

    if plot_field_component==1 % E_1 
        % (E^x)_(i,j)=-(E^x)_(i,-j), (E^x)_(i,0)=0, (E^x)_(i,j) -> 2
        
        E_1_harmonics_right=harmonics_reconstruction_y(N_X,2,E_1_harmonics);    
        E_1_harmonics_left=-fliplr(E_1_harmonics_right);    
        E_1_harmonics_zero=zeros(2*N_X+1,1);
        E_1_harmonics2=[E_1_harmonics_left,E_1_harmonics_zero,E_1_harmonics_right];
    elseif plot_field_component==2 % E_2
        % (E^y)_(i,j)=(E^y)_(i,-j), (E^y)_(i,j) -> 3, (E^y)_(i,0) -> 1
        
        E_2_harmonics_right=harmonics_reconstruction_y(N_X,3,E_2_harmonics);    
        E_2_harmonics_left=fliplr(E_2_harmonics_right);    
        E_2_harmonics_zero=harmonics_reconstruction_y(N_X,1,E_2_harmonics);    
        E_2_harmonics2=[E_2_harmonics_left,E_2_harmonics_zero,E_2_harmonics_right];   

    elseif plot_field_component==3 % E_3, S_3=1i*inv(eps_33)*(K_x*U_2-K_y*U_1) 
        
        H_1_harmonics_right=harmonics_reconstruction_y(N_X,2,H_1_harmonics);    
        H_1_harmonics_left=fliplr(H_1_harmonics_right);    
        H_1_harmonics_zero=harmonics_reconstruction_y(N_X,1,H_1_harmonics);   
        H_1_harmonics2=[H_1_harmonics_left,H_1_harmonics_zero,H_1_harmonics_right];
        
        H_2_harmonics_right=harmonics_reconstruction_y(N_X,3,H_2_harmonics);    
        H_2_harmonics_left=-fliplr(H_2_harmonics_right);    
        H_2_harmonics_zero=zeros(2*N_X+1,1);
        H_2_harmonics2=[H_2_harmonics_left,H_2_harmonics_zero,H_2_harmonics_right];
        
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
        % (H^x)_(i,j)=(H^x)_(i,-j), (H^x)_(i,j) -> 2, (H_x)_(i,0) -> 1    
    
        H_1_harmonics_right=harmonics_reconstruction_y(N_X,2,H_1_harmonics);    
        H_1_harmonics_left=fliplr(H_1_harmonics_right);    
        H_1_harmonics_zero=harmonics_reconstruction_y(N_X,1,H_1_harmonics);   
        H_1_harmonics2=[H_1_harmonics_left,H_1_harmonics_zero,H_1_harmonics_right];    
    
    elseif plot_field_component==5 % H_2    
        % (H^y)_(i,j)=-(H^y)_(i,-j), (H_x)_(i,0)=0, (H^y)_(i,j) -> 3
        
        H_2_harmonics_right=harmonics_reconstruction_y(N_X,3,H_2_harmonics);    
        H_2_harmonics_left=-fliplr(H_2_harmonics_right);    
        H_2_harmonics_zero=zeros(2*N_X+1,1);
        H_2_harmonics2=[H_2_harmonics_left,H_2_harmonics_zero,H_2_harmonics_right];
        
    elseif plot_field_component==6 % H_3, U_3=1i*inv(mu_33)*(K_x*S_2-K_y*S_1)
        
        E_1_harmonics_right=harmonics_reconstruction_y(N_X,2,E_1_harmonics);    
        E_1_harmonics_left=-fliplr(E_1_harmonics_right);    
        E_1_harmonics_zero=zeros(2*N_X+1,1);
        E_1_harmonics2=[E_1_harmonics_left,E_1_harmonics_zero,E_1_harmonics_right];
        
        E_2_harmonics_right=harmonics_reconstruction_y(N_X,3,E_2_harmonics);    
        E_2_harmonics_left=fliplr(E_2_harmonics_right);    
        E_2_harmonics_zero=harmonics_reconstruction_y(N_X,1,E_2_harmonics);    
        E_2_harmonics2=[E_2_harmonics_left,E_2_harmonics_zero,E_2_harmonics_right];
        
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