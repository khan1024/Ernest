N=number_of_orders;
d=Lambda;

data_11=data(l,2:end);

length_data=length(data_11);

data_eps11=0.*data_11;
data_inv_eps11=0.*data_11;


for i=1:1:length_data
    if mod(i,2)==1 % pozice
        data_eps11(i)=data_11(1,i);
        data_inv_eps11(i)=data_11(1,i);
        
    elseif mod(i,2)==0 % index
        data_inv_eps11(i)=1/data_11(1,i)^2;
        data_eps11(i)=data_11(1,i)^2;
        
    end
end

inv_eps11=function_toeplitz_matrix(data_inv_eps11,N,d);
eps11=function_toeplitz_matrix(data_eps11,N,d);

first_iteration=exist('Delta','var');
first_iteration=0;


if first_iteration==0
step=0.001;
x=0:step:Lambda;
stred_x=step/2:step:Lambda-step/2;

sirka_celek=thickness*1e6;

if grating==10
    jmenovatel_=sqrt(1+(pi*sirka_celek/Lambda)^2*(cos(K/1e6*stred_x)).^2);
    N_z_=-pi*sirka_celek/Lambda*cos(K/1e6*stred_x)./jmenovatel_; %přehozeno má být N_x_
    N_x_=1./jmenovatel_;
elseif mrizka==4
    %{
    prepona=sqrt(Lambda*Lambda/4+sirka_celek*sirka_celek);
    cos_psi=cos(Lambda/2/prepona);
    sin_psi=sin(sirka_celek/prepona);

    N_x_=0.*stred_x;
    N_z_=0.*stred_x;
    
    N_x_(:)=sin_psi;
    N_z_(1:ceil(length(stred_x)/2))=-cos_psi;
    N_z_(ceil(length(stred_x)/2)+1:end)=cos_psi;
    %}
end

input_Nx=zeros(1,length(x)+length(N_x_));
input_Nx(1:2:end)=x;
input_Nx(2:2:end-1)=N_x_;

input_Nz=zeros(1,length(x)+length(N_z_));
input_Nz(1:2:end)=x;
input_Nz(2:2:end-1)=N_z_;

%{
input_Nx2=zeros(1,length(x)+length(N_x_));
input_Nx2(1:2:end)=x;
input_Nx2(2:2:end-1)=N_x_.^2;

input_Nz2=zeros(1,length(x)+length(N_z_));
input_Nz2(1:2:end)=x;
input_Nz2(2:2:end-1)=N_z_.^2;

input_NxNz=zeros(1,length(x)+length(N_z_));
input_NxNz(1:2:end)=x;
input_NxNz(2:2:end-1)=N_x_.*N_z_;
%}

%toep_Nx2=function_toeplitz_matrix(input_Nx2,N,d);
%toep_Nz2=function_toeplitz_matrix(input_Nz2,N,d);
toep_Nx=function_toeplitz_matrix(input_Nx,N,d);
toep_Nz=function_toeplitz_matrix(input_Nz,N,d);
%toep_NxNz=function_toeplitz_matrix(input_NxNz,N,d);
end

norm_eps=eps11;
inve_eps=inv(inv_eps11);

Delta=(inve_eps-norm_eps);

%{
Q_11=norm_eps*toep_Nz2+inve_eps*toep_Nx2;
Q_12=-Delta*toep_NxNz;
Q_13=zero;
Q_21=-Delta*toep_NxNz;
Q_22=norm_eps*toep_Nx2+inve_eps*toep_Nz2;
Q_23=zero;
Q_31=zero;
Q_32=zero;
Q_33=norm_eps;
%}

%%{
Q_11=toep_Nz*norm_eps*toep_Nz+toep_Nx*inve_eps*toep_Nx;
Q_12=toep_Nx*norm_eps*toep_Nz-toep_Nz*inve_eps*toep_Nx;
Q_13=zero;
Q_21=toep_Nz*norm_eps*toep_Nx-toep_Nx*inve_eps*toep_Nz;
Q_22=toep_Nx*norm_eps*toep_Nx+toep_Nz*inve_eps*toep_Nz;
Q_23=zero;
Q_31=zero;
Q_32=zero;
Q_33=norm_eps;

inv_Q_22=inv(Q_22);

gamma=0;
mu=1;

alpha=K_x.*k_0/1e6;

A_22=-Q_12*inv_Q_22*alpha;    
A_24=gamma^2/((k_0/1e6)*mu)*I-(k_0/1e6)*(Q_11-Q_12*inv_Q_22*Q_21);

A_42=1/(k_0/1e6)*(alpha*inv_Q_22*alpha-mu*(k_0/1e6)^2*I);
A_44=-alpha*inv_Q_22*Q_21;

% H3, E1
[eigvec_TM,eigval_TM]=eig([A_22,A_24;A_42,A_44]);

eigval_TM=eigval_TM*1i;

eigval_TM_plus=[];
eigval_TM_minus=[];    

eigvec_TM_plus=[];%zeros(2*N,1);
eigvec_TM_minus=[];%zeros(2*N,1);

    
for i=1:1:2*N        
    if real(eigval_TM(i,i))+imag(eigval_TM(i,i))>0            
        if real(eigval_TM(i,i))<0                
            eigval_TM_plus=[eigval_TM_plus;eigval_TM(i,i)];% nespravna faktorizace -eigval_TM(i,i)];                
            eigvec_TM_plus=[eigvec_TM_plus,eigvec_TM(:,i)];%-eigvec(:,i)];            
        else            
            eigval_TM_plus=[eigval_TM_plus;eigval_TM(i,i)];                
            eigvec_TM_plus=[eigvec_TM_plus,eigvec_TM(:,i)];            
        end        
        
    elseif real(eigval_TM(i,i))+imag(eigval_TM(i,i))<0            
        if real(eigval_TM(i,i))>0                
            eigval_TM_minus=[eigval_TM_minus;eigval_TM(i,i)];%-eigval_TM(i,i)];
            eigvec_TM_minus=[eigvec_TM_minus,eigvec_TM(:,i)];%-eigvec(:,i)];            
        else            
            eigval_TM_minus=[eigval_TM_minus;eigval_TM(i,i)];                
            eigvec_TM_minus=[eigvec_TM_minus,eigvec_TM(:,i)];            
        end        
    end    
end

Q(:,:,l)=diag((-eigval_TM_minus)/(k_0/1e6));
W(:,:,l)=eigvec_TM_minus(1:N,1:N);
V(:,:,l)=1i*eigvec_TM_minus(N+1:end,1:N);