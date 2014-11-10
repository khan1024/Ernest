% Author: Pavel Kwiecien, pavel.kwiecien@seznam.cz
% Czech Technical University in Prage, Optical Physics Group, Czech Republic

function [result_asr]=asr_fourier(f,disc_x,disc_u,M,d,G,use_grid_asr_speed_up,asr_function_only)

length_disc_x=length(disc_x);

% nastaveni x a x_m1
x_m1=zeros(1,length_disc_x-1);
x=zeros(1,length_disc_x-1);

% nastaveni u a u_m1
u_m1=zeros(1,length_disc_x-1);
u=zeros(1,length_disc_x-1);

%{
for s=2:1:length_disc_x
    x_m1(1,s-1)=disc_x(1,s-1);
    x(1,s-1)=disc_x(1,s);    
    u_m1(1,s-1)=disc_u(1,s-1);
    u(1,s-1)=disc_u(1,s);
end
%}
x_m1(1,1:length_disc_x-1)=disc_x(1,1:length_disc_x-1);
x(1,1:length_disc_x-1)=disc_x(1,2:length_disc_x);    
u_m1(1,1:length_disc_x-1)=disc_u(1,1:length_disc_x-1);
u(1,1:length_disc_x-1)=disc_u(1,2:length_disc_x);


% nastaveni asr parametru a_1, a_2, a_3
%a_1=zeros(1,length_disc_x-1);
a_2=zeros(1,length_disc_x-1);
a_3=zeros(1,length_disc_x-1);

%{
for s=1:1:length_disc_x-1
    a_1(1,s)=(u(1,s)*x_m1(1,s)-u_m1(1,s)*x(1,s))/(u(1,s)-u_m1(1,s));
    a_2(1,s)=(x(1,s)-x_m1(1,s))/(u(1,s)-u_m1(1,s));    
    a_3(1,s)=G*(u(1,s)-u_m1(1,s))-(x(1,s)-x_m1(1,s));
end
%}

%a_1(1,:)=(u(1,:)*x_m1(1,:)-u_m1(1,:)*x(1,s))/(u(1,:)-u_m1(1,:));
a_2(1,:)=(x(1,:)-x_m1(1,:))./(u(1,:)-u_m1(1,:));    
a_3(1,:)=G*(u(1,:)-u_m1(1,:))-(x(1,:)-x_m1(1,:));

result_asr=zeros(1,2*M+1);
result_asr_temp=zeros(length(x),2*M+1);

if use_grid_asr_speed_up==1
    %%{
        %modified_asr=zeros(length(x),2*M+1);
        modified_asr=conj((reshape(asr_function_only,[],length(x)))');
        %for i=1:1:length(x)
        %    modified_asr(i,:)=asr_function_only(1,:,i);
        %end
        
        [m,n]=size(f);
        if m<=n
            f=conj(f');
        end
        
        result_asr=sum(bsxfun(@times,f,modified_asr),1);
    %%}
        
        %{
        for i=1:1:length(x)
            result_asr=result_asr+f(i)*asr_function_only(1,:,i);
        end
        %}
else
    for i=1:1:length(x)    
        a_dif=u(i)+u_m1(i);    
        result_asr=result_asr+f(i)*(-1./2).*(-2.*u(i)+a_dif)./d./(u(i)-u_m1(i)).*(2.*exp(-1i.*pi.*(-M:M).*a_dif./d).*a_2(i).*sinc(-(-M:M).*(-2.*u(i)+...
            a_dif)./d).*u(i)-2.*exp(-1i.*pi.*(-M:M).*a_dif./d).*a_2(i).*sinc(-(-M:M).*(-2.*u(i)+a_dif)./d).*u_m1(i));
    end  
end